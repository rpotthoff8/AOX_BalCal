function [customMatrix_rec, FLAGS]= modelOpt_backward(VIFthresh, customMatrix_permit, customMatrix_req, loaddimFlag, nterms, comIN0, anova_pct, targetMatrix0, high, FLAGS)
% Function searches for the 'recommended equation' using the approach
% Balfit reference B16 describes as "Backward Elimination".
% Function determines optimal combination of terms to include between
% possible upper and lower bounds.
% Upper Bounds: provided current permitted customMatrix
% Lower Bounds: Linear voltage from channel and tares

%INPUTS:
%  VIFthresh = Threshold for max allowed VIF (Balfit Search constraint 2)
%  customMatrix = Current matrix of what terms should be included in Eqn Set.
%  customMatrix_req = Minimum terms that must be included in Eqn Set.
%  loaddimFlag = Dimension of load data (# channels)
%  nterms = Number of predictor terms in regression model
%  comIN0 = Matrix of predictor variables
%  anova_pct = ANOVA percent confidence level for determining significance
%  targetMatrix0 = Matrix of target load values
%  high = Matrix of term hierarchy
%  FLAGS.high_con = Flag for if term hierarchy constraint is enforced (0=
%       Off, 1= Enforced after search (model is optimized, terms required to
%       support final model are added in), 2=Enforced during search (model is
%       optimized, terms are only added if they are supported by existing model)
%  FLAGS.search_metric_flag = Flag for search metric to return: 1=PRESS (minimize),
%       2=Square Root of Residual Mean Square (minimize), 3=F-Value (maximize)
%  FLAGS.VIF_stop_flag = Flag to terminate search once VIF threshold is
%       exceeded

%OUTPUTS:
%  customMatrix_rec = Optimized recommended custom matrix
% FLAGS = Structure containing flags for user preferences

high_con=FLAGS.high_con;
quick_search_flag=0;
opt_channel=FLAGS.opt_channel; %Output tracker of what channels were able to calculate permitted equation

FLAGS.VIF_stop=0;
VIFthresh_elev=max([50,VIFthresh]);

fprintf('\nCalculating Recommended Eqn Set with Backward Elimination Method....')

% Normalize the data for a better conditioned matrix (Copy of what is done
% in calc_xcalib.m)
scale = max(abs(comIN0));
scale(scale==0)=1; %To avoid NaN for channels where RBFs have self-terminated
comIN0 = comIN0./scale;

%Define number of series included:
nseries=size(customMatrix_permit,1)-nterms;

%Initialize optimized math model as permitted (upper end)
customMatrix_opt=customMatrix_permit;

num_permit=sum(customMatrix_permit); %Max number of terms
num_req=sum(customMatrix_req); %Min number of terms
num_terms=sum(customMatrix_opt); %Current count of terms
num_test=(num_permit-num_req)+1; %Number of models to find. 1 for every number of terms from current # to max #

%Initialize variables
VIF_met=false(max(num_test),loaddimFlag); %Matrix for storing if all terms meet VIF constraint
VIF_elev_met=false(max(num_test),loaddimFlag); %Matrix for storing if all terms meet elevated VIF constraint
VIF_max=zeros(max(num_test),loaddimFlag); %Initialize variable for storing max VIF at each # terms
sig_all=false(max(num_test),loaddimFlag); %Matrix for storing if all terms are sig at each # terms
P_max=zeros(max(num_test),loaddimFlag); %Matrix for storing max coefficeint p_value at each # terms
search_metric=zeros(max(num_test), loaddimFlag); %Matrix for storing seach metric value at each # terms

% if VIF_stop_flag==1
%     VIF_blacklist=zeros(size(customMatrix_opt)); %Matrix for tracking terms that violate VIF limit
% end

if high_con==1
    high_sup=false(max(num_test),loaddimFlag); %Initialize variable for tracking if model is hierarchically supported
end

for i=1:loaddimFlag %Loop through all channels
    if opt_channel(i)==1 %If optimization is turned on
        %Check initial (required) math model
        VIF_vect=zeros(size(customMatrix_opt,1),1);
        sig_vect=zeros(size(customMatrix_opt,1),1);
        P_vect=zeros(size(customMatrix_opt,1),1);
        
        [VIF_met(1,i),VIF_max(1,i),sig_all(1,i),P_max(1,i),search_metric(1,i),VIF_vect(boolean(customMatrix_opt(:,i))), sig_vect(boolean(customMatrix_opt(:,i))),P_vect(boolean(customMatrix_opt(:,i)))]=test_combo(comIN0(:,boolean(customMatrix_opt(:,i))), targetMatrix0(:,i), anova_pct, VIFthresh, nseries, FLAGS);
        VIF_elev_met(1,i)=VIF_max(1,i)<=VIFthresh_elev; %Check if model meets elevated VIF threshold
        
        if high_con==1
            sup_terms_mat=high(boolean(customMatrix_opt(1:nterms,i)),:); %Rows from hierarchy matrix for included terms. columns with '1' are needed to support variable
            sup_terms=any(sup_terms_mat,1); %Row vector with 1s for terms needed to support included terms
            high_sup(1,i)=all(customMatrix_opt(boolean(sup_terms),i));
        end
        
        customMatrix_hist=zeros(size(customMatrix_opt,1),num_test(i)); %Matrix for storing custom matrix used
        customMatrix_hist(:,1)=customMatrix_opt(:,1); %First model is required model
        for j=2:num_test(i) %Loop through each # terms from current number to max number
            
            if VIF_elev_met(j-1,i)==0
                cur_constraint=1; %Current constraint is elevated VIF
            elseif sig_all(j-1,i)==0
                cur_constraint=2; %Current constraint is term significance
            elseif VIF_met(j-1,i)==0
                cur_constraint=3; %Current constraint is final VIF
            else
                cur_constraint=0; %Current model meets both constraints
            end
            %             if VIF_elev_met(j-1,i)==0 %Elevated VIF constraint has priority
            %                 elim_terms=VIF_vect>VIFthresh_elev;
            %             elseif sig_all(j-1,i)==0 %Significance constraint has priority
            %                 elim_terms=~sig_vect;
            %             elseif VIF_met(j-1,i)==0 %Final VIF constraint has priority
            %                 elim_terms=VIF_vect>VIFthresh;
            %             else
            %                 elim_terms=true(size(VIF_vect));
            %             end
            elim_terms=true(size(VIF_vect));
            
            
            %Possible terms to be subtracted are those in current model that are not in required model
            pos_sub=all([customMatrix_opt(:,i), ~customMatrix_req(:,i),elim_terms],2);
            
            %             pos_sub=zeros(size(customMatrix_opt,1),1); %Initialize as zeros
            %             pos_sub(~boolean(customMatrix_req(:,i)))=customMatrix_opt((~boolean(customMatrix_req(:,i)) && boolean(elim_terms)),i); %Vector of terms that can be subtracted
            
            %             if VIF_stop_flag==1 %If terminating search based on VIF limit
            %                 pos_sub(boolean(VIF_blacklist(:,i)))=0; %Don't add terms that are on 'blacklist' for exceeding VIF threshold
            %             end
            
            
            if high_con==2 %If enforcing hierarchy constraint during search
                %Terms are only possible for subtraction if they are not
                %supporting any other term
                sup_Terms_mat=high(boolean(customMatrix_opt(1:nterms,i)),:); %Matrix of terms that are needed to support current terms
                %                 sup_diff=sup_Terms-customMatrix_opt(1:nterms,i)';
                sup_Terms_vec=any(sup_Terms_mat,1); %Find terms that are supporting any other terms
                pos_sub(boolean(sup_Terms_vec))=0; %Remove terms that are supporting from possibilities to subtract
            end
            pos_sub_idx=find(pos_sub); %Index of possible terms for subtracting from model
            
            if isempty(pos_sub_idx) %If no possible terms to subtract
                break; %Exit for loop testing decreasing # of math models
            end
            
            num_pos_sub=numel(pos_sub_idx); %Count of terms to be tested for subtracting from model
            
            %Initialize
            VIF_met_temp=false(num_pos_sub,1); %Matrix for storing if all terms meet VIF constraint
            VIF_elev_met_temp=false(num_pos_sub,1); %Matrix for storing if all terms meet elevated VIF constraint
            VIF_max_temp=zeros(num_pos_sub,1); %Initialize variable for storing max VIF at each # terms
            sig_all_temp=false(num_pos_sub,1); %Matrix for storing if all terms are sig at each # terms
            P_max_temp=zeros(num_pos_sub,1); %Matrix for storing max coefficeint p_value at each # terms
            search_metric_temp=zeros(num_pos_sub, 1); %Matrix for storing seach metric value at each # terms
            VIF_vect_temp=zeros(size(customMatrix_opt,1),num_pos_sub);
            sig_vect_temp=zeros(size(customMatrix_opt,1),num_pos_sub);
            P_vect_temp=zeros(size(customMatrix_opt,1),num_pos_sub);
            Press_temp=zeros(num_pos_sub, 1); %Matrix for storing seach metric value at each # terms
            F_temp=zeros(num_pos_sub, 1); %Matrix for storing seach metric value at each # terms
            if high_con==1
                high_sup_temp=false(num_pos_sub,1);
            end
            
            for k=1:num_pos_sub %Loop through, testing each possible term to add
                customMatrix_opt_temp=customMatrix_opt(:,i); %Initialize as current custom matrix
                customMatrix_opt_temp(pos_sub_idx(k))=0; %Subtract term for test
                %Test math model with new term subtracted
                [VIF_met_temp(k),VIF_max_temp(k),sig_all_temp(k),P_max_temp(k),search_metric_temp(k),VIF_vect_temp(boolean(customMatrix_opt_temp),k),sig_vect_temp(boolean(customMatrix_opt_temp),k),P_vect_temp(boolean(customMatrix_opt_temp),k), Press_temp(k), F_temp(k)]=test_combo(comIN0(:,boolean(customMatrix_opt_temp)), targetMatrix0(:,i), anova_pct, VIFthresh, nseries, FLAGS);
                VIF_elev_met_temp(k)=VIF_max_temp(k)<=VIFthresh_elev; %Check if model meets elevated VIF threshold
                
                %                 if VIF_stop_flag==1 && VIF_met_temp(k)==0 %If adding term violates VIF limit
                %                     VIF_blacklist(pos_sub_idx(k),i)=1; %Add term to blacklist.  Will not try to add again
                %                 end
                if high_con==1 %Check if model is supported
                    sup_terms_mat=high(boolean(customMatrix_opt_temp(1:nterms,i)),:); %Rows from hierarchy matrix for included terms. columns with '1' are needed to support variable
                    sup_terms=any(sup_terms_mat,1); %Row vector with 1s for terms needed to support included terms
                    high_sup_temp(k)=all(customMatrix_opt_temp(boolean(sup_terms),i));
                end
                
            end
            
            
            cMeet_Idx=find(all([VIF_met_temp,sig_all_temp],2)); %Index of tests that met both VIF and significance constraint tests
            c2Meet_Idx=find(all([sig_all_temp,VIF_elev_met_temp],2)); %Index of tests that met both significance constraint test and elevated VIF test
            c3Meet_Idx=find(VIF_elev_met_temp); %Index of terms that meet elevated VIF test
            if ~isempty(cMeet_Idx) %If any tests met both constraints
                [~,k_best]=min(search_metric_temp(cMeet_Idx)); %Index of best search metric out of those meeting constraints
                Idx_best=cMeet_Idx(k_best); %Index of best search metric for all possible
                
                cur_constraint=0;
                
            elseif ~isempty(c2Meet_Idx) %If any tests met significance constraint and elevated VIF constraint
                if high_con==1
                    c2Meet_Idx(~boolean(high_sup_temp(c2Meet_Idx)))=[]; %Eliminate from consideration non-supported models
                end
                
%                 c2Meet_Idx(boolean(VIF_vect(pos_sub_idx(c2Meet_Idx))<=VIFthresh))=[]; %Eliminate from consideration non VIF violating terms
                
                if cur_constraint==2 %If switching from significance to final VIF
                    [~,k_best]=min(VIF_max_temp(c2Meet_Idx)); %Minimize max VIF for term removed
                    Idx_best=c2Meet_Idx(k_best); %Index of min VIF for all possible
                else
                    %Find index of terms that have max VIF after eliminating terms
                    [~,VIFmax_Idx]=max(VIF_vect_temp);
                    VIFmax_Idx=unique(VIFmax_Idx);
                    
                    c2Meet_Idx=c2Meet_Idx(ismember(pos_sub_idx(c2Meet_Idx),VIFmax_Idx));
                    
                    [~,k_best]=min(search_metric_temp(c2Meet_Idx)); %Index of best search metric out of those meeting constraints
                    Idx_best=c2Meet_Idx(k_best); %Index of best search metric for all possible
                end
                
                cur_constraint=3;
                
            elseif ~isempty(c3Meet_Idx) %If any test met elevated VIF constraint
                if high_con==1
                    c3Meet_Idx(~boolean(high_sup_temp(c3Meet_Idx)))=[]; %Eliminate from consideration non-supported models
                end
                
%                 c3Meet_Idx(boolean(sig_vect(pos_sub_idx(c3Meet_Idx))))=[]; %Eliminate from consideration significant terms
                
                if cur_constraint==1 %If switching from elevated VIF to significance constraint
                    [~,k_best]=min(P_max_temp(c3Meet_Idx)); %Minimize max P for term removed
                    Idx_best=c3Meet_Idx(k_best); %Index of min P for all possible
                else
                    %Find index of terms that have max P after eliminating terms
                    [~,Pmax_Idx]=max(P_vect_temp);
                    Pmax_Idx=unique(Pmax_Idx);
                    
                    c3Meet_Idx=c3Meet_Idx(ismember(pos_sub_idx(c3Meet_Idx),Pmax_Idx));
                    
                    [~,k_best]=min(search_metric_temp(c3Meet_Idx)); %Index of best search metric out of those meeting constraints
                    Idx_best=c3Meet_Idx(k_best); %Index of best search metric for all possible
                end
                
                cur_constraint=2;
            else %If no tests meet any constraints
                %Find index of terms that have max VIF after eliminating terms
                [~,VIFmax_Idx]=max(VIF_vect_temp);
                VIFmax_Idx=unique(VIFmax_Idx);
                %Find index of test elimination terms that correspond to these terms
                probVIF=find(ismember(pos_sub_idx,VIFmax_Idx));
                
                if high_con==1
                    probVIF(~boolean(high_sup_temp(probVIF)))=[]; %Eliminate from consideration non-supported models
                end
                
                %Find minimum of search metric in test removal terms that
                %are problem terms
                [~,k_best]=min(search_metric_temp(probVIF)); %Just find minimum of search metric
                Idx_best=probVIF(k_best);
                
                cur_constraint=1;
            end
            
            %Pick best term to add to move forward
            customMatrix_opt(pos_sub_idx(Idx_best),i)=0; %Remove term from customMatrix
            
            if high_con==1 %Check if model is supported
                sup_terms_mat=high(boolean(customMatrix_opt(1:nterms,i)),:); %Rows from hierarchy matrix for included terms. columns with '1' are needed to support variable
                sup_terms=any(sup_terms_mat,1); %Row vector with 1s for terms needed to support included terms
                high_sup(j,i)=all(customMatrix_opt(boolean(sup_terms),i));
            end
            
            %Store results
            customMatrix_hist(:,j)=customMatrix_opt(:,i); %Store custom matrix
            VIF_met(j,i)=VIF_met_temp(Idx_best);
            VIF_elev_met(j,i)=VIF_elev_met_temp(Idx_best);
            VIF_max(j,i)=VIF_max_temp(Idx_best);
            sig_all(j,i)=sig_all_temp(Idx_best);
            P_max(j,i)=P_max_temp(Idx_best);
            search_metric(j,i)=search_metric_temp(Idx_best);
            VIF_vect=VIF_vect_temp(:,Idx_best);
            sig_vect=sig_vect_temp(:,Idx_best);
            P_vect=P_vect_temp(:,Idx_best);
        end
        %Now select math model that minimizes search metric and meets both
        %constraints:
        if high_con==1 %If enforcing hierarchy constraint after, only consider models that are also supported hierarchially
            cMeet_Idx=find(all([VIF_met(:,i),sig_all(:,i),high_sup(:,i)],2)); %Index of tests that met both VIF and significance constraint tests
        else
            cMeet_Idx=find(all([VIF_met(:,i),sig_all(:,i)],2)); %Index of tests that met both VIF and significance constraint tests
        end
        
        if ~isempty(cMeet_Idx) %If any tests met both constraints
            [~,k_best]=min(search_metric(cMeet_Idx,i)); %Index of best search metric out of those meeting constraints
            Idx_best=cMeet_Idx(k_best); %Index of best search metric for all possible
            %Pick best term to add to move forward
            customMatrix_opt(:,i)=customMatrix_hist(:,Idx_best); %Add term to customMatrix
            
            %             if high_con==1 %If enforcing hierarchy constraint after, add in terms needed to support model
            %                 sup_terms_mat=high(boolean(customMatrix_opt(1:nterms,i)),:); %Rows from hierarchy matrix for included terms. columns with '1' are needed to support variable
            %                 sup_terms=any(sup_terms_mat,1); %Row vector with 1s for terms needed to support included terms
            %                 customMatrix_opt(boolean(sup_terms),i)=1; %Include all terms needed to support currently included terms
            %             end
        else %No math models met both constraints: Return error message and do not optimize channel
            fprintf('\nERROR: Unable to find math model that meets constraints for channel '); fprintf(num2str(i)); fprintf('.\n');
            opt_channel(i)=0;
        end
        
    end
end
%Output final model:
customMatrix_rec=customMatrix_permit; %Initialize recommended custom matrix as provided custom matrix
customMatrix_rec(:,boolean(opt_channel))=customMatrix_opt(:,boolean(opt_channel)); %Set optimized channels to optimal Results

FLAGS.opt_channel=opt_channel;
fprintf('\nRecommended Equation Search Complete. \n ')
end

