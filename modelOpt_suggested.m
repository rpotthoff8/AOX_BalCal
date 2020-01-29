function [customMatrix_rec, FLAGS]= modelOpt_suggested(VIFthresh, customMatrix, customMatrix_req, loaddimFlag, nterms, comIN0, anova_pct, targetMatrix0, high, FLAGS)
% Function searches for the 'suggested math model' using the approach
% described in Balfit Reference B29.  This approach is a simplified and
% faster search process than searching for the "recommended" math model.
% It does not optimize as search metric, but instead just searches for an
% equation where there are no near-linear dependencies between math terms,
% the hierarchy rule is satisfied (optional), and all terms are
% statistically significant.  It searches between the upper and lower
% bounds for terms in the model:
% Upper Bounds: provided current customMatrix
% Lower Bounds: Linear voltage from channel and tares

%The following process is used in this function:
% 1) Start with the provided 'permitted math model', obtained by applying
%       SVD
% 2) PHASE 1: Apply VIF constraint: Compute VIF for all terms, remove term with
%       highest VIF, recompute VIF, repeat until max VIF is below threshold
% 3) PHASE 2: Apply Hierarchy Rule (optional): Add missing lower order terms
% 4) PHASE 3: Remove statistically 'insignificant' math terms: Find P values of all
%       regression coefficients from terms from step 3.  Remove all insignificant
%       terms as long as math model remains "hierarchical" (remove high order to
%       low order)

%INPUTS:
%  VIFthresh = Threshold for max allowed VIF (Balfit Search constraint 2)
%  customMatrix = Current matrix of what terms should be included in Eqn Set.
%  customMatrix_req = custom matrix of minimum terms that must be included in model
%  loaddimFlag = Dimension of load data (# channels)
%  nterms = Number of predictor terms in regression model
%  comIN0 = Matrix of predictor variables
%  anova_pct = ANOVA percent confidence level for determining significance
%  targetMatrix0 = Matrix of target load values
%  high = Matrix of term hierarchy
%  FLAGS = Structure containing flags for user preferences

%OUTPUTS:
%  customMatrix_rec = Optimized recommended custom matrix
%  FLAGS = Structure containing flags for user preferences

fprintf('\nCalculating Suggested Eqn Set with legacy efficient search method....')

high_con=FLAGS.high_con;
opt_channel=FLAGS.opt_channel; %Flag for if each channel should be optimized

% Normalize the data for a better conditioned matrix (Copy of what is done
% in calc_xcalib.m)
scale = max(abs(comIN0));
scale(scale==0)=1; %To avoid NaN for channels where RBFs have self-terminated
comIN0 = comIN0./scale;

%% Phase 1: Sequentially remove terms with highest VIF until max VIF is <
% threshold. Sequence will never remove tare intercept terms or
%linear voltage from corresponding channel
fprintf('\n     Starting Phase 1 of optimization.... \n')

[customMatrix_P1, P1_maxVIF_hist, P1_maxVIF_init, opt_channel] = VIF_constraint(VIFthresh, customMatrix, loaddimFlag, nterms, comIN0, opt_channel, high_con, high, customMatrix_req, FLAGS); %Remove terms based on VIF constraint

%% Phase 2: Apply 'Hierarchy Rule': Add missing lower order terms to customMatrix from Phase 1
fprintf('\n     Starting Phase 2 of optimization.... \n')
customMatrix_P2=customMatrix_P1; %Initialize as phase 1 result

if high_con==1 %If enforcing hierarchy constraint
    
    identical_custom=all(~diff(customMatrix_P2,1,2),'all'); %Check if customMatrix is identical for all channels
    if identical_custom==1 %If custom equation identical for each channel
        calcThrough=1; %only necessary to remove enforce hierarchy rule once
    else
        calcThrough=loaddimFlag; %Necessary to calculate for each channel seperate
    end
    
    for i=1:calcThrough %Loop through channels
        if opt_channel(i)==1 %Check if channel model should be optimized
            %         inc_idx=find(customMatrix_P2(:,i)); %Find index of included terms for channel
            sup_terms_mat=high(boolean(customMatrix_P2(1:nterms,i)),:); %Rows from hierarchy matrix for included terms. columns with '1' are needed to support variable
            sup_terms=any(sup_terms_mat,1); %Row vector with 1s for terms needed to support included terms
            customMatrix_P2(boolean(sup_terms),i)=1; %Include all terms needed to support currently included terms
        end
    end
    
    if calcThrough==1 %If final Phase 2 equation identical for each channel
        customMatrix_P2(:,2:loaddimFlag)=repmat(customMatrix_P2(:,1),1,loaddimFlag-1); %Duplicate for each column
    end
end



%% Phase 3: Remove all insignificant terms (not required model terms).
%Phase will never remove required term from a channel. If tare intercept terms are not significant,
%warning will be provided. Never remove terms if they are needed to support
%the model (if hierarchy selected)

%Unlike Balfit, after removing all 'insignificant' terms, significance is
%recalculated to ensure all remaining terms are now significant.  Process
%is repeated until all included terms are significant. (Repeat process is
%ONLY performed in hierarchy rule is not being enforced since by nature
%some insignificicant support terms will be included if hierarchy is turned
%on)
fprintf('\n     Starting Phase 3 of optimization.... \n')

customMatrix_P3=customMatrix_P2; %Initialize as result from phase 2

sig=zeros(size(customMatrix_P3)); %initialize matrix for significance
p_val=zeros(size(customMatrix_P3)); %initialize matrix for p-values
for i=1:loaddimFlag
    if opt_channel(i)==1 %Check if channel model should be optimized
        
        %Find significance of terms in phase 3 customMatrix
        [sig(boolean(customMatrix_P3(:,i)),i),p_val(boolean(customMatrix_P3(:,i)),i)]=sig_check(comIN0(:,boolean(customMatrix_P3(:,i))),targetMatrix0(:,i),anova_pct); %Calculate significance for all initial included terms
        
        if any(~boolean(sig(boolean(customMatrix_P3(1:nterms,i)),i))) %Check if all included terms are not significant
            sig_all=0; %Not all included terms are significant, need to remove insignificant terms
        else
            sig_all=1; %All all included terms are significant
        end
        
        while sig_all==0 %Repeat until all terms are significant (if not enforcing hierarchy)
            sig_all=1; %CHANGE, for testing
            customMatrix_iter=zeros(size(customMatrix_P3,1),1); %Initialize as zeros
            
            customMatrix_iter(boolean(sig(:,i)))=1; %Include significant terms
            customMatrix_iter(boolean(customMatrix_req(:,i)))=1; %Include terms from minimum equation
            
            customMatrix_P3(:,i)=customMatrix_iter; %Set new custom matrix for channel
            
            if high_con==1 %If enforcing hierarchy constraint: Include terms needed to support significant terms
                sup_terms_sig_mat=high(boolean(sig(1:nterms,i)),:); %Rows from hierarchy matrix for included significant terms. columns with '1' are needed to support variable
                sup_terms_sig=any(sup_terms_sig_mat,1); %Row vector with 1s for terms needed to support significant included terms
                customMatrix_P3(boolean(sup_terms_sig),i)=1; %Include all terms needed to support significant terms
                sig_all=1; %Flip flag, do not repeate loop for significance since including hierarchy terms that may not be significant
            end
            
            %check for term significance with new model
            sig(:,i)=0; %Reset
            p_val(:,i)=0; %Reset
            %Find significance of terms in new phase 3 customMatrix
            [sig(boolean(customMatrix_P3(:,i)),i),p_val(boolean(customMatrix_P3(:,i)),i)]=sig_check(comIN0(:,boolean(customMatrix_P3(:,i))),targetMatrix0(:,i),anova_pct); %Calculate significance for all initial included terms
            
            if all(boolean(sig(boolean(customMatrix_P3(1:nterms,i)),i))) %Check if all included terms are significant
                sig_all=1; %All included terms are significant, exit loop
            elseif all(boolean(sig(boolean(customMatrix_P3(setdiff(1:nterms,customMatrix_req(1:nterms,i)),i)),i))) %Check if all included terms are significant except required terms from channel
                fprintf('       Possible insignificant term from required equation in channel: '); fprintf(num2str(i)); fprintf('.\n'); %Display error message
                sig_all=1; %Exit loop
            end
        end
        
        if any(~boolean(sig(boolean(customMatrix_req(nterms+1:end,i)),i))) %If any tare intercept terms are insignificant
            fprintf('       Possible insignificant tare intercept in channel: '); fprintf(num2str(i)); fprintf('.\n'); %Display error message
        end
    end
end

customMatrix_rec=customMatrix; %Initialize recommended custom matrix as provided custom matrix
customMatrix_rec(:,boolean(opt_channel))=customMatrix_P3(:,boolean(opt_channel)); %Set optimized channels to P4 Results

FLAGS.opt_channel=opt_channel; %Output tracker of what channels were able to calculate suggested equation
fprintf('Suggested Equation Search Complete. \n ')
end

function [customMatrix_Px, Px_maxVIF_hist, Px_maxVIF_init, optChannel] = VIF_constraint(VIFthresh, customMatrix, loaddimFlag, nterms, comIN0, optChannel, high_con, high, min_eqn, FLAGS)
%% Phase 1: Sequentially remove terms with highest VIF until max VIF is <
%threshold. Sequence will never remove tare intercept terms or
%linear voltage from corresponding channel
customMatrix_Px =customMatrix; %Initialize as current custom matrix.  Will trim down terms from there

identical_custom=all(~diff(customMatrix,1,2),'all'); %Check if customMatrix is identical for all channels
if identical_custom==1 %If custom equation identical for each channel
    calcThrough=1; %only necessary to remove VIF terms for 1 channel in Phase 1
else
    calcThrough=loaddimFlag; %Necessary to calculate for each channel seperate
end

Px_maxVIF_hist=zeros(nterms,loaddimFlag); %variable for tracking history of max VIF
Px_maxVIF_init=zeros(1,loaddimFlag); %vector for initial max VIF
i=1; %Initialize counter variable
while i<=calcThrough %Loop through, calculating for each load channel if necessary
    if optChannel(i)==1 %Check if channel model should be optimized
        
        %Initial check on VIF
        vif_iter=zeros(size(customMatrix,1),1); %Initialize VIF vector
        vif_iter(boolean([0;customMatrix_Px(2:end,i)]))=vif_dl(comIN0(:,boolean([0;customMatrix_Px(2:end,i)]))); %Calculate VIF for all initial included terms except global intercept
        vif_iter(1)=1; %Set VIF for intercept=1

        Px_maxVIF_init(i)=max(vif_iter); %Store initial max VIF
        if Px_maxVIF_init(i)>VIFthresh %If initial max VIF is over elevated threshold
            vifHigh=1; %Flag for if max VIF is too high
        else
            vifHigh=0; %Elevated VIF threshold already satisfied
        end
        
        %Loop to remove terms until VIF threshold is met
        termOut_count=1;
        while vifHigh==1 %Remove terms until VIF threshold is met
            vif_iter_search=vif_iter; %initialize variable for searching for max VIF
            vif_iter_search(boolean(min_eqn(:,i)))=0; %Zero out corresponding minimum eqn so it will not be selected for removal
            [maxVIF_iter,maxI]=max(vif_iter_search); %Find index of max VIF out of allowable selections (not tare intercepts or linear voltages from corresponding channel)
            if maxVIF_iter<=VIFthresh %Check if all 'problem' terms have been removed
%                 fprintf('       Check for linear dependence in minimum required equation for Channel: '); fprintf(num2str(i)); fprintf('.\n');
%                 break
            elseif maxVIF_iter== 0 
                fprintf('       Error occured in Phase 1.');
                fprintf(' Check for severe linear dependence in minimum required equation for Channel: '); fprintf(num2str(i)); fprintf('.\n');
                optChannel(i)=0; %Mark unable to optimize for channel
                break
            end
            customMatrix_Px(maxI,i)=0; %Eliminate term with highest VIF
            %             if high_con==1 %If enforcing hierarchy constraint
            %                 customMatrix_Px(boolean(high(:,maxI)),i)=0; %Removes terms that required variable for support
            %             end
            
            vif_iter=zeros(size(customMatrix,1),1); %Reset VIF vector
            vif_iter(boolean([0;customMatrix_Px(2:end,i)]))=vif_dl(comIN0(:,boolean([0;customMatrix_Px(2:end,i)]))); %Calculate VIF for all new included terms except global intercept
            vif_iter(1)=1; %Set VIF for global intercept = 1
            Px_maxVIF_hist(termOut_count,i)=max(vif_iter); %Store new max VIF
            
            if calcThrough~=loaddimFlag && any(min_eqn(maxI,:)) %If calculating not full number of channels (for time savings) but eliminated a linear voltage term
                calcThrough=loaddimFlag; %Now necessary to calculate for each channel seperate
            end
            
            if Px_maxVIF_hist(termOut_count,i)<=VIFthresh %If new max VIF is <= elevated threshold
                vifHigh=0; %Flip flag. Constraint has been met
            else
                termOut_count=termOut_count+1; %Advance counter
            end
        end
    end
    i=i+1; %Advance channel counter
end
if calcThrough==1 %If final Phase 1 equation identical for each channel
    customMatrix_Px(:,2:loaddimFlag)=repmat(customMatrix_Px(:,1),1,loaddimFlag-1); %Duplicate for each column
end

end
