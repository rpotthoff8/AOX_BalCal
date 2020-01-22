function [customMatrix_rec]= modelOpt_backward(VIFthresh, customMatrix, loaddimFlag, nterms, comIN0, anova_pct, targetMatrix0, high, high_con)
% Function searches for the 'recommended equation' using the approach
% Balfit reference B16 describes as "Backward Elimination".
% Function determines optimal combination of terms to include between
% possible upper and lower bounds.
% Upper Bounds: provided current customMatrix
% Lower Bounds: Linear voltage from channel and tares

%See reference for more complete explanation.  Basic explanation of process from B16, pg 5:
% The alternate search strategy, i.e., backward elimination, recognizes the benefit of math term combinations
% during every step of the search as the regression model is built by removing terms during the search.
% Backward elimination, however, cannot enforce the threshold based search constraints (constraint 2 & 3 in
% Fig. 3) simultaneously as a term selected for removal may not necessarily violate all constraints. Therefore,
% it was decided to implement the search constraints in “series” if backward elimination is selected as the
% search strategy. This implementation works in three steps. First, a regression model is tested using search
% constraint 2 (near–linear dependency test) with a temporarily elevated threshold of 50. This test has priority
% over search constraint 3 (test of the significance of a term) until search constraint 2 is no longer violated
% during the search process. Then, search constraint 3 takes over. Now, search constraint 3 has priority over
% search constraint 2. Finally, search constraint 2 takes over again using the thresholds of 5 (or 10) after search
% constraint 3 can no longer be violated during the search process.

%INPUTS:
%  VIFthresh = Threshold for max allowed VIF (Balfit Search constraint 2)
%  customMatrix = Current matrix of what terms should be included in Eqn Set.
%  loaddimFlag = Dimension of load data (# channels)
%  nterms = Number of predictor terms in regression model
%  comIN0 = Matrix of predictor variables
%  anova_pct = ANOVA percent confidence level for determining significance
%  targetMatrix0 = Matrix of target load values
%  high = Matrix of term hierarchy
%  high_con = Flag for if term hierarchy constraint is enforced

%OUTPUTS:
%  customMatrix_rec = Optimized recommended custom matrix

fprintf('\nCalculating Recommended Eqn Set with Backwards Elimination Method....')

optChannel=ones(1,loaddimFlag); %Flag for if each channel should be optimized

%% Phase 1: Sequentially remove terms with highest VIF until max VIF is <
%elevated threshold. Sequence will never remove tare intercept terms or
%linear voltage from corresponding channel
fprintf('\n     Starting Phase 1 of optimization.... \n')

VIFthresh_elev=max(50, VIFthresh); %Temporary elevated VIF threshold
[customMatrix_P1, P1_maxVIF_hist, P1_maxVIF_init, optChannel] = VIF_constraint(VIFthresh_elev, customMatrix, loaddimFlag, nterms, comIN0, optChannel, high_con, high); %Remove terms based on VIF constraint

%% Phase 2: Sequentially remove insignificant terms until all terms (not necessarily including tare intercepts) are significant
%Sequence will never remove tare intercept terms or linear voltage from
%corresponding channel. If tare intercept terms are not significant,
%warning will be provided
fprintf('\n     Starting Phase 2 of optimization.... \n')

customMatrix_P2=customMatrix_P1; %Initialize as final Phase 1 matrix.  Will trim down terms from there
P2_maxP_hist=zeros(nterms,loaddimFlag); %variable for tracking history of max p
P2_maxP_init=zeros(1,loaddimFlag); %vector for initial max p


for i=1:loaddimFlag
    if optChannel(i)==1 %Check if channel model should be optimized
        
        %Initial check on term significance
        sig_iter=zeros(size(customMatrix,1),1); %Initialize sig vector
        p_iter=zeros(size(customMatrix,1),1); %Initialize coefficeint p value vector
        [sig_iter(boolean(customMatrix_P2(:,i))),p_iter(boolean(customMatrix_P2(:,i)))]=sig_check(comIN0(:,boolean(customMatrix_P2(:,i))),targetMatrix0(:,i),anova_pct); %Calculate significance for all initial included terms
        if all(sig_iter(boolean(customMatrix_P2(:,i))))
            insig=0; %All terms already significant
        else
            insig=1; %Insignificant terms are present in model
        end
        
        %Loop to remove terms until all terms are significant
        termOut_count=1;
        while insig==1 %Remove terms until all are significant
            p_iter_search=p_iter(1:nterms); %initialize variable for searching for max p
            p_iter_search(i)=0; %Zero out corresponding linear voltage so it will not be selected for removal
            [maxp_iter,maxI]=max(p_iter_search); %Find index of max coefficient p value out of allowable selections (not tare intercepts or linear voltages from corresponding channel)
            
            if all(sig_iter(boolean(customMatrix_P2(1:nterms,i)))) && ~all(sig_iter(boolean(customMatrix_P2(:,i)))) %Check if all allowed terms have been removed (possible insignificant tare intercepts)
                fprintf('       Error occured in Phase 2.');
                fprintf(' Possible insignificant tare intercept in channel: '); fprintf(num2str(i)); fprintf('.\n');
                break
            elseif maxp_iter==0 && ~all(sig_iter(boolean(customMatrix_P2(:,i)))) %Check if all allowed terms have been removed
                fprintf('       Error occured in Phase 2.');
                fprintf(' Unable to optimize channel: '); fprintf(num2str(i)); fprintf('.\n');
                optChannel(i)=0; %Mark unable to optimize for channel
                break
            end
            
            customMatrix_P2(maxI,i)=0; %Eliminate term with highest p value
            if high_con==1 %If enforcing hierarchy constraint
                customMatrix_P2(boolean(high(:,maxI)),i)=0; %Removes terms that required variable for support
            end
            
            sig_iter=zeros(size(customMatrix,1),1); %Reset sig vector
            p_iter=zeros(size(customMatrix,1),1); %Reset coefficeint p value vector
            [sig_iter(boolean(customMatrix_P2(:,i))),p_iter(boolean(customMatrix_P2(:,i)))]=sig_check(comIN0(:,boolean(customMatrix_P2(:,i))),targetMatrix0(:,i),anova_pct); %Calculate significance for all initial included terms
            P2_maxP_hist(termOut_count,i)=max(p_iter); %Store new max p value
            
            if all(sig_iter(boolean(customMatrix_P2(:,i)))) %If all model terms are significant now
                insig=0; %Flip flag. Constraint has been met
            else
                termOut_count=termOut_count+1; %Advance counter
            end
            
        end
    end
end

%% Phase 3: Sequentially remove terms with highest VIF until max VIF is <
% threshold. Sequence will never remove tare intercept terms or
%linear voltage from corresponding channel
fprintf('\n     Starting Phase 3 of optimization.... \n')

[customMatrix_P3, P3_maxVIF_hist, P3_maxVIF_init, optChannel] = VIF_constraint(VIFthresh, customMatrix_P2, loaddimFlag, nterms, comIN0, optChannel, high_con, high); %Remove terms based on VIF constraint

%% Phase 4: Search for model within constraints that optimizes search criteria
fprintf('\n     Starting Phase 4 of optimization....\n ')

customMatrix_P4=customMatrix_P3; %Initialize as final Phase 3 matrix.  Will trim down terms from there
P3NTerms=sum(customMatrix_P3,1); %Current number of terms included in each channel
minNTerms=(size(customMatrix_P3,1)-nterms)+1; %Min number of terms in each channel (must leave tare intercepts and linear voltage)

cTerms_idx_hist=cell(max(P3NTerms),loaddimFlag);
searchMetric_hist=zeros(max(P3NTerms),loaddimFlag);
test_count_best=zeros(1,loaddimFlag); 
for i=1:loaddimFlag %Loop through each load channel
    if optChannel(i)==1 %If channel should be optimized
        nterms_iter=sum(customMatrix_P4(:,i)); %Count of current number of terms included
        cTerms_idx=find(customMatrix_P4(:,i)); %Index of current included terms
        searchMetric_hist(1,i)=press_check(comIN0(:,boolean(customMatrix_P4(:,i))),targetMatrix0(:,i)); %Check PRESS for initial model
        cTerms_idx_hist{1,i}=cTerms_idx; %Store initial index of included terms
        
        test_count=2;
        while nterms_iter>minNTerms %Test models down to minimum number of terms
            cTerms_idx_test=cTerms_idx; %Duplicate for variable to test removing terms
            cTerms_idx_test(cTerms_idx_test>nterms)=[]; %Don't test removing tare intercepts
            cTerms_idx_test(cTerms_idx_test==i)=[]; %Don't test removing corresponding channel linear voltage
            
            searchMetric_comp=zeros(numel(cTerms_idx_test),1); %Initialize variable for storing search metric for each iteration
            customMatrix_comp=cell(numel(cTerms_idx_test),1); %Initialize cells for storing channel custom matrix for each iteration
            for j=1:numel(cTerms_idx_test) %Loop through, testing each term possible to remove
                customMatrix_iter=customMatrix_P4(:,i); %Initialize as current P4 matrix
                customMatrix_iter(cTerms_idx_test(j))=0; %Eliminate selected term
                if high_con==1 %If enforcing hierarchy constraint
                    customMatrix_iter(boolean(high(:,cTerms_idx_test(j))))=0; %Removes terms that required variable for support
                end
                customMatrix_comp{j}=customMatrix_iter; %Store custom matrix
                searchMetric_comp(j)=press_check(comIN0(:,boolean(customMatrix_iter)),targetMatrix0(:,i)); %Check PRESS for term combination
            end
            
            [best_iter, j_best]=min(searchMetric_comp); %Find Best (min) value from search
            customMatrix_P4(:,i)=customMatrix_comp{j_best}; %Store custom matrix for new best result
            cTerms_idx=find(customMatrix_P4(:,i)); %Index of current included terms

            searchMetric_hist(test_count,i)=best_iter; %Store search metric value in history variable
            cTerms_idx_hist{test_count,i}=cTerms_idx; %Store new index of included terms
            
            nterms_iter=numel(cTerms_idx); %Current count of included terms
            test_count=test_count+1; %Advance count
        end
        [best_mod,test_count_best(i)]=min(searchMetric_hist(1:test_count-1,i)); %Find best overall model
        customMatrix_P4(:,i)=zeros(size(customMatrix_P4,1),1); %Reset channel custom matrix
        customMatrix_P4(cTerms_idx_hist{test_count_best(i),i},i)=1; %Select terms to include
        
    end
end

customMatrix_rec=customMatrix; %Initialize recommended custom matrix as provided custom matrix
customMatrix_rec(:,boolean(optChannel))=customMatrix_P4(:,boolean(optChannel)); %Set optimized channels to P4 Results

fprintf('Recommended Equation Search Complete. \n ')
end

function [customMatrix_Px, Px_maxVIF_hist, Px_maxVIF_init, optChannel] = VIF_constraint(VIFthresh, customMatrix, loaddimFlag, nterms, comIN0, optChannel, high_con, high)
%% Phase 1 and 3: Sequentially remove terms with highest VIF until max VIF is <
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
        vif_iter(boolean(customMatrix_Px(:,i)))=vif_dl(comIN0(:,boolean(customMatrix_Px(:,i)))); %Calculate VIF for all initial included terms
        Px_maxVIF_init(i)=max(vif_iter); %Store initial max VIF
        if Px_maxVIF_init(i)>VIFthresh %If initial max VIF is over elevated threshold
            vifHigh=1; %Flag for if max VIF is too high
        else
            vifHigh=0; %Elevated VIF threshold already satisfied
        end
        
        %Loop to remove terms until VIF threshold is met
        termOut_count=1;
        while vifHigh==1 %Remove terms until VIF threshold is met
            vif_iter_search=vif_iter(1:nterms); %initialize variable for searching for max VIF
            vif_iter_search(i)=0; %Zero out corresponding linear voltage so it will not be selected for removal
            [maxVIF_iter,maxI]=max(vif_iter_search); %Find index of max VIF out of allowable selections (not tare intercepts or linear voltages from corresponding channel)
            
            if maxVIF_iter==0 %Check if all terms have been removed
                fprintf('       Error occured in Phase.');
                fprintf(' Check for linear dependence between series and gage output for Channel: '); fprintf(num2str(i)); fprintf('.\n');
                optChannel(i)=0; %Mark unable to optimize for channel
                break
            end
            customMatrix_Px(maxI,i)=0; %Eliminate term with highest VIF
            if high_con==1 %If enforcing hierarchy constraint
                customMatrix_Px(boolean(high(:,maxI)),i)=0; %Removes terms that required variable for support
            end
            
            vif_iter=zeros(size(customMatrix,1),1); %Reset VIF vector
            vif_iter(boolean(customMatrix_Px(:,i)))=vif_dl(comIN0(:,boolean(customMatrix_Px(:,i)))); %Calculate VIF for all new included terms
            Px_maxVIF_hist(termOut_count,i)=max(vif_iter); %Store new max VIF
            
            if calcThrough~=loaddimFlag && maxI<= loaddimFlag %If calculating not full number of channels (for time savings) but eliminated a linear voltage term
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
