function [customMatrix_rec]= modelOpt_suggestedNew(VIFthresh, customMatrix, loaddimFlag, nterms, comIN0, anova_pct, targetMatrix0, high, high_con)
% Searches for the 'suggested math model'.  This is an improved version of
% the approach described in Balfit Reference B29.  The "suggested" math
% model search is a faster search process than searching for the
% "recommended" math model. It does not optimize as search metric, but instead just searches for an
% equation where there are no near-linear dependencies between math terms,
% the hierarchy rule is satisfied (optional), and all terms are
% statistically significant.  It searches between the upper and lower
% bounds for terms in the model:
% Upper Bounds: provided current customMatrix
% Lower Bounds: Linear voltage from channel and tares

% Hierarchy can be enforced in 2 ways:
% high_con==1: Enforced following search: Terms are removed to meet VIF and
%       term significance constraints, then terms required to support remaining
%       terms are added back in.  This may result in a final math model where all
%       terms do not satisfy the VIF and/or significance constraints
% high_con==2: Enforced during search: Terms are removed to meet VIF and
%       term significance constraints. As terms are removed, all terms
%       requiring term for support are also removed. Results in final math
%       model where all terms satisfy VIF and significance constraints,
%       although generally fewer terms will be included and may result in
%       higher residuals 

%The following process is used in this function:
% 1) Start with the provided 'permitted math model', obtained by applying
%       SVD. This model is fed into the current function
% 2) PHASE 1: Apply VIF constraint with elevated threshold: Compute VIF for all terms, remove term with
%       highest VIF, recompute VIF, repeat until max VIF is below elevated
%       threshold (50)
% 3) (Optional): Add missing lower order terms (high_con==1)
% 4) PHASE 2: Remove statistically 'insignificant' math terms: Find P values of all
%       regression coefficients from terms from previous step.  Remove term with
%       highest P-value.  Recompute, repeat until all all terms are
%       significant based on a lowered percent confidence (75%)
% 5) (Optional): Add missing lower order terms (high_con==1)
% 6) PHASE 3: Apply VIF constraint with final threshold: Compute VIF for all terms, remove term with
%       highest VIF, recompute VIF, repeat until max VIF is below final threshold
% 7) (Optional): Add missing lower order terms (high_con==1)
% 8) PHASE 4: Remove statistically 'insignificant' math terms: Find P values of all
%       regression coefficients from terms from previous step.  Remove term with
%       highest P-value.  Recompute, repeat until all all terms are
%       significant based on a final percent confidence
% 9) (Optional): Add missing lower order terms (high_con==1)

%INPUTS:
%  VIFthresh = Threshold for max allowed VIF (Balfit Search constraint 2)
%  customMatrix = Current matrix of what terms should be included in Eqn Set.
%  loaddimFlag = Dimension of load data (# channels)
%  nterms = Number of predictor terms in regression model
%  comIN0 = Matrix of predictor variables
%  anova_pct = ANOVA percent confidence level for determining significance
%  targetMatrix0 = Matrix of target load values
%  high = Matrix of term hierarchy
%  high_con = Flag for if term hierarchy constraint is enforced. 0==off, 1==after search, 2==during search

%OUTPUTS:
%  customMatrix_rec = Optimized recommended custom matrix

fprintf('\nCalculating Suggested Eqn Set with updated efficient search method....')

optChannel=ones(1,loaddimFlag); %Flag for if each channel should be optimized

%Set Minimum equation (lower bound of search)
min_eqn=zeros(size(customMatrix));
min_eqn(nterms+1:end,:)=1; %Must include series intercepts
min_eqn(1:loaddimFlag,1:loaddimFlag)=eye(loaddimFlag); %Must include linear voltage from channel

% Normalize the data for a better conditioned matrix (Copy of what is done
% in calc_xcalib.m)
scale = max(abs(comIN0));
scale(scale==0)=1; %To avoid NaN for channels where RBFs have self-terminated
comIN0 = comIN0./scale;

%% Phase 1: Sequentially remove terms with highest VIF until max VIF is <
%elevated threshold. Sequence will never remove tare intercept terms or
%linear voltage from corresponding channel
fprintf('\n     Starting Phase 1 of optimization.... \n')

VIFthresh_elev=max([50, VIFthresh]); %Temporary elevated VIF threshold
[customMatrix_P1, P1_maxVIF_hist, P1_maxVIF_init, optChannel] = VIF_constraint(VIFthresh_elev, customMatrix, loaddimFlag, nterms, comIN0, optChannel, high_con, high, min_eqn); %Remove terms based on VIF constraint

if high_con==1 %If enforcing hierarchy constraint after
    [customMatrix_P1]= high_add(customMatrix_P1, high, optChannel, loaddimFlag, nterms); %Call function to add in missing terms
end

%% Phase 2: Sequentially remove insignificant terms until all terms (not necessarily including tare intercepts) are significant
%Sequence will never remove tare intercept terms or linear voltage from
%corresponding channel. If tare intercept terms are not significant,
%warning will be provided
fprintf('\n     Starting Phase 2 of optimization.... \n')

anova_pct_low=min([75, anova_pct]);
[customMatrix_P2, P2_maxP_hist, optChannel] = sig_constraint(anova_pct_low, customMatrix_P1, loaddimFlag, nterms, comIN0, targetMatrix0, optChannel, high_con, high, min_eqn);

if high_con==1 %If enforcing hierarchy constraint after
    [customMatrix_P2]= high_add(customMatrix_P2, high, optChannel, loaddimFlag, nterms); %Call function to add in missing terms
end

%% Phase 3: Sequentially remove terms with highest VIF until max VIF is <
% threshold. Sequence will never remove tare intercept terms or
%linear voltage from corresponding channel
fprintf('\n     Starting Phase 3 of optimization.... \n')

[customMatrix_P3, P3_maxVIF_hist, P3_maxVIF_init, optChannel] = VIF_constraint(VIFthresh, customMatrix_P2, loaddimFlag, nterms, comIN0, optChannel, high_con, high, min_eqn); %Remove terms based on VIF constraint

if high_con==1 %If enforcing hierarchy constraint after
    [customMatrix_P3]= high_add(customMatrix_P3, high, optChannel, loaddimFlag, nterms); %Call function to add in missing terms
end

%% Phase 4: Sequentially remove insignificant terms until all terms (not necessarily including tare intercepts) are significant
%Sequence will never remove tare intercept terms or linear voltage from
%corresponding channel. If tare intercept terms are not significant,
%warning will be provided
fprintf('\n     Starting Phase 4 of optimization.... \n')

[customMatrix_P4, P2_maxP_hist, optChannel] = sig_constraint(anova_pct, customMatrix_P3, loaddimFlag, nterms, comIN0, targetMatrix0, optChannel, high_con, high, min_eqn);

if high_con==1 %If enforcing hierarchy constraint after
    [customMatrix_P4]= high_add(customMatrix_P4, high, optChannel, loaddimFlag, nterms); %Call function to add in missing terms
end

customMatrix_rec=customMatrix; %Initialize recommended custom matrix as provided custom matrix
customMatrix_rec(:,boolean(optChannel))=customMatrix_P4(:,boolean(optChannel)); %Set optimized channels to P4 Results

fprintf('Suggested Equation Search Complete. \n ')
end

function [customMatrix_Px, Px_maxVIF_hist, Px_maxVIF_init, optChannel] = VIF_constraint(VIFthresh, customMatrix, loaddimFlag, nterms, comIN0, optChannel, high_con, high, min_eqn)
%% Phase 1 and 3: Sequentially remove terms with highest VIF until max VIF is <
%threshold. Sequence will never remove minimum equation terms  from corresponding channel
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
            vif_iter_search=vif_iter; %initialize variable for searching for max VIF
            vif_iter_search(boolean(min_eqn(:,i)))=0; %Zero out corresponding minimum eqn so it will not be selected for removal
            [maxVIF_iter,maxI]=max(vif_iter_search); %Find index of max VIF out of allowable selections (not tare intercepts or linear voltages from corresponding channel)
            
            if maxVIF_iter<=VIFthresh %Check if all 'problem' terms have been removed
                fprintf('       Check for linear dependence between minimum required equation for Channel: '); fprintf(num2str(i)); fprintf('.\n');
                break
            elseif maxVIF_iter== 0 
                fprintf('       Error occured in Phase.');
                fprintf(' Check for severe linear dependence minimum required equation for Channel: '); fprintf(num2str(i)); fprintf('.\n');
                optChannel(i)=0; %Mark unable to optimize for channel
                break
            end
            customMatrix_Px(maxI,i)=0; %Eliminate term with highest VIF
            if high_con==2 %If enforcing hierarchy constraint
                customMatrix_Px(boolean(high(:,maxI)),i)=0; %Removes terms that required variable for support
            end
            
            vif_iter=zeros(size(customMatrix,1),1); %Reset VIF vector
            vif_iter(boolean(customMatrix_Px(:,i)))=vif_dl(comIN0(:,boolean(customMatrix_Px(:,i)))); %Calculate VIF for all new included terms
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

function [customMatrix_Px, Px_maxP_hist, optChannel] = sig_constraint(anova_pct, customMatrix, loaddimFlag, nterms, comIN0, targetMatrix0, optChannel, high_con, high, min_eqn)
%% Phase 2 and 4: Sequentially remove insignificant terms until all terms (not necessarily including tare intercepts) are significant
%Sequence will never remove tare intercept terms or linear voltage from
%corresponding channel. If tare intercept terms are not significant,
%warning will be provided

customMatrix_Px=customMatrix; %Initialize as final Phase 1 matrix.  Will trim down terms from there
Px_maxP_hist=zeros(nterms,loaddimFlag); %variable for tracking history of max p
Px_maxP_init=zeros(1,loaddimFlag); %vector for initial max p


for i=1:loaddimFlag
    if optChannel(i)==1 %Check if channel model should be optimized
        
        %Initial check on term significance
        sig_iter=zeros(size(customMatrix,1),1); %Initialize sig vector
        p_iter=zeros(size(customMatrix,1),1); %Initialize coefficeint p value vector
        [sig_iter(boolean(customMatrix_Px(:,i))),p_iter(boolean(customMatrix_Px(:,i)))]=sig_check(comIN0(:,boolean(customMatrix_Px(:,i))),targetMatrix0(:,i),anova_pct); %Calculate significance for all initial included terms
        if all(sig_iter(boolean(customMatrix_Px(:,i))))
            insig=0; %All terms already significant
        else
            insig=1; %Insignificant terms are present in model
        end
        
        %Loop to remove terms until all terms are significant
        termOut_count=1;
        while insig==1 %Remove terms until all are significant
            p_iter_search=p_iter; %initialize variable for searching for max p
            p_iter_search(boolean(min_eqn(:,i)))=0; %Zero out corresponding min eqn terms it will not be selected for removal
            [maxp_iter,maxI]=max(p_iter_search); %Find index of max coefficient p value out of allowable selections (not tare intercepts or linear voltages from corresponding channel)
            
            if all(sig_iter(boolean(customMatrix_Px(1:nterms,i)))) && ~all(sig_iter(boolean(customMatrix_Px(:,i)))) %Check if all allowed terms have been removed (possible insignificant tare intercepts)
                fprintf('       Error occured in Phase 2.');
                fprintf(' Possible insignificant tare intercept in channel: '); fprintf(num2str(i)); fprintf('.\n');
                break
            elseif maxp_iter==0 && ~all(sig_iter(boolean(customMatrix_Px(:,i)))) %Check if all allowed terms have been removed
                fprintf('       Error occured in Phase 2.');
                fprintf(' Unable to optimize channel: '); fprintf(num2str(i)); fprintf('.\n');
                optChannel(i)=0; %Mark unable to optimize for channel
                break
            end
            
            customMatrix_Px(maxI,i)=0; %Eliminate term with highest p value
            if high_con==2 %If enforcing hierarchy constraint
                customMatrix_Px(boolean(high(:,maxI)),i)=0; %Removes terms that required variable for support
            end
            
            sig_iter=zeros(size(customMatrix,1),1); %Reset sig vector
            p_iter=zeros(size(customMatrix,1),1); %Reset coefficeint p value vector
            [sig_iter(boolean(customMatrix_Px(:,i))),p_iter(boolean(customMatrix_Px(:,i)))]=sig_check(comIN0(:,boolean(customMatrix_Px(:,i))),targetMatrix0(:,i),anova_pct); %Calculate significance for all initial included terms
            Px_maxP_hist(termOut_count,i)=max(p_iter); %Store new max p value
            
            if all(sig_iter(boolean(customMatrix_Px(:,i)))) %If all model terms are significant now
                insig=0; %Flip flag. Constraint has been met
            else
                termOut_count=termOut_count+1; %Advance counter
            end
            
        end
    end
end

end

function [customMatrix_Px]= high_add(customMatrix, high, optChannel, loaddimFlag, nterms)
% Apply 'Hierarchy Rule': Add missing lower order terms to customMatrix
% from previous phase
customMatrix_Px=customMatrix; %Initialize as previous phase result

identical_custom=all(~diff(customMatrix_Px,1,2),'all'); %Check if customMatrix is identical for all channels
if identical_custom==1 %If custom equation identical for each channel
    calcThrough=1; %only necessary to remove enforce hierarchy rule once
else
    calcThrough=loaddimFlag; %Necessary to calculate for each channel seperate
end

for i=1:calcThrough %Loop through channels
    if optChannel(i)==1 %Check if channel model should be optimized
        %         inc_idx=find(customMatrix_P2(:,i)); %Find index of included terms for channel
        sup_terms_mat=high(boolean(customMatrix_Px(1:nterms,i)),:); %Rows from hierarchy matrix for included terms. columns with '1' are needed to support variable
        sup_terms=any(sup_terms_mat,1); %Row vector with 1s for terms needed to support included terms
        customMatrix_Px(boolean(sup_terms),i)=1; %Include all terms needed to support currently included terms
    end
end

if calcThrough==1 %If final equation identical for each channel
    customMatrix_Px(:,2:loaddimFlag)=repmat(customMatrix_Px(:,1),1,loaddimFlag-1); %Duplicate for each column
end

end