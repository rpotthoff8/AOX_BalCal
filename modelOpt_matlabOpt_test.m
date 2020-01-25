function [customMatrix_rec]= modelOpt_matlabOpt_test(VIFthresh, customMatrix_permit, loaddimFlag, nterms, comIN0, anova_pct, targetMatrix0, high, FLAGS)
% Function searches for the 'recommended equation' using the approach
% Balfit reference B16 describes as "Forward Selection".
% Function determines optimal combination of terms to include between
% possible upper and lower bounds.
% Upper Bounds: provided current permitted customMatrix
% Lower Bounds: Linear voltage from channel and tares

%See Balfit reference B9 for more complete explanation.
%Flowchart on page 10/21

%INPUTS:
%  VIFthresh = Threshold for max allowed VIF (Balfit Search constraint 2)
%  customMatrix = Current matrix of what terms should be included in Eqn Set.
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

high_con=FLAGS.high_con;
search_metric_flag=FLAGS.search_metric;
VIF_stop_flag=FLAGS.VIF_stop;

fprintf('\nCalculating Recommended Eqn Set with Forward Selection Method....')

optChannel=ones(1,loaddimFlag); %Flag for if each channel should be optimized

% Normalize the data for a better conditioned matrix (Copy of what is done
% in calc_xcalib.m)
scale = max(abs(comIN0));
scale(scale==0)=1; %To avoid NaN for channels where RBFs have self-terminated
comIN0 = comIN0./scale;

%Set lower bound of search (Required Math Model)
customMatrix_req=zeros(size(customMatrix_permit));
customMatrix_req(nterms+1:end,:)=1; %Must include series intercepts
customMatrix_req(1:loaddimFlag,1:loaddimFlag)=eye(loaddimFlag); %Must include linear voltage from channel

%Define number of series included:
nseries=size(customMatrix_permit,1)-nterms;

%Initialize optimized math model as required (lower end)
customMatrix_opt=customMatrix_req;

num_permit=sum(customMatrix_permit); %Max number of terms
num_req=sum(customMatrix_req); %Min number of terms
num_terms=sum(customMatrix_opt); %Current count of terms
num_test=(num_permit-num_req)+1; %Number of models to find. 1 for every number of terms from current # to max #

%Initialize variables


for i=1:loaddimFlag %Loop through all channels
    if optChannel(i)==1 %If optimization is turned on

        x=optimvar('x',num_test(i)-1,1,'Type','integer','LowerBound',0,'UpperBound',1);
        
        %Possible terms to be added are those not in current model that are in permitted model
        pos_add=zeros(size(customMatrix_opt,1),1); %Initialize as zeros
        pos_add(~boolean(customMatrix_opt(:,i)))=customMatrix_permit(~boolean(customMatrix_opt(:,i)),i); %Vector of terms that can be added
        pos_add_idx=find(pos_add); %Index of possible terms for adding to model
        
        objFunc_iter=@(x_in) objFunc(x_in, comIN0,customMatrix_req(:,i),pos_add_idx, targetMatrix0(:,i), anova_pct, VIFthresh, nseries, search_metric_flag);
        
        [searchMet,vif_c, sig_c] = fcn2optimexpr(objFunc_iter,x,'ReuseEvaluation',true,'OutputSize',[1 1 1]);
        
        vif_cons= vif_c==1;
        sig_cons= sig_c==1;
        
        prob = optimproblem('Objective', searchMet);
        prob.Constraints.vif_cons=vif_cons;
        prob.Constraints.sig_cons=sig_cons;
        
%         show(prob);
        
%         x0.x=zeros(size(x));
%         [sol,fval,exitflag,output] = solve(prob);
        
        sol=ga(objFunc_iter, num_test(i)-1,[],[],[],[],zeros(num_test(i)-1,1),ones(num_test(i)-1,1),[],([1: num_test(i)-1]));

        %Store new custom equation
        customMatrix_temp=customMatrix_req(:,i); %Initialize as required equation
        customMatrix_temp(pos_add_idx(boolean(sol)))=1; %Turn on terms based on optimization
        customMatrix_opt(:,i)=customMatrix_temp; %Store as optimal 
    end
end
%Output final model:
customMatrix_rec=customMatrix_permit; %Initialize recommended custom matrix as provided custom matrix
customMatrix_rec(:,boolean(optChannel))=customMatrix_opt(:,boolean(optChannel)); %Set optimized channels to optimal Results

fprintf('\nRecommended Equation Search Complete. \n ')
end

function [search_metric, VIF_met, sig_all]= objFunc(x, comIN0,customMatrix_req,pos_add_idx, targetMatrix0_chan, anova_pct, VIFthresh, nseries, search_metric_flag)

customMatrix_temp=customMatrix_req;
customMatrix_temp(pos_add_idx(boolean(x)))=1;
[VIF_met, ~, sig_all, ~ ,search_metric]=test_combo_matlabOpt(comIN0(:,boolean(customMatrix_temp)),targetMatrix0_chan, anova_pct, VIFthresh, nseries, search_metric_flag, 2);


end
