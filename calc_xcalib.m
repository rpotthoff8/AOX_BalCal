function [xcalib,ANOVA]=calc_xcalib(comIN,targetMatrix,series,nterms,nseries0,dimFlag,FLAGS,customMatrix, anova_pct, labels,method,calc_channel)
%Function calculates coefficient matrix (xcalib)
% calc_channel is used as a flag to determine if coefficients in that
% channel should be calculated, used for RBF self termination

%INPUTS:
%  comIN = Matrix of predictor variables. Each row is observation, each column is predictor variable
%  targetMatrix = Matrix of target values (loads)
%  series = Series labels for each point
%  nterms = Number of predictor terms in regression model
%  nseries0 = Number of series
%  dimFlag = Dimension (# of Channels) for data
%  FLAGS = Structure containing flags for user preferences
%  customMatrix = Matrix of 1's and 0's for which predictor variables should be included in regression model for each channel
%  anova_pct = Percent confidence level for ANOVA calculations
%  labels = Load labels for each channel
%  method = String for current section (ALG or RBF)
%  calc_channel = Boolean vector for which channels coefficients should be solved for

%OUTPUTS:
%  xcalib = Coefficient Matrix
%  ANOVA = Results of ANOVA calculations 

if exist('calc_channel','var')==0 %If no variable provided for which channels to calculate
    calc_channel=ones(1,dimFlag); %Calculate all channels
end

%Orders data by series
[series,sortI]=sort(series);
comIN=comIN(sortI,:);
targetMatrix=targetMatrix(sortI,:);

% Normalize the data for a better conditioned matrix
scale = max(abs(comIN));
scale(scale==0)=1; %To avoid NaN for channels where RBFs have self-terminated
comIN = comIN./scale;

% Characterizes the series in the subsamples
[~,s_1st,~] = unique(series);
nseries = length(s_1st);

xcalib = zeros(nterms+nseries0,dimFlag);
% Solves for the coefficient one column at a time.
% This is to account for Custom Models, where the terms may be
% different depending on the channel.
for k = 1:dimFlag
    if calc_channel(k)==1
        comIN_k = comIN;
        scale_k = scale;
        
        if FLAGS.model == 4
            comIN_k(:,customMatrix(:,k)==0) = [];
            scale_k(customMatrix(:,k)==0) = [];
        end
        
        % SOLUTION
        xcalib_k = comIN_k\targetMatrix(:,k);
        
        % De-normalize the coefficients to be used with raw data
        xcalib_k = xcalib_k./scale_k';
        
        if FLAGS.model == 4
            xcalib(customMatrix(:,k)==1,k) = xcalib_k;
        else
            xcalib(:,k) = xcalib_k;
        end
        
        %Call Anova
        if FLAGS.anova==1
            %test_FLAG used to 'turn off' VIF when iterating to recommended
            %equation for time saving
            if isfield(FLAGS,'test_FLAG')==0
                FLAGS.test_FLAG=0;
            end
            if FLAGS.test_FLAG==0
                fprintf(['\nCalculating ', method,' ANOVA statistics for channel ', num2str(k), ' (',labels{k},')....\n'])
            end
            ANOVA(k)=anova(comIN_k,targetMatrix(:,k),nseries0,FLAGS.test_FLAG,anova_pct);
            
            % There are several ANOVA metrics that also must be denormalized
            ANOVA(k).beta    = ANOVA(k).beta./scale_k';
            ANOVA(k).beta_CI = ANOVA(k).beta_CI./scale_k';
            
            % Prediction interval calculation does not include tares, so scale
            % vector has to be truncated
            scale_PI = scale_k(1:end-nseries);
            ANOVA(k).PI.invXtX = ANOVA(k).PI.invXtX./(scale_PI'*scale_PI);
            if FLAGS.test_FLAG==0
                fprintf('Complete\n')
            end
        end
    end
end
% fprintf('\n')

if FLAGS.anova==0 || all(~calc_channel) %If ANOVA was not calculated or no channels were calculated
    ANOVA='ANOVA NOT PERFORMED';
else
    if FLAGS.model==4 && any(calc_channel) %If custom equation, expand ANOVA statistics to standard full term matrix
        ANOVA_exp=ANOVA;
        ExpandList=["beta","beta_CI","T","p_T","VIF","sig"]; %List of ANOVA structure elements that should be expanded
        for i=1:size(ExpandList,2)
            for j=1:dimFlag
                if calc_channel(j)==1
                    eval(strcat('ANOVA_exp(',num2str(j),').',ExpandList(i),'=zeros(size(xcalib,1),1);')); %initialize zeros
                    eval(strcat('ANOVA_exp(',num2str(j),').',ExpandList(i),'(customMatrix(:,j)==1,:)=ANOVA(',num2str(j),').',ExpandList(i),';')); %fill with ANOVA statistics
                end
            end
        end
        %Expand to full matrix for invXtX
        for j=1:dimFlag
            if calc_channel(j)==1
                ANOVA_exp(j).PI.invXtX=zeros(nterms,nterms);
                ANOVA_exp(j).PI.invXtX(customMatrix((1:nterms),j)==1,customMatrix((1:nterms),j)==1)=ANOVA(j).PI.invXtX;
            end
        end
        
        ANOVA=ANOVA_exp;
    end
    
end
end
