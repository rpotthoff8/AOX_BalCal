function [xcalib,ANOVA]=calc_xcalib(comIN,targetMatrix,series,nterms,nseries0,dimFlag,FLAGS,customMatrix, anova_pct, labels,method)
%Function calculates coefficient matrix (xcalib)

%Orders data by series (needed for bootstrap)
[series,sortI]=sort(series);
comIN=comIN(sortI,:);
targetMatrix=targetMatrix(sortI,:);

% Characterizes the series in the subsamples
[~,s_1st,~] = unique(series);
nseries = length(s_1st);

% Indirect approach uses the modeled voltage
%     if approach_FLAG == 1
%         if balCal_FLAG == 2
%             excessVec0 = qtaprxINminGZ2 + globalZerosAllPoints;
%         else
%             excessVec0 = qtaprxINminGZ + globalZerosAllPoints;
%         end
%     end

xcalib = zeros(nterms+nseries0,dimFlag);
% Solves for the coefficient one column at a time.
% This is to account for Custom Models, where the terms may be
% different depending on the channel.
for k = 1:dimFlag
    comIN_k = comIN;

    if FLAGS.model == 4
        comIN_k(:,customMatrix(:,k)==0) = [];
    end

    % SOLUTION
    if nseries==nseries0
        lastwarn('');
        xcalib_k = comIN_k\targetMatrix(:,k); %CHANGED
        [~,warnID]=lastwarn;
        if strcmp(warnID,'MATLAB:rankDeficientMatrix')==1 || strcmp(warnID,'MATLAB:nearlySingularMatrix')==1
            fprintf('Nearly Singular Result solving for coefficients using ''\\'': Using ''pinv'' instead.  Expect longer computation times. \n');
            xcalib_k = pinv(comIN_k)*targetMatrix(:,k); %CHANGED
%             xcalib_k_3=pinv(comIN_k,10^-1)*targetMatrix(:,k);
%             xcalib_k_4=lsqminnorm(comIN_k,targetMatrix(:,k));
        end
    else
        xcalib_k = pinv(comIN_k)*targetMatrix(:,k);
    end
    if FLAGS.model == 4
        xcalib(customMatrix(:,k)==1,k) = xcalib_k;
    else
        xcalib(:,k) = xcalib_k;
    end

    %Call Anova
    if FLAGS.anova==1
        fprintf(['\nCalculating ', method,' ANOVA statistics for channel ', num2str(k), ' (',labels{k},')....\n'])
        ANOVA(k)=anova(comIN_k,targetMatrix(:,k),nseries0,0,anova_pct);
        fprintf('Complete')
    end

end
fprintf('\n')

if FLAGS.anova==0
    ANOVA='ANOVA NOT PERFORMED';
else
    if FLAGS.model==4 %If custom equation, expand ANOVA statistics to standard 96 term matrix
        ANOVA_exp=ANOVA;
        ExpandList=["beta","beta_CI","T","p_T","VIF","sig"]; %List of ANOVA structure elements that should be expanded
        for i=1:size(ExpandList,2)
            for j=1:dimFlag
                eval(strcat('ANOVA_exp(',num2str(j),').',ExpandList(i),'=zeros(size(xcalib,1),1);')); %initialize zeros
                eval(strcat('ANOVA_exp(',num2str(j),').',ExpandList(i),'(customMatrix(:,j)==1,:)=ANOVA(',num2str(j),').',ExpandList(i),';')); %fill with ANOVA statistics
            end
        end
        %Expand to 96x96 matrix for invXtX
        for j=1:dimFlag
            ANOVA_exp(j).PI.invXtX=zeros(nterms,nterms);
            ANOVA_exp(j).PI.invXtX(customMatrix((1:nterms),j)==1,customMatrix((1:nterms),j)==1)=ANOVA(j).PI.invXtX;
        end
   
        ANOVA=ANOVA_exp;
    end

end
end
