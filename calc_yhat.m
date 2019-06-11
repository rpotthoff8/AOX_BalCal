function [y_hat,ANOVA]=calc_yhat(comIN,targetMatrix,series,nterms,nseries0,dimFlag,model_FLAG,customMatrix, anova_Flag,comINfull)
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
        
        if model_FLAG == 4
            comIN_k(:,customMatrix(:,k)==0) = [];
        end

        % SOLUTION
        if nseries==nseries0
            xcalib_k = comIN_k\targetMatrix(:,k);
        else
            xcalib_k = pinv(comIN_k)*targetMatrix(:,k);
        end
        if model_FLAG == 4
            xcalib(customMatrix(:,k)==1,k) = xcalib_k;
        else
            xcalib(:,k) = xcalib_k;
        end
        
        %Call Anova
        if anova_Flag==1
            ANOVA(k)=anova(comIN_k,targetMatrix(:,k));
        end
    end
        if anova_Flag==0
            ANOVA='ANOVA NOT PERFORMED';
        end
    
        y_hat=comINfull*xcalib;
end 