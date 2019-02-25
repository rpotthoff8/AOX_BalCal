%ADDED 9 Jan 19 JRP: Monte carlo
carlo=1;
if carlo==1
    nCarlo=11;
    
    % %Define magnitude of noise:
    % percentVoltage=0.01; %percent of max recorded voltage that is noise
    % Maxnoise=percentVoltage*max(abs(dainputs2));
    Maxnoise=1;
    
    %equally distributed noise
    % sigma=Maxnoise./2; % 95% of the noise points will be within the defined Maxnoise
    
    sigma=1/2; %Trust all channels down to 1 microvolt, sigma is trust/2
    wHist_carlo=cell(numBasis,1);
    cHist_carlo=cell(numBasis,1);
    centerIndexHist_carlo=cell(nCarlo,1);
    etaHist_carlo = cell(nCarlo,1);
    aprxINminGZ_Hist_carlo = cell(nCarlo,1);
    tareHist = cell(nCarlo,1);
    tareGRBFHist_carlo=cell(nCarlo,1);
    resSquareHist_carlo=cell(nCarlo,1);
    for carloCount=1:nCarlo
        
        % Normally distributed noise
        noise=sigma.*randn(size(excessVec,1),size(excessVec,2));
        
        % % Equally distributed Noise
        % noise=-Maxnoise+(2*Maxnoise).*rand(nCarlo,size(localZeros_da,2));
        
        excessVecCarlo=excessVec+noise;
        
        for i=1:dimFlag
            dainputscalib(:,i) = excessVec(:,i)-globalZeros(i);
            dalzcalib(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
        end
        
        %    localZeroMatrix = localZerosAllPoints;
        globalZerosAllPoints = zeros(length(excessVec(:,1)),dimFlag); % ajm 6_2_18
        
        etaLZ = dot(dalzcalib-dainputscalib,dalzcalib-dainputscalib);
        etaGZ = dot(globalZerosAllPoints-dainputscalib,globalZerosAllPoints-dainputscalib);
        
        for u=1:numBasis
            for s=1:dimFlag
                [goopLoop(s),centerIndexLoop(s)] = max(abs(targetRes2(:,s)));
                
                for r=1:length(excessVec(:,1))
                    eta(r,s) = dot(dainputscalib(r,:)-dainputscalib(centerIndexLoop(s),:),dainputscalib(r,:)-dainputscalib(centerIndexLoop(s),:));
                end
                
                %find widths 'w' by optimization routine
                w(s) = fminbnd(@(w) balCal_meritFunction2(w,targetRes2(:,s),eta(:,s)),0,1 );
                
                rbfINminLZ(:,s)=exp(eta(:,s)*log(abs(w(s)))) - exp(etaLZ(:,s)*log(abs(w(s))));
                rbfINminGZ(:,s)=exp(eta(:,s)*log(abs(w(s))));
                rbfLZminGZ(:,s)=exp(etaLZ(:,s)*log(abs(w(s))));%to find tares AAM042016
                
                coeff(s) = dot(rbfINminGZ(:,s),targetRes2(:,s)) / dot(rbfINminGZ(:,s),rbfINminGZ(:,s));
                
                rbfc_INminLZ(:,s) = coeff(s)*rbfINminLZ(:,s);
                rbfc_INminGZ(:,s) = coeff(s)*rbfINminGZ(:,s);
                rbfc_LZminGZ(:,s) = coeff(s)*rbfLZminGZ(:,s); %to find tares AAM042016
            end
            
            wHist(u,:) = w;
            cHist(u,:) = coeff;
            centerIndexHist(u,:) = centerIndexLoop;
            etaHist{u} = eta;
            
            %update the approximation
            aprxINminGZ2 = aprxINminGZ2+rbfc_INminGZ;
            aprxINminGZ_Hist{u} = aprxINminGZ2;
            
            % SOLVE FOR TARES BY TAKING THE MEAN
            taretalGRBF = meantare(series,aprxINminGZ2-targetMatrix);
            
            taresGRBF = taretalGRBF(s_1st,:);
            
            tareGRBFHist{u} = taresGRBF;
            
            targetRes2 = targetMatrix-aprxINminGZ2+taretalGRBF;      %0=b-Ax
            newRes2 = targetRes2'*targetRes2;
            resSquare2 = diag(newRes2);
            resSquareHist(u,:) = resSquare2;
        end
        
        for j=1:numBasis
            wHist_carlo{j}(carloCount,:)=wHist(j,:);
            cHist_carlo{j}(carloCount,:)=cHist(j,:);
        end
%         wHist_carlo{carloCount}=wHist;
%         cHist_carlo{carloCount}=cHist;
        centerIndexHist_carlo{carloCount}=centerIndexHist;
        etaHist_carlo{carloCount}=etaHist;
        aprxINminGZ_Hist_carlo{carloCount}=aprxINminGZ_Hist;
        tareGRBFHist_carlo{carloCount}=tareGRBFHist;
        resSquareHist_carlo{carloCount}=resSquareHist;
    end
    
    
    
    
end
%END added for monte-carlo