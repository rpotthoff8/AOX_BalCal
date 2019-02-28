function[]=RBFcarlo(nCarlo,excessVec, dalzcalib,globalZerosAllPoints,centerIndexHist,wHist,cHist,aprxINminGZ,series,targetMatrix,s_1st)
%ADDED 9 Jan 19 JRP: Monte carlo
    
     
    % %Define magnitude of noise:
    % percentVoltage=0.01; %percent of max recorded voltage that is noise
    % Maxnoise=percentVoltage*max(abs(dainputs2));

    
    %equally distributed noise
    % sigma=Maxnoise./2; % 95% of the noise points will be within the defined Maxnoise
    
    sigma=1/2; %Trust all channels down to 1 microvolt, sigma is trust/2
    
    
    
    aprxINminGZ2_storage=zeros(size(aprxINminGZ,1),size(aprxINminGZ,2),nCarlo);
    taresGRBF_storage=zeros(size(s_1st,1),size(targetMatrix,2),nCarlo);
    
    for carloCount=1:nCarlo
        aprxINminGZ2 = aprxINminGZ;
        % Normally distributed noise
        noise=sigma.*randn(size(excessVec,1),size(excessVec,2));
        
        % % Equally distributed Noise
        % noise=-Maxnoise+(2*Maxnoise).*rand(nCarlo,size(localZeros_da,2));
        
        excessVecCarlo=excessVec+noise;
        
        for i=1:dimFlag
            dainputscalib(:,i) = excessVecCarlo(:,i)-globalZeros(i);
%             dalzcalib(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
        end
        
        %    localZeroMatrix = localZerosAllPoints;
%         globalZerosAllPoints = zeros(length(excessVec(:,1)),dimFlag); % ajm 6_2_18
        
        etaLZ = dot(dalzcalib-dainputscalib,dalzcalib-dainputscalib);
        etaGZ = dot(globalZerosAllPoints-dainputscalib,globalZerosAllPoints-dainputscalib);
        
        for u=1:numBasis
            w=wHist(:,u);
            centerIndexLoop= centerIndexHist(u,:);
            coeff=cHist(u,:);
            for s=1:dimFlag
                
%                 [goopLoop(s),centerIndexLoop(s)] = max(abs(targetRes2(:,s)));
                
                for r=1:length(excessVec(:,1))
                    eta(r,s) = dot(dainputscalib(r,:)-dainputscalib(centerIndexLoop(s),:),dainputscalib(r,:)-dainputscalib(centerIndexLoop(s),:));
                end
                
%                 %find widths 'w' by optimization routine
%                 w(s) = fminbnd(@(w) balCal_meritFunction2(w,targetRes2(:,s),eta(:,s)),0,1 );
%                 
                rbfINminLZ(:,s)=exp(eta(:,s)*log(abs(w(s)))) - exp(etaLZ(:,s)*log(abs(w(s))));
                rbfINminGZ(:,s)=exp(eta(:,s)*log(abs(w(s))));
                rbfLZminGZ(:,s)=exp(etaLZ(:,s)*log(abs(w(s))));%to find tares AAM042016
                
%                 coeff(s) = dot(rbfINminGZ(:,s),targetRes2(:,s)) / dot(rbfINminGZ(:,s),rbfINminGZ(:,s));
                
                rbfc_INminLZ(:,s) = coeff(s)*rbfINminLZ(:,s);
                rbfc_INminGZ(:,s) = coeff(s)*rbfINminGZ(:,s);
                rbfc_LZminGZ(:,s) = coeff(s)*rbfLZminGZ(:,s); %to find tares AAM042016
            end
            
%         etaHist{u} = eta;

        %update the approximation
        aprxINminGZ2 = aprxINminGZ2+rbfc_INminGZ;
        
        % SOLVE FOR TARES BY TAKING THE MEAN
        taretalGRBF = meantare(series,aprxINminGZ2-targetMatrix);

        taresGRBF = taretalGRBF(s_1st,:);

%         tareGRBFHist{u} = taresGRBF;

%         targetRes2 = targetMatrix-aprxINminGZ2+taretalGRBF;      %0=b-Ax
%         newRes2 = targetRes2'*targetRes2;
%         resSquare2 = diag(newRes2);
%         resSquareHist(u,:) = resSquare2;

        end
        aprxINminGZ2_storage(:,:,carloCount)=aprxINminGZ2;
        taresGRBF_storage(:,:,carloCount)=taresGRBF;
           
    end
end
    %END added for monte-carlo