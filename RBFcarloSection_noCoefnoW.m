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
    
    for carloCount=1:nCarlo
        
        % Normally distributed noise
        noise=sigma.*randn(size(excessVec,1),size(excessVec,2));
        
        % % Equally distributed Noise
        % noise=-Maxnoise+(2*Maxnoise).*rand(nCarlo,size(localZeros_da,2));
        
        excessVecCarlo=excessVec+noise;
        
        for i=1:dimFlag
            dainputscalib(:,i) = excessVecCarlo(:,i)-globalZeros(i);
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
                
%                 %find widths 'w' by optimization routine
%                 w(s) = fminbnd(@(w) balCal_meritFunction2(w,targetRes2(:,s),eta(:,s)),0,1 );
%                 
                rbfINminLZ(:,s)=exp(eta(:,s)*log(abs(w(s)))) - exp(etaLZ(:,s)*log(abs(w(s))));
                rbfINminGZ(:,s)=exp(eta(:,s)*log(abs(w(s))));
                rbfLZminGZ(:,s)=exp(etaLZ(:,s)*log(abs(w(s))));%to find tares AAM042016
                
                coeff(s) = dot(rbfINminGZ(:,s),targetRes2(:,s)) / dot(rbfINminGZ(:,s),rbfINminGZ(:,s));
                
                rbfc_INminLZ(:,s) = coeff(s)*rbfINminLZ(:,s);
                rbfc_INminGZ(:,s) = coeff(s)*rbfINminGZ(:,s);
                rbfc_LZminGZ(:,s) = coeff(s)*rbfLZminGZ(:,s); %to find tares AAM042016
            end
            
        end
        
        
        
        
    end
end
    %END added for monte-carlo