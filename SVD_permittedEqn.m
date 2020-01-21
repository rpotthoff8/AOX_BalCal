function customMatrix_permitted=SVD_permittedEqn(customMatrix, voltdimFlag, loaddimFlag, dainputs0, FLAGS, targetMatrix0, series0, voltagelist, zero_threshold, loadCapacities, nterms, nseries0)
% Method uses SVD to determine the "permitted" set of terms for math model
% by enforcing the constraint that no terms (columns of predictor variable matrix)
% are linearly dependant. Mirrors approach used by BalFit to ensure a
% non-singular solution always exists for the global regression problem. 
% Basic flow off process:
% 1) Perform 2 term regression (intercept and voltage) for simple voltage to load model
% 2) Use model to calculate gage output capacities from channel load capacities
% 3) Set threshold for "zero" as a percentage of gage capacity. Temporarily set voltages below threshold to zero
% 4) Assemble predictor term matrix with zeroed original voltage outputs
% 5) Start with first column of predictor variables. Use 'rank' command to determine rank with second column added. If full rank, add column to permitted equation. If not, discard column.
% 6) Repeat in order until all columns have been tested. This results in vector of which terms are not linear dependant and can be included in the model without resulting in a near-singlular matrix
% 7) Repeat for other channels as necessary (if custom equation is different between channels)
% See BalFit references 'B13' (pg 3), 'ReferenceGuide-2019' (pg 17, 40),
% 'B2' (pg 3) for explanation of BalFit's approach

%INPUTS:
%  customMatrix = Current matrix of what terms should be included in Eqn Set. 
%  voltdimFlag = Dimension of voltage data (# channels)
%  loaddimFlag = Dimension of load data (# channels)
%  dainputs0 = Matrix of voltage outputs
%  FLAGS = Structure containing flags for user preferences
%  targetMatrix0 = Matrix of target values (loads)
%  series0 = Series labels for each point
%  voltagelist = Labels for voltage columns
%  zero_threshold = In Balfit: MATH MODEL SELECTION THRESHOLD IN % OF CAPACITY. This variable is used in performing SVD. Datapoints where the gage output (voltage) is less
...then the threshold as a percentage of gage capacity are set to zero for constructing comIN and performing SVD
%  loadCapacities = Maximum load capacity for each channel
%  nterms = Number of predictor terms in regression model
%  nseries0 = Number of series

%OUTPUTS:
% customMatrix_permitted = New custom Eqn matrix for which terms are
% supported by the dataset

fprintf('\nCalculating Permitted Eqn Set with SVD....')

%Perform linear regression to solve for gage capacities
    gageCapacities=zeros(1,voltdimFlag); %initialize
    for i=1:loaddimFlag
        A=[ones(size(dainputs0,1),1),dainputs0(:,i)]; %Predictor variables are just channel voltage and intercept
        B=[targetMatrix0(:,i)]; %Target is loads from channel
        X_lin=A\B; %coefficients for linear model
        gageCapacities(i)=(loadCapacities(i)-X_lin(1))/X_lin(2); %Find gage capacity from load capacity and linear regression model
    end
    if voltdimFlag>loaddimFlag %If greater number of voltage channels than load channels
        gageCapacities(i+1:end)=max(abs(dainputs0(:,i+1:end)),[],1); %In remaining channels where linear regression not possible, set gage capacity as max absolute value gage output
    end

    %Set voltages below threshold to zero:
    volt_theshold=zero_threshold*gageCapacities; %Zero threshold in gage output limits
    dainputs_svd=dainputs0; %Initialize as dainputs
    dainputs_svd(abs(dainputs0)<volt_theshold)=0; %Set voltages less than threshold to 0

    %Construct comIN with thresholded voltages
    comIN_svd = balCal_algEqns(FLAGS.model,dainputs_svd,series0,1,voltagelist); %Matrix of predictor variables from thresholded voltages
    % Normalize the data for a better conditioned matrix
    scale = max(abs(comIN_svd),[],1);
    scale(scale==0)=1; %To avoid NaN
    comIN_svd = comIN_svd./scale;

    identical_custom=all(~diff(customMatrix,1,2),'all'); %Check if customMatrix is identical for all channels
    if identical_custom==1 %If custom equation identical for each channel
        calcThrough=1; %only 1 run through SVD necessary
    else
        calcThrough=loaddimFlag; %Necessary to calculate for each channel seperate
    end

    svd_include=zeros(nterms,loaddimFlag); %Initialize vector for tracking which terms are supported
    svd_include(1,:)=1; %Initialize first term is included
    for j=1:calcThrough
        for i=2:nterms %Loop through all possible terms
            if customMatrix(i,j)==1 %If term is included according to customMatrix for eqn
                svd_include_test=svd_include(:,j); %Initialize test variable for iteration
                svd_include_test(i)=1; %Include new term for this iteration
                rankIter=rank(comIN_svd(:,boolean(svd_include_test))); %Using rank command (SVD) find rank of predictor variable matrix
                if rankIter==sum(svd_include(:,j))+1 %If rank is equal to number of terms
                    svd_include(i,j)=1; %Add term to supported terms
                end
            end
        end
    end

    if identical_custom==1 %If custom equation identical for each channel
        svd_include(:,2:loaddimFlag)=repmat(svd_include(:,1),1,loaddimFlag-1); %Duplicate for each column
    end
    customMatrix_permitted = [svd_include; ones(nseries0,loaddimFlag)]; %customMatrix for permitted eqn set.  Always include intercepts
    fprintf(' Complete. \n')
end