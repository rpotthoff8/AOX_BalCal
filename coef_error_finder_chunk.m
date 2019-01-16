%Take Confidence interval found for coefficients using bootstrap, and store
%the absolute value of the larger error between +/- in a matrix
for i=1:size(xcalib,1)
    for j=1:size(xcalib,2)
        for k=1:2
            error(k,i,j)=abs(fxcalib_ci(k,i,j)-xcalib(i,j));
        end
        xcalib_error(i,j)=max(error(:,i,j));
    end
end
comIN_square=comIN.^2; %Square input matrix
xcalib_error_square=xcalib_error.^2; %square error matrix
coeff_error=sqrt((xcalib_error_square'*comIN_square)'); %Matrix for error in each calibration point from uncertainty in coefficients(+/-)
