clear all;
%Suppose we have a dataset perfectly described by the following
%relationship:
%   F = a*V+b*V^2+d*V^3
%The true coefficient values are:
a=2;
b=0.5;
d=0.3;
coeff_true=[a;b;d];

%   This dataset is collected in 4 series of 20 datapoints each:
series=[repmat(1,20,1);repmat(2,20,1);repmat(3,20,1); repmat(4,20,1)];
[~,s_1st,~] = unique(series);

numpts=length(series); %total number of datapoints

% V values are distributed across all 4 series:
V_in=10*rand(numpts,1);

%The input variables can be assembled in a matrix:
comIN=[V_in,V_in.^2,V_in.^3];

%The true F values can be calculated:
F_true=comIN*coeff_true;

%Plot this true relationship
figure(1)
subplot(3,2,1);
scatter(V_in,F_true,20,series,'filled');
xlabel('V');
ylabel('F');
title('True relationship V>F');

disp('See Figure 1 for true relationship between V and F (global load).');
disp('Press ENTER to continue.');
fprintf('\n');
pause();


%There are tare loads present, however so the reported F value (F') is:
%   F'=F-T_k 
% where T_k is the series tare load

%The true tare values are:
tares_true=[20;-2;-17;8]; 
F_prime=F_true-tares_true(series); %Subtract tares from true global F for F'
%Plot the reported data:
subplot(3,2,2);
scatter(V_in,F_prime,20,series,'filled');
xlabel('V');
ylabel('F prime');
title('V>F prime');

disp('See Figure 1 for true relationship between V and reported load, F prime (global load - tares).');
disp('Press ENTER to continue.');
fprintf('\n');
pause();

%Our goal is to solve for the true coefficients for the relationship V>F.
%We only have access to V and F' to do so, however:

%% Our first approach is regressing for the coefficients, then taking the
%mean of the difference in each series as the tare load:

% Normalize the data for a better conditioned matrix
scale1 = max(abs(comIN));
scale1(scale1==0)=1; %To avoid NaN for channels where RBFs have self-terminated
comIN_1 = comIN./scale1;

xcalib_1 = comIN_1\F_prime; %Solve for coefficients:
xcalib_1 = xcalib_1./scale1'; % De-normalize the coefficients to be used with raw data
coeff_1 = xcalib_1;

%Calculate our approximation without tare correcting:
aprxIN1 = comIN*coeff_1;

%Find residual between F' and approximation:
checkit1 = aprxIN1-F_prime;
taresAllPoints1 = meantare(series,checkit1);
tares1=taresAllPoints1(s_1st);
%Tare subtracted approximation to compare to F_prime
aprxIN1_mTares=aprxIN1-taresAllPoints1;

%Calculate true RMS and tare corrected RMS
RMS1_true=sqrt(mean((aprxIN1-F_true).^2));
RMS1_tareC=sqrt(mean((aprxIN1_mTares-F_prime).^2));
%Calculate Tare RMS
tare_RMS1=sqrt(mean((tares1-tares_true).^2));

%Plot the approximation:
subplot(3,2,3);
scatter(V_in,aprxIN1,20,series,'filled');
xlabel('V');
ylabel('F approx');
title('True load (F) Approximation 1');

subplot(3,2,4);
scatter(V_in,aprxIN1_mTares,20,series,'filled');
xlabel('V');
ylabel('F approx');
title('Tare subracted (F prime) Approximation 1');

disp('Approach 1 calculates coefficients without intercepts, then calculates tares as mean difference between approximation and F prime in each series')
fprintf('RMS between global approximation and true global load: '); fprintf(num2str(RMS1_true)); fprintf('\n');
fprintf('RMS between tare corrected approximation and F_prime: '); fprintf(num2str(RMS1_tareC)); fprintf('\n');
fprintf('RMS between calculated and true tares: '); fprintf(num2str(tare_RMS1)); fprintf('\n');
disp('Press ENTER to continue.');
fprintf('\n');
pause();
%% Our second approach is regressing for the coefficients including a specific intercept term in each series

%We will generate series specific intercepts for each series
comIN2_intercepts=zeros(numpts,max(series));
for i=1:max(series)
    comIN2_intercepts(series==i,i)=1;
end
comIN_2=[comIN,comIN2_intercepts]; %Our predictor variables for regression include the series specific tare intercepts

% Normalize the data for a better conditioned matrix
scale2 = max(abs(comIN_2));
scale2(scale2==0)=1; %To avoid NaN for channels where RBFs have self-terminated
comIN_2 = comIN_2./scale2;

xcalib_2 = comIN_2\F_prime; %Solve for coefficients:
xcalib_2 = xcalib_2./scale2'; % De-normalize the coefficients to be used with raw data

coeff_2 = xcalib_2(1:length(coeff_true)); %Extract coefficients
tares2 = -xcalib_2(length(coeff_true)+1:end); % Tares were solved for with intercept terms
taresAllPoints2 = tares2(series); %Expand out to length of input V

%Calculate our approximation without tare correcting:
aprxIN2 = comIN*coeff_2;

%Tare subtracted approximation to compare to F_prime
aprxIN2_mTares=aprxIN2-taresAllPoints2;

%Calculate true RMS and tare corrected RMS
RMS2_true=sqrt(mean((aprxIN2-F_true).^2));
RMS2_tareC=sqrt(mean((aprxIN2_mTares-F_prime).^2));
%Calculate Tare RMS
tare_RMS2=sqrt(mean((tares2-tares_true).^2));

%Plot the approximation:
subplot(3,2,5);
scatter(V_in,aprxIN2,20,series,'filled');
xlabel('V');
ylabel('F approx');
title('True load (F) Approximation 2');

subplot(3,2,6);
scatter(V_in,aprxIN2_mTares,20,series,'filled');
xlabel('V');
ylabel('F approx');
title('Tare subracted (F prime) Approximation 2');

disp('Approach 2 calculates coefficients with series specific intercepts.  This provides the tares along with coefficients')
fprintf('RMS between global approximation and true global load: '); fprintf(num2str(RMS2_true)); fprintf('\n');
fprintf('RMS between tare corrected approximation and F_prime: '); fprintf(num2str(RMS2_tareC)); fprintf('\n');
fprintf('RMS between calculated and true tares: '); fprintf(num2str(tare_RMS2)); fprintf('\n');
disp('Press ENTER to continue.');
fprintf('\n');
pause()
%% Our third approach is from the Tare formulation provided by Dr. Meade on 4 FEB 20
%% Approach 3A: Calculate coefficients as in approach 1, calculate tares from mean inputs (V) in each series and coefficients

% Normalize the data for a better conditioned matrix
scale3A = max(abs(comIN));
scale3A(scale3A==0)=1; %To avoid NaN for channels where RBFs have self-terminated
comIN_3A = comIN./scale3A;

xcalib_3A = comIN_3A\F_prime; %Solve for coefficients:
xcalib_3A = xcalib_3A./scale3A'; % De-normalize the coefficients to be used with raw data
coeff_3A = xcalib_3A;

%Calculate our approximation without tare correcting:
aprxIN3A = comIN*coeff_3A;

%Determine mean value for each input variable in series for calculating Tares
mean_Term_series=zeros(max(series),size(comIN,2));
comIN_3A_tare=zeros(max(series),size(comIN,2));
for i=1:max(series)
    mean_Term_series(i,:)=mean(comIN(series==i,:),1); %Determine mean for each input variable in series
    comIN_3A_tare(i,:)=mean_Term_series(i,:);
end

%Calculate tares using coefficients and mean input values
tares3A=comIN_3A_tare*coeff_3A;

%Tare subtracted approximation to compare to F_prime
aprxIN3A_mTares=aprxIN3A-tares3A(series);

%Calculate true RMS and tare corrected RMS
RMS3A_true=sqrt(mean((aprxIN3A-F_true).^2));
RMS3A_tareC=sqrt(mean((aprxIN3A_mTares-F_prime).^2));
%Calculate Tare RMS
tare_RMS3A=sqrt(mean((tares3A-tares_true).^2));

%Plot the approximation:
figure(2)
subplot(3,2,1);
scatter(V_in,aprxIN3A,20,series,'filled');
xlabel('V');
ylabel('F approx');
title('True load (F) Approximation 3A');

subplot(3,2,2);
scatter(V_in,aprxIN3A_mTares,20,series,'filled');
xlabel('V');
ylabel('F approx');
title('Tare subracted (F prime) Approximation 3A');

disp('Approach 3A calculates coefficients without intercepts as in Approach 1, then calculates tares using mean input in each series and coefficients')
fprintf('RMS between global approximation and true global load: '); fprintf(num2str(RMS3A_true)); fprintf('\n');
fprintf('RMS between tare corrected approximation and F_prime: '); fprintf(num2str(RMS3A_tareC)); fprintf('\n');
fprintf('RMS between calculated and true tares: '); fprintf(num2str(tare_RMS3A)); fprintf('\n');
disp('Press ENTER to continue.');
fprintf('\n');
pause();

%% Approach 3B: Calculate coefficients using formulation subtracting mean of input in each series (zero mean each series), calculate tares from mean inputs in each series and coefficients
%Determine mean value for each input variable in series for calculating Tares
mean_Term_series=zeros(max(series),size(comIN,2));
comIN_3B_tare=zeros(max(series),size(comIN,2));
for i=1:max(series)
    mean_Term_series(i,:)=mean(comIN(series==i,:),1); %Determine mean for each input variable in series
    comIN_3B_tare(i,:)=mean_Term_series(i,:);
end
comIN_3B=comIN-mean_Term_series(series,:); %Calculate coefficients using comIN corrected with mean voltages

% Normalize the data for a better conditioned matrix
scale3B = max(abs(comIN_3B));
scale3B(scale3B==0)=1; %To avoid NaN
comIN_3B = comIN_3B./scale3B;

xcalib_3B = comIN_3B\F_prime; %Solve for coefficients:
xcalib_3B = xcalib_3B./scale3B'; % De-normalize the coefficients to be used with raw data
coeff_3B = xcalib_3B;

%Calculate our approximation without tare correcting:
aprxIN3B = comIN*coeff_3B;

%Calculate tares using coefficients and mean input values
tares3B=comIN_3B_tare*coeff_3B;

%Tare subtracted approximation to compare to F_prime
aprxIN3B_mTares=aprxIN3B-tares3B(series);

%Calculate true RMS and tare corrected RMS
RMS3B_true=sqrt(mean((aprxIN3B-F_true).^2));
RMS3B_tareC=sqrt(mean((aprxIN3B_mTares-F_prime).^2));
%Calculate Tare RMS
tare_RMS3B=sqrt(mean((tares3B-tares_true).^2));

%Plot the approximation:
subplot(3,2,3);
scatter(V_in,aprxIN3B,20,series,'filled');
xlabel('V');
ylabel('F approx');
title('True load (F) Approximation 3B');

subplot(3,2,4);
scatter(V_in,aprxIN3B_mTares,20,series,'filled');
xlabel('V');
ylabel('F approx');
title('Tare subracted (F prime) Approximation 3B');

disp('Approach 3B calculates coefficients without intercepts using comIN corrected with mean inputs in each series, then calculates tares using mean input in each series and coefficients')
fprintf('RMS between global approximation and true global load: '); fprintf(num2str(RMS3B_true)); fprintf('\n');
fprintf('RMS between tare corrected approximation and F_prime: '); fprintf(num2str(RMS3B_tareC)); fprintf('\n');
fprintf('RMS between calculated and true tares: '); fprintf(num2str(tare_RMS3B)); fprintf('\n');
disp('Press ENTER to continue.');
fprintf('\n');