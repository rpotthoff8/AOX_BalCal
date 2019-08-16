% Last Updated: 6/04/19

function [err,xc,w,c,err_hist,add,data]=wahba5_2D_DQ_Driver(numRBF,data,options)
%
%INPUT VARIABLES:
%numRBF: number of RBFs to be used
%data: data to use in function, includes following:
%       data.x: x positions, COLUMN vector, MUST BE PROVIDED
%       data.f0: y values for x positions, COLUMN vector, MUST BE PROVIDED
%       data.x_test: test x positions, COLUMN vector, OPTIONAL
%       data.f0a_test: test y values, COLUMN vector, OPTIONAL
%       data.noise_std: guess for standard deviation of noise, used for
%       self termination, Single value, OPTIONAL
%       data.f_true: actual y values (without noise), COLUMN vector, OPTIONAL

%options: options for function, includes the following:
%       data.plot: if plotter should be used: 1 for yes, 0 for no
%       (default is to not plot)
%       data.terminate_percent: if data.noise_std is provided, code will
%       terminate when this percentage of residuals are <= 2*noise_std,
%       OPTIONAL

%Variables for plotter function
method=options.label;

%Check if test data is provided, set 'test_data=1' if it is
if isfield(data,'x_test')
    test_data=1 ;
    [~, I_xOrder]=sortrows([data.x;data.x_test]); %Sort x in order, save indices of order
else
    test_data=0;
end

%Check if true data is provide, set 'true_data=1' if it is
if isfield(data,'f_true')
    true_data=1;
else
    true_data=0;
end

%Check if plotter option is provide, if not, default is to turn plotter off
if isfield(options,'plot')
    plot=options.plot;
    fignum=options.fignum;
else
    plot=0;
end

% %Check if termination percent is provided, if not default is to set to 0.92
% if isfield(options,'terminate_percent')
%     cutoff_percent=options.terminate_percent;
% else
%     cutoff_percent=0.92;
% end

%variables determined from x data provided
Ndata = size(data.x,1);
xmin = min(data.x);
xmax = max(data.x);
della = (xmax-xmin)/(Ndata-1);

%Set convergence limit
converge_limit = 0.001;
%converge_limit = 0.25;

% END initialize section

% Initial residual
g0=mean(data.f0);
add = g0; %add is constant that is added to phi (RBF approximation)
%Initial approximation (zero)
f_hat.train = zeros(size(data.x,1),1)+add;
R.train = data.f0 -( f_hat.train);

%Store starting error
err.train = sqrt(mean(R.train.^2));
err_hist.train = err.train;

if test_data==1 %If test data is provided, determine intitial residual/error for test data
    f_hat.test=zeros(size(data.x_test,1),1)+add;
    R.test = data.f0_test - f_hat.test;
    err.test=sqrt(mean(R.test.^2));
    err_hist.test = err.test;
    f_hat.full=[f_hat.train;f_hat.test];
    f_hat.full=f_hat.full(I_xOrder,:);
else
    f_hat.full=f_hat.train;
end

if true_data==1 %If true data is provided, determine actual residual/error
    R.true=data.f_true-f_hat.full;
    err.true=sqrt(mean(R.true.^2));
    err_hist.true = err.true;
end

if plot==1 %Call plotter function if turned on
    plotter_2D(data,f_hat.full,err_hist,0,R,method,fignum);
end

%optLoop outputs the most optimizing pair of center and width.
if isfield(data,'DQcoeff')
    DQcoeff=data.DQcoeff;
else
    [DQcoeff,data]=DQcoeff_solver(data);
    data.DQcoeff=DQcoeff;
end

[xc(:,1),w(:,1), kstar(1),data,Iplace(1)] = DQ_2D_vW_maxR_optLoop(data,R.train,DQcoeff,options);
options.Iplace=Iplace;
% eval(['[xc(:,1),w(:,1), kstar(1)] = ' char(options.opt_funct) '(data.x,R.train,DQcoeff,options);']) %Run optloop
% [xc(:,1),w(:,1), kstar(1)] = optLoop(data.x,R.train,DQcoeff,options);

phi.train = basisFunction(data.x,xc(:,1),w(:,1));
add=g0-mean(phi.train);

% In the first iteration, the coefficient is solved for analytically.
% c = dot((phi.train+add),R.train)/dot((phi.train+add),(phi.train+add));
c = dot((phi.train-mean(phi.train)),R.train)/dot((phi.train-mean(phi.train)),(phi.train-mean(phi.train)));


%The approximation is updated.
% f_hat.train =c*(phi.train+add);
f_hat.train =c*(phi.train-mean(phi.train))+g0;
R.train = data.f0 - f_hat.train;
err.train = sqrt(mean(R.train.^2));
%The error history and the current approximation are plotted together every
%iteration.
err_hist.train = [err_hist.train;err.train];

if test_data==1 %If test data provided, determine approximation/residuals for test datapoints
    phi.test = basisFunction(data.x_test,xc(:,1),w(:,1));
    %The approximation is updated.
%     f_hat.test =c*(phi.test+add);
    f_hat.test =c*(phi.test-mean(phi.train))+g0;

    R.test = data.f0_test - f_hat.test;
    err.test = sqrt(mean(R.test.^2));
    err_hist.test = [err_hist.test;err.test];
    f_hat.full=[f_hat.train;f_hat.test];
    f_hat.full=f_hat.full(I_xOrder,:);
else
    f_hat.full=f_hat.train;
end

if true_data==1 %If true data is provided, determine actual residual/error
    R.true=(data.f_true)-f_hat.full;
    err.true=sqrt(mean(R.true.^2));
    err_hist.true = [err_hist.true;err.true];
end

if plot==1 %Call plotter function if turned on
    plotter_2D(data,f_hat.full,err_hist,c,R,method,fignum)
end


%NOTE: For the time being, the code is run a predetermined number of
%iterations, the convergence criteria has not been defined.

for i=2:numRBF %iterations
    
    %optLoop is run again, finding the "optimal" center/width pair for the
    %next basis function.
    %NOTE: the optLoop merit function to find the optimal width is derived
    %under the assumption that the coefficient will be calculated
    %analytically, though this will not be the case starting with the
    %second iteration. In this sense, the center/width pair aren't truly
    %optimal
    
    % Optimization
    %     [xc(:,i),w(:,i), kstar(i)] = optLoop(data.x,R.train,DQcoeff,options);
    %     eval(['[xc(:,i),w(:,i), kstar(i)] = ' char(options.opt_funct) '(data.x,R.train,DQcoeff,options);']) %Run optloop
    [xc(:,i),w(:,i), kstar(i),data,Iplace(i)] = DQ_2D_vW_maxR_optLoop(data,R.train,DQcoeff,options);
    options.Iplace=Iplace;

    phi.train = basisFunction(data.x,xc,w);
    add=g0-mean(phi.train);
    
    %Starting witht the second iteration and until the end of the program,
    %the coefficients are all recalculated using an inverse matrix
    %calculation, to find the best set of coefficients in a least-squares
    %sense.
    c = lsqminnorm(phi.train-mean(phi.train),data.f0);
%     c = pinv(phi.train-mean(phi.train))*data.f0;
    %    c = phi\f;
    %The aproximation is updated.
    f_hat.train = (phi.train-mean(phi.train))*c+g0;
    R.train = (data.f0) - f_hat.train;
    err.train = sqrt(mean(R.train.^2));
    err_hist.train = [err_hist.train;err.train];
    
    if test_data==1 %If test data provided, determine approximation/residuals for test datapoints
        phi.test = basisFunction(data.x_test,xc,w);
        %The aproximation is updated.
        f_hat.test = (phi.test-mean(phi.train))*c+g0;
        R.test = (data.f0_test) - f_hat.test;
        err.test = sqrt(mean(R.test.^2));
        err_hist.test = [err_hist.test;err.test];
        f_hat.full=[f_hat.train;f_hat.test];
        f_hat.full=f_hat.full(I_xOrder,:);
        
        if i>10
            period_change(i)=min(err_hist.test(i-9:i))-min(err_hist.test(i-10));
        else
            period_change(i)=min(err_hist.test(2:i))-min(err_hist.test(1));
        end
        
    else
        f_hat.full=f_hat.train;
    end
    
    if true_data==1 %If true data is provided, determine actual residual/error
        R.true=(data.f_true)-f_hat.full;
        err.true=sqrt(mean(R.true.^2));
        err_hist.true = [err_hist.true;err.true];
    end
    
    if plot==1 %Call plotter function if turned on
        plotter_2D(data,f_hat.full,err_hist,c,R,method,fignum)
    end
    
    converg = -log(err.train/err_hist.train(end-1));
    convergence_rate = converg;
    
    %% Termination Criteria %%%%%
    %    if converg < converge_limit
    %        break
    %    end
    %     if isfield(data,'noise_std') %If guess for noise standard deviation is provided
    %         percent_in_noise.train=sum(abs(R.train)<=2*data.noise_std)/numel(R.train); %Calculate percentage of data points that are within 2*(noise std)
    %         if test_data==1
    %             percent_in_noise.test=sum(abs(R.test)<=2*data.noise_std)/numel(R.test); %Calculate percentage of data points that are within 2*(noise std)
    %         end
    %         if percent_in_noise.train>=cutoff_percent
    %             disp(strcat('Reached noise std cutoff, # RBF=',string(numel(xc))))
    %             break
    %         end
    %     end
    if test_data==1 && options.self_terminate==1
        if period_change(i)>0
            disp(strcat('Reached period change termination criteria, # RBF=',string(size(xc,2))))
            break
        end
    end

end

%Set coefficeints to zero for RBFs past min training error
if test_data==1 && options.best_test==1
    [~,minItest]=min(err_hist.test);
    if minItest<i
        w(:,minItest+1:end)=[];
        xc(:,minItest+1:end)=[];
        phi.train = basisFunction(data.x,xc,w);
        add=mean(data.f0-phi.train);
        
        %Starting witht the second iteration and until the end of the program,
        %the coefficients are all recalculated using an inverse matrix
        %calculation, to find the best set of coefficients in a least-squares
        %sense.
        c = pinv(phi.train+add)*data.f0;
        %    c = phi\f;
        %The aproximation is updated.
        f_hat.train = (phi.train+add)*c;
        R.train = (data.f0) - f_hat.train;
        err.train = sqrt(mean(R.train.^2));
        
        if test_data==1 %If test data provided, determine approximation/residuals for test datapoints
            phi.test = basisFunction(data.x_test,xc,w);
            %The aproximation is updated.
            f_hat.test = (phi.test+add)*c;
            R.test = (data.f0_test) - f_hat.test;
            err.test = sqrt(mean(R.test.^2));
%             err_hist.test = [err_hist.test;err.test];
            f_hat.full=[f_hat.train;f_hat.test];
            f_hat.full=f_hat.full(I_xOrder,:);
        else
            f_hat.full=f_hat.train;
        end
        
        if true_data==1 %If true data is provided, determine actual residual/error
            R.true=(data.f_true)-f_hat.full;
            err.true=sqrt(mean(R.true.^2));
        end
    end
end

end

%% subroutine functions past this point.
%
%************************


%*********************