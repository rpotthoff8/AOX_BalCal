%Function plots correlation between variables
function correlationPlot(ind_var,dep_var,ind_label,dep_label)
%Function plots correlation between variables

%INPUTS:
%  ind_var  = Matrix of Independent variable
%  dep_var  =  Matrix of Dependent variable
%  ind_label = Label for Independent variable
%  dep_label  =  Label for Dependent variable

%OUTPUTS:
%  comIN_RBF  =  Input for load approximation, to be multiplied by
%  coefficients with algebraic comIN to determine load value


ind_n_dim = size(ind_var,2); %Dimesion of independent data (# channels)
dep_n_dim= size(dep_var,2); %Dimension of dependent variable (# channels)

%Loop through variables and generate subplots for correlation between
%independent and dependent variables
ccmat = zeros(ind_n_dim,dep_n_dim);
for j = 1:ind_n_dim 
    for k = 1:dep_n_dim
        it_in = (j-1)*dep_n_dim + k;
        cc = corrcoef(ind_var(:,j),dep_var(:,k));
        subplot(ind_n_dim,dep_n_dim,it_in)
        plot(ind_var(:,j),dep_var(:,k),'.')
        title(cc(2), 'Interpreter', 'none')
        xlabel(ind_label{j}, 'Interpreter', 'none'); ylabel(dep_label{k}, 'Interpreter', 'none');
        set(gca,'xtick',[],'ytick',[])
        ccmat(j,k) = cc(2);
    end
end