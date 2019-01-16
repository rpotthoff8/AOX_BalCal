function correlationPlot(ind_var,dep_var,ind_label,dep_label)
n_dim = size(ind_var,2);

ccmat = zeros(n_dim,n_dim);
for j = 1:n_dim
    for k = 1:n_dim
        it_in = (j-1)*n_dim + k;
        cc = corrcoef(ind_var(:,j),dep_var(:,k));
        subplot(n_dim,n_dim,it_in)
        plot(ind_var(:,j),dep_var(:,k),'.')
        title(cc(2))
        xlabel(ind_label{j}); ylabel(dep_label{k});
        set(gca,'xtick',[],'ytick',[])
        ccmat(j,k) = cc(2);
    end
end