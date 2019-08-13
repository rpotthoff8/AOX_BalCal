function err = meritFunction(x,R,xc,w)
if size(xc,1)~=size(w,1)
    first=num2cell(ones(1,numel(size(w))));
    w=repmat(w(first{:}),size(xc));
end
phi = basisFunction(x,xc,w);
add=mean(phi);
f_hat=phi-add;
err = dot(R,R)-dot(f_hat,R)^2/(dot(f_hat,f_hat));


% #####################
end