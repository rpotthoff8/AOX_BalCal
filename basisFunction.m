%For multi Dimensions
function phi = basisFunction(x,xc,w)
n = size(xc,2);
numdata = size(x,1);
bee = zeros(n,1);
stuffer = ones(numdata,1);

for i = 1:n
    zee = exp(  (  (x-xc(:,i)').^2 )*w(:,i));
    bee(i) = dot(stuffer,zee)/(numdata+0.0);
%     phi(:,i) = exp(  (  ([x]-xc(:,i)').^2 )*w(:,i))  - bee(i)*stuffer;
    phi(:,i) = exp(  (  ([x]-xc(:,i)').^2 )*w(:,i));

end
end
