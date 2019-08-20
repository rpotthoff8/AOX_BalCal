function [x_opt,x_val] = max_218(func, bounds, options)
%John Potthoff
%Speedy Bee II
% Artificial Bee Colony Algorithm
% Modified for speed (hopefully)
if nargin==3
    if isfield(options,'maxIter')
        itermax = options.maxIter ;
    else
        itermax=10^8;
    end
else
    itermax=10^8;
end

fminfunc = @(x) 0-func(x) ;
nvar=size(bounds,1);

npop = 45; %Bee Population size

xmin = bounds(:,1);
xmax = bounds(:,2);
a = 1;

lim = round(nvar*npop); %Limit for abandoning point
globest_cost = inf;
globest_loc=zeros(nvar,1);

cost=zeros(npop,1);
% Scout Bees initiate Food Source
loc=  (xmax-xmin).*rand(nvar,npop)+xmin;
for i = 1:npop
    cost(i) = fminfunc(loc(:,i));
end

[minCost,minI]=min(cost);
if minCost < globest_cost
    globest_cost = minCost;
    globest_loc=loc(:,minI);
end

C = zeros(npop,1);
abandonTrack=zeros(npop,1);


Kmatrix=repmat(1:npop,npop,1);
Kmatrix = Kmatrix(~eye(size(Kmatrix)));
Kmatrix=flipud(reshape(Kmatrix,npop,npop-1));

actualIter=0;
outerLoopCount=0;
printerCount=round(1000/(npop*3));
while actualIter<=itermax
    outerLoopCount=outerLoopCount+1;
    
    % Employed Bees
    k=Kmatrix(sub2ind(size(Kmatrix),(1:size(Kmatrix,1))',randi(npop-1,npop,1)));
    phi = a*unifrnd(-1,1,nvar,npop);
    new_loc = min(max(loc + phi.*(loc-loc(:,k)), xmin),xmax);
    
    new_cost=zeros(npop,1);
    for i = 1:npop
        new_cost(i) = fminfunc(new_loc(:,i));
    end
    
    %Replace entries where new cost is lower, otherwise add 1 to counter
    Replace=new_cost<cost;
    cost(Replace)=new_cost(Replace);
    loc(:,Replace)=new_loc(:,Replace);
    C=C+(1-Replace);
    
    
    % Onlooker Bees
    F = zeros(npop,1);
    F(cost>=0)=1./(1+cost(cost>=0));
    F(cost<0)= 1+abs(cost(cost<0));
    
    P = F/sum(F);
    
    randmatrix=rand(npop,1);
    sumP=cumsum(P);
    randiVec=randi(npop-1,npop,1);
    for j = 1:npop
        i=find(randmatrix(j)<=sumP,1,'first');
        K = [1:i-1 i+1:npop];
        k(j) = K(randiVec(j));
    end
    
    phi = a*unifrnd(-1,1,nvar,npop);
    new_loc = min(max(loc + phi.*(loc-loc(:,k)), xmin),xmax);
    new_cost=zeros(npop,1);
    for i = 1:npop
        new_cost(i) = fminfunc(new_loc(:,i));
    end
    
    %Replace entries where new cost is lower, otherwise add 1 to counter
    Replace=new_cost<cost;
    cost(Replace)=new_cost(Replace);
    loc(:,Replace)=new_loc(:,Replace);
    C=C+(1-Replace);
    
    
    % Scout Bees
    Replace=C>=lim;
    ReplaceI=find(Replace);
    loc(:,ReplaceI)=(xmax-xmin).*rand(nvar,numel(ReplaceI))+xmin;
    C(ReplaceI)=0;
    abandonTrack(ReplaceI)=1;
    for i = 1:numel(ReplaceI)
        cost(ReplaceI(i)) = fminfunc(loc(:,ReplaceI(i)));
    end
    actualIter=actualIter+2*npop+numel(ReplaceI); %Actual iterations is defined as the number of function evaluations performed
    
    
    [roundMin,minI]=min(cost);
    if roundMin < globest_cost
        globest_cost = roundMin;
        globest_loc=loc(:,minI);
        abandonTrack=zeros(npop,1);
    end
    
    if sum(abandonTrack)>=.75*npop
        break
    end
    if nargin==3
        if isfield(options,'display') && rem(outerLoopCount,printerCount)==0
            if options.display == 1
                fprintf(['Iteration ' num2str(actualIter) ': the optimal value is ' num2str(-globest_cost) '.\n'])
            end
        end
    end
end

x_opt=globest_loc;
x_val=globest_cost;
x_val = 0-x_val ;
end
