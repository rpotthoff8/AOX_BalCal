%Run with "No intercept" option selected
%Insert Breakpoint at Line 408:
%Run Following:
mean_Term_series=zeros(nseries0,size(comIN0,2));
comIN_tare=zeros(nseries0,size(comIN0,2));
for i=1:nseries0
    mean_Term_series(i,:)=mean(comIN0(series0==i,:),1);
    comIN_tare(i,:)=mean_Term_series(i,:);

end
comIN_TareC=comIN-mean_Term_series(series0);

%Insert Breakpoint at Line 488:
%Run Following:
tares=comIN_tare*xcalib;
[tares2AllPointsvalid,taretalstdvalid] = meantare(series0,targetRes);% SOLVE FOR TARES BY TAKING THE MEAN

targetRes = (targetMatrix0+tares(series0,:))-aprxIN;

