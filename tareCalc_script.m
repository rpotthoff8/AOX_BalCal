%Run with "No intercept" option selected
%Insert Breakpoint at Line 408:
%Run Following:
mean_Term_series=zeros(nseries0,size(comIN0,2));
comIN_tare=zeros(nseries0,size(comIN0,2));
numpts_series=zeros(nseries0,1);
for i=1:nseries0
    mean_Term_series(i,:)=mean(comIN0(series0==i,:),1);
    numpts_series(i)=sum(series0==i);
%     comIN_tare(i,:)=mean_Term_series(i,:)/numpts_series(i);
    comIN_tare(i,:)=mean_Term_series(i,:);

end
comIN_tare_tal=comIN_tare(series0,:);
comIN_TareC=comIN-comIN_tare_tal;

%Insert Breakpoint at Line 488:
%Run Following:
tares=comIN_tare*xcalib;
[tares2AllPointsvalid,taretalstdvalid] = meantare(series0,targetRes);% SOLVE FOR TARES BY TAKING THE MEAN

targetRes = (targetMatrix0+tares(series0,:))-aprxIN;

