mu=-16.6;
sigma=2.2;
kurt=3.3;
m=5000;
n=1;
skew=(-1:.01:1);
for i=1:numel(skew)
r=pearsrnd(mu,sigma,skew(i),kurt,m,n);
meanr=r-mean(r);
medianr=r-median(r);
meanrms(i)=rms(meanr);
medianrms(i)=rms(medianr);
end
diff=medianrms-meanrms;
plot(skew,diff)
xlabel('Skewness')
ylabel('(Median rms)-(Mean rms)')
