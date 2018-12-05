count=0;
for i=1:size(taresALGB,1)
    for j=1:size(taresALGB,2)
        if zapvalid(i,j)>=zap_ci(2*i-1,j) && zapvalid(i,j)<=zap_ci(2*i,j)
            in_range(i,j)=1;
        else
            in_range(i,j)=0;
            disp(strcat('Series: ',string(i),'; Channel: ',string(j), '; validation tare: ',string(zapvalid(i,j)),'; tare: ',string(taresALGB(i,j)), '; bound: ', string(zap_ci(2*i-1,j)), ":",string(zap_ci(2*i,j))))
            count=count+1;
        end
    end
end
disp(strcat(string(count),'/',string(numel(taresALGB)), ' Validation Tares outside of uncertainty range.'))