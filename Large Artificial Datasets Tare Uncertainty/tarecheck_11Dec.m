load('artificial_tares.mat')
count=0;
for i=1:size(artificial_tares,1)
    for j=1:size(artificial_tares,2)
        if artificial_tares(i,j)>=fzap_ci(1,i,j) && artificial_tares(i,j)<=fzap_ci(2,i,j)
            in_range(i,j)=1;
        else
            in_range(i,j)=0;
            disp(strcat('Series: ',string(i),'; Channel: ',string(j), '; actual tare: ',string(artificial_tares(i,j)),'; tare: ',string(taresALGB(i,j)), '; bound: ', string(fzap_ci(1,i,j)), ":",string(fzap_ci(2,i,j))))
            count=count+1;
        end
    end
end
disp(strcat(string(count),'/',string(numel(artificial_tares)), ' Actual Tares outside of uncertainty range.'))