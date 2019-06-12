checkit=aprxINminGZ-targetMatrix0;
[taretal_mean,tarestdevsub_mean] = meantare(series,checkit);

% SOLVE TARES BY TAKING THE MEAN
[s,~,s_id] = unique(series0);
nseries = length(s);
for i = 1:nseries
    check_mean(i,:) = mean(checkit(s_id==i,:));
    boot_mean_CI_result(:,:,i)=bootci(20000,{@mean, checkit(s_id==i,:)});
    boot_mean_error(:,:,i)=abs(boot_mean_CI_result(:,:,i)-check_mean(i,:));
    boot_mean_CI(i,:)=max(boot_mean_error(:,:,i));
    check_std(i,:) = std(checkit(s_id==i,:));
    series_count(i)=sum(s_id==i);
    check_standard_error(i,:)=check_std(i,:)./sqrt(series_count(i));
    t_stat(i)=tinv(0.975,series_count(i)-1);
    check_mean_CI(i,:)=check_standard_error(i,:).*t_stat(i);
end

%% scatter plotting
figure(1)
 sz=8;
 for i =1:6
     subplot(2,3,i)
scatter(1:max(series0),beta_CI_comb(size(beta_CI_comb,1)-max(series0)+1:size(beta_CI_comb,1),i),sz,'filled')
hold on
scatter(1:max(series0),xcalib_error(size(xcalib_error,1)-max(series0)+1:size(xcalib_error,1),i),sz,'filled')
scatter(1:max(series0),check_mean_CI(:,i),sz,'filled')
title(strcat('Channel: ', string(i)))
ylabel('Tare 95% CI (+/- lbs)')
xlabel('Series')
legend('ANOVA (Beta CI)','Bootstrap','CI of checkit mean')
sgtitle('Comparison in Tare Uncertainty')
grid on
hold off
 end

 figure(2)
  for i =1:6
     subplot(2,3,i)
scatter(1:max(series0),beta_CI_comb(size(beta_CI_comb,1)-max(series0)+1:size(beta_CI_comb,1),i)-check_mean_CI(:,i),sz,'filled')
title(strcat('Channel: ', string(i)))
ylabel('Difference Tare 95% CI (+/- lbs)')
xlabel('Series')
% legend('ANOVA (Beta CI)','Bootstrap','CI of checkit mean')
sgtitle('(beta CI)-(CI of checkit mean) Tare Uncertainty Difference')
grid on
hold off
  end
 
  %% Bar plotting
  figure(3)
  sz=8;
 for i =1:6
     subplot(2,3,i)
 bar([beta_CI_comb(size(beta_CI_comb,1)-max(series0)+1:size(beta_CI_comb,1),i) xcalib_error(size(xcalib_error,1)-max(series0)+1:size(xcalib_error,1),i) check_mean_CI(:,i) boot_mean_CI(:,i)  ])
 hold on
 title(strcat('Channel: ', string(i)))
ylabel('Tare 95% CI (+/- lbs)')
xlabel('Series')
sgtitle('Comparison in Tare Uncertainty')
grid on
 hold off
 end
 hL= legend('ANOVA (Beta CI)','Bootstrap (Beta CI)','CI of checkit mean','Bootstrap of checkit mean');
 % Programatically move the Legend
newPosition = [0.42 0.47 0.2 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'NumColumns',4);

 figure(4)
  sz=8;
 for i =1:6
     subplot(2,3,i)
 bar([beta_CI_comb(size(beta_CI_comb,1)-max(series0)+1:size(beta_CI_comb,1),i) check_mean_CI(:,i)])
 hold on
 title(strcat('Channel: ', string(i)))
ylabel('Tare 95% CI (+/- lbs)')
xlabel('Series')
sgtitle('Comparison in Tare Uncertainty')
grid on
 hold off
 end
 hL= legend('ANOVA (Beta CI)','CI of checkit mean');
 % Programatically move the Legend
newPosition = [0.42 0.47 0.2 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(hL,'NumColumns',4);
