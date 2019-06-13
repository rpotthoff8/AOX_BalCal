Res_dif=y_hat_PI_comb-abs(targetRes);
Res_within=Res_dif>=0;
Percent_within=(sum(Res_within))./(size(targetRes,1));

%% Plotting
figure(1)
  sz=8;
 for i =1:6
     subplot(2,3,i)
 scatter(targetMatrix(:,i),Res_dif(:,i),sz,'filled')
 hold on
  title(strcat('Channel: ', string(i),'; Residuals within PI: ',string(round(100* Percent_within(i),1)),'%'))
ylabel('PI-abs(Residual) (lbs)')
xlabel('Target Load (lbs)')
sgtitle('Comparison in Prediction Interval and Residual')
grid on
grid minor
 hold off
 end
 hL= legend('(y hat PI)-abs(targetRes); Positive value indicates residual is within 95% PI');
 % Programatically move the Legend
newPosition = [0.42 0.47 0.2 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);