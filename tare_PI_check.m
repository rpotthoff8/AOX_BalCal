tare_load_approx=aprxINvalid(s_1stV,:);
tare_diff=abs(tare_load_approx-taresvalid);
tare_PI=loadPI_valid(s_1stV,:);
tare_in_PI=tare_diff<=tare_PI;
num_in_PI=sum(sum(tare_in_PI));
valid_percent_in_PI=num_in_PI/numel(tare_in_PI)

tare_load_approx=aprxINminGZ(s_1stV,:);
tare_diff=abs(tare_load_approx-tares);
tare_PI=y_hat_PI_comb(s_1stV,:);
tare_in_PI=tare_diff<=tare_PI;
num_in_PI=sum(sum(tare_in_PI));
calib_percent_in_PI=num_in_PI/numel(tare_in_PI)