max_per_out=cell(16,17);
if FLAGS.model==1
    set={'Full'};
elseif FLAGS.model==2
    set={'Truncated'};
elseif FLAGS.model==3
    set={'Linear'};
else
    set={'Custom'};
end
max_mult=1;
run_info=[set,numBasis,eps_min,eps_max,max_mult];
max_per_out(:,1:5)=repmat(run_info,16,1);
extract_num=[1;5;10;25;35;50;100;250;400;500;750;1000;1250;1500;1750;2000];
max_per_out(:,2)=num2cell(extract_num);
max_per_out(:,6:17)=num2cell([resStdHist(extract_num,:),resStdHistvalid(extract_num,:)]);
% max_per_out=[set,numBasis,eps_min,eps_max,max_mult,num2cell(std(targetRes2)),num2cell(std(targetRes2valid))]
