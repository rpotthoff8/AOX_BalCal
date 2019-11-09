eps_out=cell(1,16);
if FLAGS.model==1
    set={'Full'};
elseif FLAGS.model==2
    set={'Truncated'};
elseif FLAGS.model==3
    set={'Linear'};
else
    set={'Custom'};
end
eps_out=[set,numBasis,eps_min,eps_max,num2cell(std(targetRes2)),num2cell(std(targetRes2valid))]