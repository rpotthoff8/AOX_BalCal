max_per_out=cell(1,16);
if FLAGS.model==1
    set={'Full'};
elseif FLAGS.model==2
    set={'Truncated'};
elseif FLAGS.model==3
    set={'Linear'};
else
    set={'Custom'};
end
max_per_out=[set,numBasis,eps_min,eps_max,max_mult,num2cell(std(targetRes2_resolve)),num2cell(std(targetRes2valid))]