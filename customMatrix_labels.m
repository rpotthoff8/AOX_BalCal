function [recTable]=customMatrix_labels(loadlist,voltagelist,dimFlag,RECOMM_ALG_EQN,FLAGS)
%Variable labels are voltages
topRow=loadlist(1:dimFlag);

%Initialize counter and empty variables
count5=1;
block1=cell(dimFlag,1);
block2=cell(dimFlag,1);
block3=cell(dimFlag,1);
block4=cell(dimFlag,1);
block5=cell(((dimFlag-1)*dimFlag)/2,1);
block6=cell(((dimFlag-1)*dimFlag)/2,1);
block7=cell(((dimFlag-1)*dimFlag)/2,1);
block8=cell(((dimFlag-1)*dimFlag)/2,1);
block9=cell(dimFlag,1);
block10=cell(dimFlag,1);

%write text for variable names and combinations
for i=1:dimFlag
    block1(i)=voltagelist(i);
    block2(i)=strcat('|',voltagelist(i),'|');
    block3(i)=strcat(voltagelist(i),'^2');
    block4(i)=strcat(voltagelist(i),'*|',voltagelist(i),'|');
    
    for j=i+1:dimFlag
        block5(count5)=strcat(voltagelist(i),'*',voltagelist(j));
        block6(count5)=strcat('|',voltagelist(i),'*',voltagelist(j),'|');
        block7(count5)=strcat(voltagelist(i),'*|',voltagelist(j),'|');
        block8(count5)=strcat('|',voltagelist(i),'|*',voltagelist(j));
        count5=count5+1;
    end
    block9(i)=strcat(voltagelist(i),'^3');
    block10(i)=strcat('|',voltagelist(i),'^3|');
end

%Select Terms based on model type selected
if FLAGS.model==3
    leftColumn =block1;
elseif FLAGS.model==2
    leftColumn=[block1;block3;block5];
else
    leftColumn=[block1;block2;block3;block4;block5;block6;block7;block8;block9;block10];
end

%Combine in table
recTable=table(RECOMM_ALG_EQN,'VariableNames',topRow,'RowNames',leftColumn);

end
