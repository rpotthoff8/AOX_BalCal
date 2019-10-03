function [leftColumn, topRow]=customMatrix_labels(loadlist,voltagelist,dimFlag,model,combined_terms)

%Determine what terms should be combined: voltages or loads;
if strcmp(combined_terms,{'voltages'})==1
    toplist=loadlist;
    leftlist=voltagelist;
else
    toplist=voltagelist;
    leftlist=loadlist;
end

%Variable labels are voltages
topRow=toplist(1:dimFlag)';

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
    block1(i)=leftlist(i);
    block2(i)=strcat('|',leftlist(i),'|');
    block3(i)=strcat(leftlist(i),'*',leftlist(i));
    block4(i)=strcat(leftlist(i),'*|',leftlist(i),'|');
    
    for j=i+1:dimFlag
        block5(count5)=strcat(leftlist(i),'*',leftlist(j));
        block6(count5)=strcat('|',leftlist(i),'*',leftlist(j),'|');
        block7(count5)=strcat(leftlist(i),'*|',leftlist(j),'|');
        block8(count5)=strcat('|',leftlist(i),'|*',leftlist(j));
        count5=count5+1;
    end
    block9(i)=strcat(leftlist(i),'*',leftlist(i),'*',leftlist(i));
    block10(i)=strcat('|',leftlist(i),'*',leftlist(i),'*',leftlist(i),'|');
end

%Select Terms based on model type selected
if model==3
    leftColumn =block1;
elseif model==2
    leftColumn=[block1;block3;block5];
else
    leftColumn=[block1;block2;block3;block4;block5;block6;block7;block8;block9;block10];
end

end
