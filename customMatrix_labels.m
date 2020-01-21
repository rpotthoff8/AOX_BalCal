function [leftColumn, topRow]=customMatrix_labels(loadlist,voltagelist,voltdimFlag,loaddimFlag,model,combined_terms,numRBF)
%Function generates text labels for terms used in model.  Labels are used
%in annotating ANOVA outputs for properties on each coefficient

%INPUTS:
%  loadlist = Labels for each load channel
%  voltagelist = Labels for each voltage channel
%  dimFlag = Number of data channels
%  model = Model type selected (Full, Truncated, Linear, Custom)
%  combined_terms = Character array for which term type is algebraicly combined (voltage or loads)
%  numRBF = Number of RBFs placed (if any)

%OUTPUTS:
%  leftColumn = Labels for combined terms
%  topRow = Labels for target terms

%Determine what terms should be combined: voltages or loads;
if strcmp(combined_terms,{'voltages'})==1
    toplist=loadlist;
    leftlist=voltagelist;
    %Variable labels are loads
    topRow=toplist(1:loaddimFlag)';
else
    toplist=voltagelist;
    leftlist=loadlist;
    %Variable labels are loads
    topRow=toplist(1:voltdimFlag)';
end



%Initialize counter and empty variables
count5=1;
block1=cell(voltdimFlag,1);
block2=cell(voltdimFlag,1);
block3=cell(voltdimFlag,1);
block4=cell(voltdimFlag,1);
block5=cell(((voltdimFlag-1)*voltdimFlag)/2,1);
block6=cell(((voltdimFlag-1)*voltdimFlag)/2,1);
block7=cell(((voltdimFlag-1)*voltdimFlag)/2,1);
block8=cell(((voltdimFlag-1)*voltdimFlag)/2,1);
block9=cell(voltdimFlag,1);
block10=cell(voltdimFlag,1);

%write text for variable names and combinations
for i=1:voltdimFlag
    block1(i)=leftlist(i);
    block2(i)=strcat('|',leftlist(i),'|');
    block3(i)=strcat(leftlist(i),'*',leftlist(i));
    block4(i)=strcat(leftlist(i),'*|',leftlist(i),'|');
    
    for j=i+1:voltdimFlag
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

if nargin>=7
    if numRBF>0
        channel=repmat([1:loaddimFlag]',numRBF,1);
        rbf_leftColumn=cellstr(strcat(reshape(toplist(channel),numel(channel),1), repmat(' RBF ',numRBF*loaddimFlag,1), num2str(repelem([1:numRBF]',loaddimFlag,1))));
        leftColumn=[leftColumn;rbf_leftColumn];
    end
end

end
