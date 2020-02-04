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

% Term order for full equation set:
% (INTERCEPT), F, |F|, F*F, F*|F|, F*G, |F*G|, F*|G|, |F|*G, F*F*F, |F*F*F|, F*G*G, F*G*H

%Initialize counter and empty variables
count5=1;
count13=1;
block1=cellstr('INTERCEPT');
block2=cell(voltdimFlag,1);
block3=cell(voltdimFlag,1);
block4=cell(voltdimFlag,1);
block5=cell(voltdimFlag,1);
block6=cell(((voltdimFlag-1)*voltdimFlag)/2,1);
block7=cell(((voltdimFlag-1)*voltdimFlag)/2,1);
block8=cell(((voltdimFlag-1)*voltdimFlag)/2,1);
block9=cell(((voltdimFlag-1)*voltdimFlag)/2,1);
block10=cell(voltdimFlag,1);
block11=cell(voltdimFlag,1);
block13=cell(factorial(voltdimFlag)/(factorial(3)*factorial(voltdimFlag-3)),1);
%write text for variable names and combinations for terms 1:11, 13
for i=1:voltdimFlag
    block2(i)=leftlist(i);
    block3(i)=strcat('|',leftlist(i),'|');
    block4(i)=strcat(leftlist(i),'*',leftlist(i));
    block5(i)=strcat(leftlist(i),'*|',leftlist(i),'|');
    
    for j=i+1:voltdimFlag
        block6(count5)=strcat(leftlist(i),'*',leftlist(j));
        block7(count5)=strcat('|',leftlist(i),'*',leftlist(j),'|');
        block8(count5)=strcat(leftlist(i),'*|',leftlist(j),'|');
        block9(count5)=strcat('|',leftlist(i),'|*',leftlist(j));
        count5=count5+1;
        for k=j+1:voltdimFlag
            block13(count13)=strcat(leftlist(i),'*',leftlist(j),'*',leftlist(k));
            count13=count13+1;
        end
    end
    block10(i)=strcat(leftlist(i),'*',leftlist(i),'*',leftlist(i));
    block11(i)=strcat('|',leftlist(i),'*',leftlist(i),'*',leftlist(i),'|');
end

%write text for variable names and combinations for term 12 (F*G*G)
block12=cell(factorial(voltdimFlag)/factorial(voltdimFlag-2),1);
count=1;
for i=1:voltdimFlag
    j_ind=setdiff([1:voltdimFlag],i); %Indices for inner loop
    for j=1:length(j_ind)
        block12(count)=strcat(leftlist(i),'*',leftlist(j_ind(j)),'*',leftlist(j_ind(j)));
        count=count+1;
    end
end

%Select Terms based on model type selected
if model==3
    leftColumn =[block1;block2];
elseif model==2
    leftColumn=[block1;block2;block4;block6];
else
    leftColumn=[block1;block2;block3;block4;block5;block6;block7;block8;block9;block10;block11;block12;block13];
end

if nargin>=7
    if numRBF>0
        channel=repmat([1:loaddimFlag]',numRBF,1);
        rbf_leftColumn=cellstr(strcat(reshape(toplist(channel),numel(channel),1), repmat(' RBF ',numRBF*loaddimFlag,1), num2str(repelem([1:numRBF]',loaddimFlag,1))));
        leftColumn=[leftColumn;rbf_leftColumn];
    end
end
end
