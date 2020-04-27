%Program is a GUI to split data into calibration/validation data files.  2
%modes for splitting data are available:
%
% LHS: Data split using Latin Hypercube Sampling.  First, datapoints with
% voltages in the bottom or top 5% are automatically put into calibration
% dataset. Next, remaining points are split based on set percentages using
% LHS. This may result in not every series being represented in calibration or 
% validation datasets. The first datapoints in each series are included in both 
% calibration and validation datasets
%
% Random Split: Data split randomly. This mode is primarily intended for
% testing AOX_BalCal's ability to generalize to new data. In each series,
% set percentage of points are randomly assigned to calibration, remainder
% go to validation.  The first datapoints in each series are included in both calibration
% and validation datasets
%
% Calibration and Validation Files are output with pointIDs, series labels,
% and the top header section of the original data file included. A
% description of the data splitting process is also included in the top row


function varargout = Cal_Val_dataSplitter_GUI(varargin)
% Cal_Val_dataSplitter_GUI MATLAB code for Cal_Val_dataSplitter_GUI.fig
%      Cal_Val_dataSplitter_GUI, by itself, creates a new Cal_Val_dataSplitter_GUI or raises the existing
%      singleton*.
%
%      H = Cal_Val_dataSplitter_GUI returns the handle to a new Cal_Val_dataSplitter_GUI or the handle to
%      the existing singleton*.
%
%      Cal_Val_dataSplitter_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Cal_Val_dataSplitter_GUI.M with the given input arguments.
%
%      Cal_Val_dataSplitter_GUI('Property','Value',...) creates a new Cal_Val_dataSplitter_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Cal_Val_dataSplitter_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Cal_Val_dataSplitter_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Cal_Val_dataSplitter_GUI

% Last Modified by GUIDE v2.5 18-Jan-2020 17:40:37
%Add path for subfolder of functions
addpath(genpath('AOX_RequiredSupportFiles'))

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Cal_Val_dataSplitter_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Cal_Val_dataSplitter_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Cal_Val_dataSplitter_GUI is made visible.
function Cal_Val_dataSplitter_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Cal_Val_dataSplitter_GUI (see VARARGIN)

% Choose default command line output for Cal_Val_dataSplitter_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

stateFile = 'state.split';
if exist(stateFile)
    try
    load(stateFile,'-mat');
    set(handles.filePath,'String',state.filePath);
    set(handles.r11,'String',state.r11);
    set(handles.r12,'String',state.r12);
    set(handles.r21,'String',state.r21);
    set(handles.r22,'String',state.r22);
    set(handles.r31,'String',state.r31);
    set(handles.r32,'String',state.r32);
    set(handles.saveCalPath,'String',state.saveCalPath);
    set(handles.saveValPath,'String',state.saveValPath);
    set(handles.pEdit,'String',state.pEdit);
    set(handles.split_type,'Value',state.split_type);
    end
end

% UIWAIT makes Cal_Val_dataSplitter_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Cal_Val_dataSplitter_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function filePath_Callback(hObject, eventdata, handles)
% hObject    handle to filePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filePath as text
%        str2double(get(hObject,'String')) returns contents of filePath as a double


% --- Executes during object creation, after setting all properties.
function filePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fileFind.
function fileFind_Callback(hObject, eventdata, handles)
% hObject    handle to fileFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName_ext, PathName] = uigetfile('*.csv');
if FileName_ext ~= 0
    [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
    CurrentPath = [CurrentPath,filesep];
    if strcmp(CurrentPath,PathName)
        FullPath = FileName_ext;
    else
        FullPath = [PathName,FileName_ext];
    end
    set(handles.filePath,'String',FullPath)
end
[~,FileName,~] = fileparts(FileName_ext);
set(handles.saveCalPath,'String',[PathName,FileName,'_cal.csv']);
set(handles.saveValPath,'String',[PathName,FileName,'_val.csv']);



function r11_Callback(hObject, eventdata, handles)
% hObject    handle to r11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r11 as text
%        str2double(get(hObject,'String')) returns contents of r11 as a double


% --- Executes during object creation, after setting all properties.
function r11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r12_Callback(hObject, eventdata, handles)
% hObject    handle to r12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r12 as text
%        str2double(get(hObject,'String')) returns contents of r12 as a double


% --- Executes during object creation, after setting all properties.
function r12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r21_Callback(hObject, eventdata, handles)
% hObject    handle to r21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r21 as text
%        str2double(get(hObject,'String')) returns contents of r21 as a double


% --- Executes during object creation, after setting all properties.
function r21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r22_Callback(hObject, eventdata, handles)
% hObject    handle to r22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r22 as text
%        str2double(get(hObject,'String')) returns contents of r22 as a double


% --- Executes during object creation, after setting all properties.
function r22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r31_Callback(hObject, eventdata, handles)
% hObject    handle to r31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r31 as text
%        str2double(get(hObject,'String')) returns contents of r31 as a double


% --- Executes during object creation, after setting all properties.
function r31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r32_Callback(hObject, eventdata, handles)
% hObject    handle to r32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r32 as text
%        str2double(get(hObject,'String')) returns contents of r32 as a double


% --- Executes during object creation, after setting all properties.
function r32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pEdit as text
%        str2double(get(hObject,'String')) returns contents of pEdit as a double


% --- Executes during object creation, after setting all properties.
function pEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function saveCalPath_Callback(hObject, eventdata, handles)
% hObject    handle to saveCalPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveCalPath as text
%        str2double(get(hObject,'String')) returns contents of saveCalPath as a double


% --- Executes during object creation, after setting all properties.
function saveCalPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveCalPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveCalFind.
function saveCalFind_Callback(hObject, eventdata, handles)
% hObject    handle to saveCalFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uiputfile('*.csv');
if FileName ~= 0
    FullPath = [PathName,FileName];
    set(handles.saveCalPath,'String',FullPath)
end


function saveValPath_Callback(hObject, eventdata, handles)
% hObject    handle to saveValPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveValPath as text
%        str2double(get(hObject,'String')) returns contents of saveValPath as a double


% --- Executes during object creation, after setting all properties.
function saveValPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveValPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveValFind.
function saveValFind_Callback(hObject, eventdata, handles)
% hObject    handle to saveValFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uiputfile('*.csv');
if FileName ~= 0
    FullPath = [PathName,FileName];
    set(handles.saveValPath,'String',FullPath)
end


% --- Executes on button press in runButton.
function runButton_Callback(hObject, eventdata, handles)
% hObject    handle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure1,'Visible','off');
disp('Splitting Data...')
pause(0.1);

uiresume(handles.figure1);

filePath = get(handles.filePath,'String');
Range{1} = [get(handles.r11,'String'),'..',get(handles.r12,'String')];
Range{2} = [get(handles.r21,'String'),'..',get(handles.r22,'String')];
Range{3} = [get(handles.r31,'String'),'..',get(handles.r32,'String')];
CSV(1,:) = a12rc(get(handles.r11,'String'));
CSV(2,:) = a12rc(get(handles.r21,'String'));
CSV(3,:) = a12rc(get(handles.r31,'String'));
  
%Read series labels using 'readtable': JRP 19 June 19
        A=extractAfter(Range{1},'..');
        bottom=str2double(regexp(A,'\d*','Match'));
        opts_data=delimitedTextImportOptions('DataLines',[CSV(1,1)+1 bottom]);
        data_bulk=readtable(filePath,opts_data);
        opts_header=delimitedTextImportOptions('DataLines',[1 CSV(1,1)]);
        header_bulk=readtable(filePath,opts_header);
        s=str2double(table2array(data_bulk(:,CSV(1,2)+1)));
        series2=table2array(data_bulk(:,CSV(1,2)+2));
        pointID=table2array(data_bulk(:,CSV(1,2)));
        leftCol=table2array(data_bulk(:,CSV(1,2)-1));
        
        
        
% s = csvread(filePath,CSV(1,1),CSV(1,2),Range{1});
y = csvread(filePath,CSV(2,1),CSV(2,2),Range{2});
x = csvread(filePath,CSV(3,1),CSV(3,2),Range{3});

split_selection=get(handles.split_type, 'Value');
if split_selection==1
    [X_cal,X_val]=LHS(s,y,x,handles);
else 
    [X_cal,X_val]=rand_split(s,y,x,handles);
end
writeout_CSV(X_cal,X_val,header_bulk,leftCol,pointID,s,series2,y,x,handles);
saveState(handles);
disp('Complete');
close(get(hObject,'Parent'))

function a = a12rc(a1)
%converts spreadsheet notation "A1" to row and column numbers  (0-based)
alpha_ind = find(isletter(a1));
alpha = abs(upper(a1(alpha_ind)))-65;
if length(alpha) == 2
    c = 26*(alpha(1)+1) + alpha(2);
else
    c = alpha;
end
a1(alpha_ind) = [];
r = str2num(a1)-1;
a = [r c];

function [X_cal, X_val]=LHS(s,y,x,handles)
%x= voltages
%y= loads
p = str2num(get(handles.pEdit,'String'))/100; %percent that is sampled
thr = 0.05;

[n, dim] = size(x); %n= number of datapoint, dim=# of channels
XI = [1:n]'; %count of 1 to n

[srs, i_s1st, i_s] = unique(s); %srs= series numbers, i_s1st = starting count of each series, i_s= column of series
n_s = length(srs); %total number of series

tM_tmp = x - ones(n,1)*min(x); %tM_tmp=subtract the smallest value in each column from each value: min value is 0 in each column
x_norm = tM_tmp./( ones(n,1)*max(tM_tmp) + eps ); %normalizes so all values are between 0-1

LZ = XI(i_s1st); %same as XI ??

XI_tmp  = XI;        i_s_tmp = i_s; %XI_tmp and i_s_tmp start out as parent variables
XI_tmp(i_s1st) = []; i_s_tmp(i_s1st) = []; %removes counts that each new series starts at
x_norm_tmp = x_norm(XI_tmp,:); %removes first line of each series from normalized voltages

boundary = (x_norm_tmp<thr) + (x_norm_tmp>(1-thr)); %marks normalized voltages that are < or .05 or > .95
bnd = sum(boundary,2) > 0; %marks if each row has a voltage <.05 or >.95
BD = XI_tmp(bnd); %list of row numbers that were 1 in bnd

XI_tmp(bnd) = []; i_s_tmp(bnd) = []; %removes rows stored in BD (bottom or top .05)
x_norm_tmp = x_norm(XI_tmp,:); %removes those rows from normalized voltages

x_tmp = x_norm_tmp - ones(size(x_norm_tmp,1),1)*min(x_norm_tmp); %Resubtract so min in each column is 0
x_renorm = x_tmp./( ones(size(x_tmp,1),1)*max(x_tmp) + eps ); %Renormalize so max in each column is 0

%At this point x_renorm contains voltages without first voltage from each
%series, or test points where any channel voltage was in the bottom or top
% 5% of the range of voltage
n_sample = round(p*size(x_renorm,1)); %round to how many data points for calibration based on input %
x_lhs = lhsdesign(n_sample,dim); %returns (# dim) columns with (# samples) values between 0-1 randomly selected from (# samples) evenly spaced intervals

D = pdist2(x_renorm,x_lhs); %caculates distances between all rows in normalized voltages and lhs numbers 
%D(i,j) is distance between row i in normalized voltages, row j in lhs
%values
XI_sample = [];
for i = 1:n_sample
    [D_min1,DI] = min(D,[],2); %D_min1= min value in each row, DI = index of where this value if sound
    [~,DI_2] = min(D_min1); %DI_2= index (row number) of absolute minimum distance from every row
    DI_1 = DI(DI_2); %Find column number of min value.  This is row in x_renorm (normalized voltages)
    XI_sample = [XI_sample; XI_tmp(DI_2)]; %What original voltage row number this min value came from
    
    D(DI_2,:) = []; D(:,DI_1) = []; %remove the row and column the min came from (all the distances involving that lhs row and normalized voltage row)  
    XI_tmp(DI_2) = []; %Remove that row from the selectable voltages
end %repeat this until n_samples are collected

X_cal = sort([LZ;BD;XI_sample]); %selected percentage of rows go to calibration- plus first row in series, values from top and bottom 5% of voltage range
X_val = sort([LZ;XI_tmp]); %remainder go to validation- plus first row in series

function [X_cal, X_val]=rand_split(s,y,x,handles)
%x= voltages
%y= loads
p = str2num(get(handles.pEdit,'String'))/100; %percent that is sampled

[n, dim] = size(x); %n= number of datapoint, dim=# of channels
XI = [1:n]'; %count of 1 to n

[srs, i_s1st, i_s] = unique(s); %srs= series numbers, i_s1st = starting count of each series, i_s= column of series
n_s = length(srs); %total number of series

series_bookends=[i_s1st;n+1]; %starting point of each series, including next series after final series
n_series=series_bookends(2:end)-series_bookends(1:end-1); %number of datapoints in each series

n_cal_series=round(p*n_series); %Number of calibration datapoints in each series
n_val_series=n_series-n_cal_series; %Number Validation datapoints in each series

LZ = XI(i_s1st); %same as XI ??

rand_I_cal=zeros(sum(n_cal_series),1);
rand_I_val=zeros(sum(n_val_series),1);
for i=1:n_s
    series_min=series_bookends(i); %Starting index for series
    series_max=series_bookends(i+1)-1; %Last index for series
    
    I_cal=(series_min-1)+sort(randperm(n_series(i),n_cal_series(i)))'; %Random integers for calibration points
    I_val=setdiff([series_min:series_max]',I_cal); %Remainder of points go to validation
    
    rand_I_cal(sum(n_cal_series(1:i-1))+1:sum(n_cal_series(1:i)))=I_cal;
    rand_I_val(sum(n_val_series(1:i-1))+1:sum(n_val_series(1:i)))=I_val;
    
end
X_cal=sort(unique([LZ;rand_I_cal]));
X_val=sort(unique([LZ;rand_I_val]));

function writeout_CSV(X_cal,X_val,header_bulk,leftCol,pointID,s,series2,y,x,handles)
%Writing output
header_bulk_cell=table2cell(header_bulk);

i_comment=find(cellfun(@isempty,header_bulk{1,:}),1); %First first empty cell in first row to insert comment
cTime=datestr(now,'yyyy-mmdd-HHMMSS'); 
filePath = get(handles.filePath,'String');
split_string=handles.split_type.String(handles.split_type.Value);
comment_str=strcat('Data Split:'," ", cTime, '. Split using:'," ", split_string, '. Calibration data percentage= ', handles.pEdit.String, '%. Original File:', " ", filePath); %Comment describing split
header_bulk_cell{1,i_comment}=comment_str;

cal_out{1}=[header_bulk_cell; leftCol(X_cal), pointID(X_cal), num2cell(s(X_cal)), series2(X_cal), num2cell(y(X_cal,:)), num2cell(x(X_cal,:))];
val_out{1}=[header_bulk_cell; leftCol(X_val), pointID(X_val), num2cell(s(X_val)), series2(X_val), num2cell(y(X_val,:)), num2cell(x(X_val,:))];
disp('Writing Calibration File...');
writetable(cell2table(cal_out{1}),get(handles.saveCalPath,'String'),'writevariablenames',0); %write to csv
disp('Writing Validation File... ');
writetable(cell2table(val_out{1}),get(handles.saveValPath,'String'),'writevariablenames',0); %write to csv
% csvwrite(get(handles.saveCalPath,'String'),[s(X_cal),y(X_cal,:),x(X_cal,:)]);
% csvwrite(get(handles.saveValPath,'String'),[s(X_val),y(X_val,:),x(X_val,:)]);

% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(get(hObject,'Parent'))

function saveState(handles)
state.filePath = get(handles.filePath,'String');
state.r11 = get(handles.r11,'String');
state.r12 = get(handles.r12,'String');
state.r21 = get(handles.r21,'String');
state.r22 = get(handles.r22,'String');
state.r31 = get(handles.r31,'String');
state.r32 = get(handles.r32,'String');
state.saveCalPath = get(handles.saveCalPath,'String');
state.saveValPath = get(handles.saveValPath,'String');
state.pEdit = get(handles.pEdit,'String');
state.split_type=get(handles.split_type,'Value');
save('state.split','state');

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on selection change in split_type.
function split_type_Callback(hObject, eventdata, handles)
% hObject    handle to split_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns split_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from split_type


% --- Executes during object creation, after setting all properties.
function split_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to split_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
