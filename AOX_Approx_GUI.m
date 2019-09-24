% Version 11: Last modified on 6/11/18
function varargout = AOX_Approx_GUI(varargin)
% AOX_Approx_GUI MATLAB code for AOX_Approx_GUI.fig
%      AOX_Approx_GUI, by itself, creates a new AOX_Approx_GUI or raises the existing
%      singleton*.
%
%      H = AOX_Approx_GUI returns the handle to a new AOX_Approx_GUI or the handle to
%      the existing singleton*.
%
%      AOX_Approx_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AOX_Approx_GUI.M with the given input arguments.
%
%      AOX_Approx_GUI('Property','Value',...) creates a new AOX_Approx_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AOX_Approx_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AOX_Approx_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AOX_Approx_GUI

% Last Modified by GUIDE v2.5 01-Aug-2019 15:53:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AOX_Approx_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @AOX_Approx_GUI_OutputFcn, ...
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


% --- Executes just before AOX_Approx_GUI is made visible.
function AOX_Approx_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AOX_Approx_GUI (see VARARGIN)
set(handles.figure1, 'units', 'normalized', 'position', [0.15 0.1 0.7 0.7])
set(handles.output_location,'String',pwd);
global VERSION
VERSION = 16;
try
    
    [nasalogo,~,aln] = imread('nasa.png','BackgroundColor',[0.941, 0.941, 0.941]);
    axes(handles.axesNASA);
    imshow(nasalogo, []);
    
    [ricelogo,~,alr] = imread('rice.png','BackgroundColor',[0.941, 0.941, 0.941]);
    axes(handles.axesRice);
    imshow(ricelogo, []);
    
end

% Choose default command line output for AOX_Approx_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes AOX_Approx_GUI wait for user response (see UIRESUME)
[CurrentPath,~,~] = fileparts(mfilename('fullpath'));
fileName = [CurrentPath,filesep,'default_approx.ini'];

if exist(fileName,'file')
    try
        load(fileName,'-mat');
        versionCheck(default);
        %Calibration Model Group
        set(handles.calib_model_path,'String',default.calib_model_path);
        set(handles.grbf,'Value',default.grbf);
        %Uncertainty Group
        set(handles.loadPI_FLAGcheck,'Value',default.loadPI);
        set(handles.PI_percent_confidence,'String',default.PI_percent_confidence)
        %Approximation Group
        set(handles.appPath,'String',default.appPath);
        set(handles.a11,'String',default.appRange{1,1});
        set(handles.a12,'String',default.appRange{1,2});
        set(handles.a21,'String',default.appRange{2,1});
        set(handles.a22,'String',default.appRange{2,2});
        set(handles.a31,'String',default.appRange{3,1});
        set(handles.a32,'String',default.appRange{3,2});
        set(handles.a41,'String',default.appRange{4,1});
        set(handles.a42,'String',default.appRange{4,2});
        %Outputs Group
        set(handles.output_location,'String',default.output_location);
        set(handles.output_to_approx_FLAG,'Value',default.output_to_approx_FLAG);
        set(handles.subfolder_FLAG,'Value',default.subfolder_FLAG);
        set(handles.excel_FLAGcheck,'Value',default.excel);
        set(handles.approx_and_PI_print,'Value',default.approx_and_PI_print);
        set(handles.PI_print,'Value',default.PI_print);
        %Callbacks
        calib_model_path_Callback(handles.calib_model_path, eventdata, handles)
        appPath_Callback(handles.appPath, eventdata, handles)
        output_to_approx_FLAG_Callback(handles.output_to_approx_FLAG, eventdata, handles)
        loadPI_FLAGcheck_Callback(handles.loadPI_FLAGcheck, eventdata, handles)
        
    catch
        disp('local default_approx.ini may be outdated or incompatible with GUI.');
    end
end

uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AOX_Approx_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.figure1);


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);

%Calibration Model Group
outStruct.calib_model_path=get(handles.calib_model_path,'String');
outStruct.grbf = 1 + get(handles.grbf,'Value');
%Uncertainty Group
outStruct.loadPI = get(handles.loadPI_FLAGcheck,'Value');
outStruct.PI_percent_confidence=str2num(get(handles.PI_percent_confidence,'String'));
%Approximation Group
app.type = 'approximate';
app.Path = get(handles.appPath,'String');
[~,~,appext] = fileparts(app.Path);
switch appext
    case '.csv'
        app.Range{1} = [get(handles.a11,'String'),'..',get(handles.a12,'String')];
        app.CSV(1,:) = a12rc(get(handles.a11,'String'));
        app.Range{2} = [get(handles.a21,'String'),'..',get(handles.a22,'String')];
        app.CSV(2,:) = a12rc(get(handles.a21,'String'));
        app.Range{3} = [get(handles.a31,'String'),'..',get(handles.a32,'String')];
        app.CSV(3,:) = a12rc(get(handles.a31,'String'));
        app.Range{4} = [get(handles.a41,'String'),'..',get(handles.a42,'String')];
        app.CSV(4,:) = a12rc(get(handles.a41,'String'));
        outStruct.savePathapp = loadCSV(app);
    case '.app'
        outStruct.savePathapp = app.Path;
end
%Outputs Group
outStruct.excel = get(handles.excel_FLAGcheck,'Value');
outStruct.approx_and_PI_print = get(handles.approx_and_PI_print,'Value');
outStruct.PI_print = get(handles.PI_print,'Value');
outStruct.output_location=get(handles.output_location,'String');
outStruct.subfolder_FLAG=get(handles.subfolder_FLAG,'Value');

%Make new subfolder if selected as option
%Default output location to current directory if empty
if isempty(outStruct.output_location)==1
    outStruct.output_location=cd;
end

outStruct.REPORT_NO=datestr(now,'yyyy-mmdd-HHMMSS');
outStruct.output_location=[outStruct.output_location,filesep];
if outStruct.subfolder_FLAG==1
    try
        new_subpath=fullfile(outStruct.output_location,['AOX_BalCal_Results_',outStruct.REPORT_NO]);
        mkdir(char(new_subpath));
        outStruct.output_location=[new_subpath,filesep];
    catch
        fprintf('Unable to create new subfolder. Saving results in: ');
        fprintf('%s',outStruct.output_location); fprintf('\n');
    end
end

outStruct.cancel = 0;

handles.output = outStruct;
guidata(hObject, handles);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    outStruct.cancel = 1;
    handles.output = outStruct;
    guidata(hObject,handles);
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
outStruct.cancel = 1;
handles.output = outStruct;
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in ini_button.
function ini_button_Callback(hObject, eventdata, handles)
% hObject    handle to ini_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global VERSION
default.version = VERSION;

%Calibration Model Group
        default.calib_model_path=get(handles.calib_model_path,'String');
        default.grbf=get(handles.grbf,'Value');
        %Uncertainty Group
        default.loadPI=get(handles.loadPI_FLAGcheck,'Value');
        default.PI_percent_confidence=get(handles.PI_percent_confidence,'String');
        %Approximation Group
        default.appPath=get(handles.appPath,'String');
        default.appRange{1,1}=get(handles.a11,'String');
        default.appRange{1,2}=get(handles.a12,'String');
        default.appRange{2,1}=get(handles.a21,'String');
        default.appRange{2,2}=get(handles.a22,'String');
        default.appRange{3,1}=get(handles.a31,'String');
        default.appRange{3,2}=get(handles.a32,'String');
        default.appRange{4,1}=get(handles.a41,'String');
        default.appRange{4,2}=get(handles.a42,'String');
        %Outputs Group
        default.output_location=get(handles.output_location,'String');
        default.output_to_approx_FLAG=get(handles.output_to_approx_FLAG,'Value');
        default.subfolder_FLAG=get(handles.subfolder_FLAG,'Value');
        default.excel=get(handles.excel_FLAGcheck,'Value');
        default.approx_and_PI_print=get(handles.approx_and_PI_print,'Value');
        default.PI_print=get(handles.PI_print,'Value');


[CurrentPath,~,~] = fileparts(mfilename('fullpath'));
fileName = [CurrentPath,filesep,'default_approx.ini'];
save(fileName,'default');
set(hObject,'String','Saving...');
pause(1);
set(hObject,'String','Saved');
pause(0.5);
set(hObject,'String','Set as Default...')

function a = a12rc(a1)
%converts spreadsheet notation "A1" to row and column numbers  (0-based)
alpha_ind = find(isletter(a1));
alpha = abs(upper(a1(alpha_ind)))-65;
c = alpha;
a1(alpha_ind) = [];
r = str2num(a1)-1;
a = [r c];

function a1 = rc2a1(a)
r = a(1); c = a(2);
alpha = char(c+65);
num = int2str(r+1);
a1 = [alpha, num];


function a11_Callback(hObject, eventdata, handles)
% hObject    handle to a11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a11 as text
%        str2double(get(hObject,'String')) returns contents of a11 as a double


% --- Executes during object creation, after setting all properties.
function a11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a12_Callback(hObject, eventdata, handles)
% hObject    handle to a12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a12 as text
%        str2double(get(hObject,'String')) returns contents of a12 as a double


% --- Executes during object creation, after setting all properties.
function a12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function appPath_Callback(hObject, eventdata, handles)
% hObject    handle to appPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of appPath as text
%        str2double(get(hObject,'String')) returns contents of appPath as a double
[~,~,ext] = fileparts(get(hObject,'String'));
try
    if strcmp(ext,'.csv')
        %set(handles.appSave, 'Enable', 'on');
        set(handles.a11, 'Enable', 'on');
        set(handles.a12, 'Enable', 'on');
        set(handles.a21, 'Enable', 'on');
        set(handles.a22, 'Enable', 'on');
        set(handles.a31, 'Enable', 'on');
        set(handles.a32, 'Enable', 'on');
        set(handles.a41, 'Enable', 'on');
        set(handles.a42, 'Enable', 'on');
    else
        %set(handles.appSave, 'Enable', 'off');
        if strcmp(ext,'.app')
            load(get(hObject,'String'), '-mat', 'app');
            splitrange = split(app.Range,'..');
        else
            splitrange=cell(1,4,2);
        end
        
        set(handles.a11, 'Enable', 'off', 'String', splitrange{1,1,1});
        set(handles.a12, 'Enable', 'off', 'String', splitrange{1,1,2});
        set(handles.a21, 'Enable', 'off', 'String', splitrange{1,2,1});
        set(handles.a22, 'Enable', 'off', 'String', splitrange{1,2,2});
        set(handles.a31, 'Enable', 'off', 'String', splitrange{1,3,1});
        set(handles.a32, 'Enable', 'off', 'String', splitrange{1,3,2});
        set(handles.a41, 'Enable', 'off', 'String', splitrange{1,4,1});
        set(handles.a42, 'Enable', 'off', 'String', splitrange{1,4,2});
    end
    output_to_approx_FLAG_Callback(handles.output_to_approx_FLAG, eventdata, handles)

catch
    disp('Problem occurred while reading Approximation file');
end

% --- Executes during object creation, after setting all properties.
function appPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to appPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in appFind.
function appFind_Callback(hObject, eventdata, handles)
% hObject    handle to appFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.csv;*.app');
if FileName ~= 0
    [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
    CurrentPath = [CurrentPath,'\'];
    if strcmp(CurrentPath,PathName)
        FullPath = FileName;
    else
        FullPath = [PathName,FileName];
    end
    set(handles.appPath,'String',FullPath)
end
appPath_Callback(handles.appPath, eventdata, handles)


% --- Executes on button press in appSave.
function appSave_Callback(hObject, eventdata, handles)
% hObject    handle to appSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function savePath = loadCSV(cva)
% This function pre-loads the csv's and saves the data as .mat for quicker
% reading.
% Input: type - String that changes depending on whether loading
% calibration, validation, or approximation data.
switch cva.type
    case 'approximate'
        app = cva;
        
        loadCapacitiesapprox =    csvread(app.Path,app.CSV(1,1),app.CSV(1,2),app.Range{1});
        natzerosapprox =          csvread(app.Path,app.CSV(2,1),app.CSV(2,2),app.Range{2});
        
       %Read series labels using 'readtable': JRP 19 June 19
        A=extractAfter(app.Range{3},'..');
        bottom=str2double(regexp(A,'\d*','Match'));
        opts=delimitedTextImportOptions('DataLines',[app.CSV(3,1)+1 bottom]);
        series_bulk=readtable(app.Path,opts);
        seriesapprox=str2double(table2array(series_bulk(:,app.CSV(3,2)+1)));
        series2approx=table2array(series_bulk(:,app.CSV(3,2)+2));
        pointIDapprox=table2array(series_bulk(:,app.CSV(3,2)));
        clear A bottom opts series_bulk 
%         seriesapprox =            csvread(app.Path,app.CSV(3,1),app.CSV(3,2),app.Range{3});
  
        excessVecapprox =         csvread(app.Path,app.CSV(4,1),app.CSV(4,2),app.Range{4});
        
        try
            file=fopen(app.Path); %open file
            all_text1 = textscan(file,'%s','Delimiter','\n'); %read in all text
            fclose(file); %close file
            %Eliminate rows with ";" in first column
            lastrow=a12rc(extractAfter(app.Range{3},".."));
            all_text_points=all_text1{1}(app.CSV(4,1)+1:lastrow(1)+1);
            for i=1:size(all_text_points,1)
                all_text_points_split(i,:)=cellstr(strsplit(string(all_text_points(i)),',','CollapseDelimiters',false)); %Extract row with labels
            end
            first_col=all_text_points_split(:,1);
            ignore_row=find(contains(first_col,';')); %Find rows with semicolons in the first column
            
            excessVecapprox(ignore_row,:)=[];
            pointIDapprox(ignore_row,:)=[];
            seriesapprox(ignore_row,:)=[];
            series2approx(ignore_row,:)=[];
        catch
            fprintf('\n UNABLE TO REMOVE ROWS FLAGGED WITH ";" FROM INPUT FILE \n')
        end
        
        [~,appName,~] = fileparts(app.Path);
        fileNameapprox = [appName,'.app'];
        [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
        savePathapprox = [CurrentPath,filesep,fileNameapprox];
        
        clear cva appName CurrentPath
        save(savePathapprox);
        savePath = savePathapprox;
end

function versionCheck(default)
global VERSION
mesg = 'An outdated version of default.ini was detected.';
if ~isfield(default, 'version')
    error(mesg);
elseif default.version ~= VERSION
    error(mesg);
end


function a21_Callback(hObject, eventdata, handles)
% hObject    handle to a21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a21 as text
%        str2double(get(hObject,'String')) returns contents of a21 as a double


% --- Executes during object creation, after setting all properties.
function a21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a22_Callback(hObject, eventdata, handles)
% hObject    handle to a22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a22 as text
%        str2double(get(hObject,'String')) returns contents of a22 as a double


% --- Executes during object creation, after setting all properties.
function a22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function a31_Callback(hObject, eventdata, handles)
% hObject    handle to a31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a31 as text
%        str2double(get(hObject,'String')) returns contents of a31 as a double


% --- Executes during object creation, after setting all properties.
function a31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a32_Callback(hObject, eventdata, handles)
% hObject    handle to a32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a32 as text
%        str2double(get(hObject,'String')) returns contents of a32 as a double


% --- Executes during object creation, after setting all properties.
function a32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a41_Callback(hObject, eventdata, handles)
% hObject    handle to a41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a41 as text
%        str2double(get(hObject,'String')) returns contents of a41 as a double


% --- Executes during object creation, after setting all properties.
function a41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a42_Callback(hObject, eventdata, handles)
% hObject    handle to a42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a42 as text
%        str2double(get(hObject,'String')) returns contents of a42 as a double


% --- Executes during object creation, after setting all properties.
function a42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in excel_FLAGcheck.
function excel_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to excel_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of excel_FLAGcheck

% --- Executes on button press in loadPI_FLAGcheck.
function loadPI_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to loadPI_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadPI_FLAGcheck
if get(hObject,'Value') == 0
    set(handles.approx_and_PI_print,'Enable','off','Value',0);
    set(handles.PI_print,'Enable','off','Value',0);
    set(handles.percent_confidence_text,'Enable','off');
    set(handles.PI_percent_confidence,'Enable','off');
else
     set(handles.approx_and_PI_print,'Enable','on');
     set(handles.PI_print,'Enable','on');
     set(handles.percent_confidence_text,'Enable','on');
     set(handles.PI_percent_confidence,'Enable','on');
end


% --- Executes on button press in approx_and_PI_print.
function approx_and_PI_print_Callback(hObject, eventdata, handles)
% hObject    handle to approx_and_PI_print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of approx_and_PI_print


% --- Executes on button press in PI_print.
function PI_print_Callback(hObject, eventdata, handles)
% hObject    handle to PI_print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PI_print


function output_location_Callback(hObject, eventdata, handles)
% hObject    handle to output_location (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_location as text
%        str2double(get(hObject,'String')) returns contents of output_location as a double


% --- Executes during object creation, after setting all properties.
function output_location_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_location (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in output_location_button.
function output_location_button_Callback(hObject, eventdata, handles)
% hObject    handle to output_location_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selpath=uigetdir;
set(handles.output_location,'String',selpath);


% --- Executes on button press in output_to_approx_FLAG.
function output_to_approx_FLAG_Callback(hObject, eventdata, handles)
% hObject    handle to output_to_approx_FLAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of output_to_approx_FLAG
if get(hObject,'Value') == 0
    set(handles.output_location_button,'Enable','on');
    set(handles.output_location,'Enable','on');
    set(handles.output_location,'String',pwd);
else
    set(handles.output_location_button,'Enable','off');
    set(handles.output_location,'Enable','off');
    approx_path=get(handles.appPath,'String');
    if isempty(approx_path)==0
        [approx_path,~,~] = fileparts(approx_path);
%         approx_path=extractBefore(approx_path,find(approx_path == '\', 1, 'last'));
    else
        approx_path=pwd;
    end
    set(handles.output_location,'String',approx_path);
end


% --- Executes on button press in subfolder_FLAG.
function subfolder_FLAG_Callback(hObject, eventdata, handles)
% hObject    handle to subfolder_FLAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subfolder_FLAG



% --- Executes on button press in calibFind.
function calibFind_Callback(hObject, eventdata, handles)
% hObject    handle to calibFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.mat');
if FileName ~= 0
    [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
    CurrentPath = [CurrentPath,'\'];
    if strcmp(CurrentPath,PathName)
        FullPath = FileName;
    else
        FullPath = [PathName,FileName];
    end
    set(handles.calib_model_path,'String',FullPath)
end
calib_model_path_Callback(handles.calib_model_path, eventdata, handles)


function calib_model_path_Callback(hObject, eventdata, handles)
% hObject    handle to calib_model_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calib_model_path as text
%        str2double(get(hObject,'String')) returns contents of calib_model_path as a double
try
    load(char(get(hObject,'String')));
end
if exist('coeff','var') == 1 && exist('ANOVA','var')==1 && exist('loadlist','var')==1 && exist('model','var')==1
    load_message=strcat(num2str(size(coeff,2)), ' Channel Calibration Model Loaded');
    set(handles.master_warn_text,'String',load_message,'ForegroundColor','blue');
    set(handles.runbutton,'Enable','on');
    if isfield(ANOVA,'PI')==1
        set(handles.loadPI_FLAGcheck,'Enable','on');
        set(handles.PI_text,'String','');
    else
        set(handles.loadPI_FLAGcheck,'Enable','off','Value',0);
        set(handles.PI_text,'String','ANOVA NOT FOUND IN CALIBRATION FILE');
    end
    loadPI_FLAGcheck_Callback(handles.loadPI_FLAGcheck, eventdata, handles)
    
    if exist('cHist','var')==1 && exist('center_daHist','var')==1 && exist('wHist','var')==1
        set(handles.grbf,'Enable','on');
        grbf_message=strcat(num2str(size(wHist,1)), ' GRBFs per channel');
        set(handles.GRBF_text,'String',grbf_message, 'ForegroundColor','blue');
    else
        set(handles.grbf,'Enable','off','Value',0);
        set(handles.GRBF_text,'String','GRBFS NOT FOUND IN CALIBRATION FILE','ForegroundColor','red');
    end
else
    set(handles.master_warn_text,'String','INVALID CALIBRATION .mat FILE','ForegroundColor','red');
    set(handles.runbutton,'Enable','off');
    set(handles.grbf,'Enable','off','Value',0);
    set(handles.GRBF_text,'String','GRBFS NOT FOUND IN CALIBRATION FILE','ForegroundColor','red');
    set(handles.loadPI_FLAGcheck,'Enable','off','Value',0);
    set(handles.PI_text,'String','ANOVA NOT FOUND IN CALIBRATION FILE');
    
end

% --- Executes during object creation, after setting all properties.
function calib_model_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calib_model_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PI_percent_confidence_Callback(hObject, eventdata, handles)
% hObject    handle to PI_percent_confidence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PI_percent_confidence as text
%        str2double(get(hObject,'String')) returns contents of PI_percent_confidence as a double


% --- Executes during object creation, after setting all properties.
function PI_percent_confidence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PI_percent_confidence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in grbf.
function grbf_Callback(hObject, eventdata, handles)
% hObject    handle to grbf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of grbf
