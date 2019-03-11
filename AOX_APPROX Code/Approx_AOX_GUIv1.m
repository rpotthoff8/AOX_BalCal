% Version 10: Last modified on 5/16/18
function varargout = Approx_AOX_GUIv1(varargin)
% APPROX_AOX_GUIV1 MATLAB code for Approx_AOX_GUIv1.fig
%      APPROX_AOX_GUIV1, by itself, creates a new APPROX_AOX_GUIV1 or raises the existing
%      singleton*.
%
%      H = APPROX_AOX_GUIV1 returns the handle to a new APPROX_AOX_GUIV1 or the handle to
%      the existing singleton*.
%
%      APPROX_AOX_GUIV1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APPROX_AOX_GUIV1.M with the given input arguments.
%
%      APPROX_AOX_GUIV1('Property','Value',...) creates a new APPROX_AOX_GUIV1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Approx_AOX_GUIv1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Approx_AOX_GUIv1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Approx_AOX_GUIv1

% Last Modified by GUIDE v2.5 31-May-2018 14:17:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Approx_AOX_GUIv1_OpeningFcn, ...
                   'gui_OutputFcn',  @Approx_AOX_GUIv1_OutputFcn, ...
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


% --- Executes just before Approx_AOX_GUIv1 is made visible.
function Approx_AOX_GUIv1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Approx_AOX_GUIv1 (see VARARGIN)

global VERSION
VERSION = 1;

[nasalogo,~,aln] = imread('nasa.png','BackgroundColor',[0.941, 0.941, 0.941]);
axes(handles.axesNASA);
imshow(nasalogo, []);

[ricelogo,~,alr] = imread('rice.png','BackgroundColor',[0.941, 0.941, 0.941]);
axes(handles.axesRice);
imshow(ricelogo, []);

% Choose default command line output for Approx_AOX_GUIv1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes Approx_AOX_GUIv1 wait for user response (see UIRESUME)
[CurrentPath,~,~] = fileparts(mfilename('fullpath'));
fileName = [CurrentPath,filesep,'approx.ini'];
if exist(fileName,'file')
    try
        
        load(fileName,'-mat');

        versionCheck(default);

        %set(handles.tares_FLAGcheck,'Value',default.tares);
        %set(handles.coeff_FLAGcheck,'Value',default.coeff);
        %coeff_FLAGcheck_Callback(handles.coeff_FLAGcheck, eventdata, handles);
        %set(handles.grbftares_FLAGcheck,'Value',default.grbftares);
        %set(handles.tables_FLAGcheck,'Value',default.tables);
        %tables_FLAGcheck_Callback(handles.tables_FLAGcheck, eventdata, handles);
        %set(handles.res_FLAGcheck,'Value',default.res);
        %set(handles.hist_FLAGcheck,'Value',default.hist);

        set(handles.calPath,'String',default.calPath);
        set(handles.c11,'String',default.calRange{1,1});
        set(handles.c12,'String',default.calRange{1,2});
        set(handles.c21,'String',default.calRange{2,1});
        set(handles.c22,'String',default.calRange{2,2});
        calPath_Callback(handles.calPath, eventdata, handles)

        set(handles.appPath,'String',default.appPath);
        set(handles.a11,'String',default.appRange{1,1});
        set(handles.a12,'String',default.appRange{1,2});
        set(handles.a21,'String',default.appRange{2,1});
        set(handles.a22,'String',default.appRange{2,2});
        set(handles.a31,'String',default.appRange{3,1});
        set(handles.a32,'String',default.appRange{3,2});
        appPath_Callback(handles.appPath, eventdata, handles)
        
        set(handles.coeffPath,'String',default.coeffPath);
        
        %set(handles.full,'Value',default.full);
        %set(handles.truncated,'Value',default.truncated);
        %set(handles.linear,'Value',default.linear);

        %set(handles.direct,'Value',default.direct);
        %set(handles.indirect,'Value',default.indirect);

        set(handles.grbf,'Value',default.grbf);
        %set(handles.grbfcoeff_FLAGcheck,'Value',default.grbf_coeff);
        %set(handles.loglog_FLAGcheck,'Value',default.loglog);
        grbf_Callback(handles.grbf, eventdata, handles);
        set(handles.grbfcoeffPath,'String',default.grbfcoeffPath);
        set(handles.grbfwPath,'String',default.grbfwPath);
        set(handles.grbfcPath,'String',default.grbfcPath);
    end
end

uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Approx_AOX_GUIv1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.figure1);


% --- Executes on button press in tares_FLAGcheck.
function tares_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tares_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tares_FLAGcheck


% --- Executes on button press in hist_FLAGcheck.
function hist_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to hist_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hist_FLAGcheck


% --- Executes on button press in loglog_FLAGcheck.
function loglog_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to loglog_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loglog_FLAGcheck


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
%outStruct.approach = get(handles.indirect,'Value');
%outStruct.tares = get(handles.tares_FLAGcheck,'Value');
%outStruct.coeff = get(handles.coeff_FLAGcheck,'Value');
%outStruct.grbftares = get(handles.grbftares_FLAGcheck,'Value');
%outStruct.tables = 2*get(handles.tables_FLAGcheck,'Value') + get(handles.coeff_FLAGcheck,'Value');
%outStruct.res = get(handles.res_FLAGcheck,'Value');
%outStruct.hist = get(handles.hist_FLAGcheck,'Value');
%outStruct.loglog = get(handles.loglog_FLAGcheck,'Value');

%switch get(get(handles.modelPanel,'SelectedObject'),'Tag');
%    case 'full', outStruct.model = 1;
%    case 'truncated', outStruct.model = 2;
%    case 'linear', outStruct.model = 3;
%end
    
outStruct.coeffPath = get(handles.coeffPath,'String');

outStruct.grbf = 1 + get(handles.grbf,'Value');
outStruct.grbfcoeffPath = get(handles.grbfcoeffPath,'String');
outStruct.grbfwPath = get(handles.grbfwPath,'String');
outStruct.grbfcPath = get(handles.grbfcPath,'String');

cal.type = 'calibrate';
cal.Path = get(handles.calPath,'String');
[~,~,calext] = fileparts(cal.Path);
switch calext
    case '.csv'
        cal.Range{1} = [get(handles.c11,'String'),'..',get(handles.c12,'String')];
        cal.CSV(1,:) = a12rc(get(handles.c11,'String'));
        cal.Range{2} = [get(handles.c21,'String'),'..',get(handles.c22,'String')];
        cal.CSV(2,:) = a12rc(get(handles.c21,'String'));
        outStruct.savePathcal = loadCSV(cal);
    case '.appcal'
        outStruct.savePathcal = cal.Path;
end

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
        outStruct.savePathapp = loadCSV(app);
    case '.approx'
        outStruct.savePathapp = app.Path;
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


% --- Executes when selected object is changed in modelPanel.
function modelPanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in modelPanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if (hObject == handles.full)
    %set(handles.numBasisIn, 'String', 'First');
else
    %set(handles.numBasisIn, 'String', 'Second');
end

function calcsv_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in modelPanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if (hObject == handles.direct)
    %set(handles.numBasisIn, 'String', 'First');
else
    %set(handles.numBasisIn, 'String', 'Second');
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

function calPath_Callback(hObject, eventdata, handles)
% hObject    handle to calPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calPath as text
%        str2double(get(hObject,'String')) returns contents of calPath as a double
[~,~,ext] = fileparts(get(hObject,'String'));
if strcmp(ext,'.csv')
    %set(handles.calSave, 'Enable', 'on');
    set(handles.c11, 'Enable', 'on');
    set(handles.c12, 'Enable', 'on');
    set(handles.c21, 'Enable', 'on');
    set(handles.c22, 'Enable', 'on');
else
    %set(handles.calSave, 'Enable', 'off');
    set(handles.c11, 'Enable', 'off');
    set(handles.c12, 'Enable', 'off');
    set(handles.c21, 'Enable', 'off');
    set(handles.c22, 'Enable', 'off');
end


% --- Executes during object creation, after setting all properties.
function calPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calFind.
function calFind_Callback(hObject, eventdata, handles)
% hObject    handle to calFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.csv;*.appcal');
if FileName ~= 0
    [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
    CurrentPath = [CurrentPath,filesep];
    if strcmp(CurrentPath,PathName)
        FullPath = FileName;
    else
        FullPath = [PathName,FileName];
    end
    set(handles.calPath,'String',FullPath)
end
calPath_Callback(handles.calPath, eventdata, handles);

% --- Executes on button press in grbf.
function grbf_Callback(hObject, eventdata, handles)
% hObject    handle to grbf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of grbf
if get(hObject,'Value') == 1
    set(handles.grbfcoeffPath,'Enable','on');
    set(handles.grbfcoeffFind,'Enable','on');
    set(handles.grbfwPath,'Enable','on');
    set(handles.grbfwFind,'Enable','on');
    set(handles.grbfcPath,'Enable','on');
    set(handles.grbfcFind,'Enable','on');
    
    set(handles.calPath,'Enable','on');
    set(handles.calFind,'Enable','on');
    set(handles.c11,'Enable','on');
    set(handles.c12,'Enable','on');
    set(handles.c21,'Enable','on');
    set(handles.c22,'Enable','on');
else
    set(handles.grbfcoeffPath,'Enable','off');
    set(handles.grbfcoeffFind,'Enable','off');
    set(handles.grbfwPath,'Enable','off');
    set(handles.grbfwFind,'Enable','off');
    set(handles.grbfcPath,'Enable','off');
    set(handles.grbfcFind,'Enable','off');
    
    set(handles.calPath,'Enable','off');
    set(handles.calFind,'Enable','off');
    set(handles.c11,'Enable','off');
    set(handles.c12,'Enable','off');
    set(handles.c21,'Enable','off');
    set(handles.c22,'Enable','off');
end

function c11_Callback(hObject, eventdata, handles)
% hObject    handle to c11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c11 as text
%        str2double(get(hObject,'String')) returns contents of c11 as a double


% --- Executes during object creation, after setting all properties.
function c11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function c21_Callback(hObject, eventdata, handles)
% hObject    handle to c21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c21 as text
%        str2double(get(hObject,'String')) returns contents of c21 as a double

% --- Executes during object creation, after setting all properties.
function c21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function c12_Callback(hObject, eventdata, handles)
% hObject    handle to c12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c12 as text
%        str2double(get(hObject,'String')) returns contents of c12 as a double


% --- Executes during object creation, after setting all properties.
function c12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c22_Callback(hObject, eventdata, handles)
% hObject    handle to c22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c22 as text
%        str2double(get(hObject,'String')) returns contents of c22 as a double


% --- Executes during object creation, after setting all properties.
function c22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in grbfcoeff_FLAGcheck.
function grbfcoeff_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to grbfcoeff_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of grbfcoeff_FLAGcheck


% --- Executes on button press in coeff_FLAGcheck.
function coeff_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to coeff_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coeff_FLAGcheck
if get(hObject,'Value') == 1
    set(handles.tables_FLAGcheck,'Enable','off');
else
    set(handles.tables_FLAGcheck,'Enable','on');
end

% --- Executes on button press in tables_FLAGcheck.
function tables_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tables_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tables_FLAGcheck
if get(hObject,'Value') == 1
    set(handles.coeff_FLAGcheck,'Enable','off');
else
    set(handles.coeff_FLAGcheck,'Enable','on');
end


% --- Executes on button press in ini_button.
function ini_button_Callback(hObject, eventdata, handles)
% hObject    handle to ini_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global VERSION
default.version = VERSION;

%default.tares = get(handles.tares_FLAGcheck,'Value');
%default.approach = get(handles.indirect,'Value');
%default.coeff = get(handles.coeff_FLAGcheck,'Value');
%default.grbftares = get(handles.grbftares_FLAGcheck,'Value');
%default.tables = get(handles.tables_FLAGcheck,'Value');
%default.res = get(handles.res_FLAGcheck,'Value');
%default.hist = get(handles.hist_FLAGcheck,'Value');

default.calPath = get(handles.calPath,'String');
default.calRange{1,1} = get(handles.c11,'String');
default.calRange{1,2} = get(handles.c12,'String');
default.calRange{2,1} = get(handles.c21,'String');
default.calRange{2,2} = get(handles.c22,'String');

default.appPath = get(handles.appPath,'String');
default.appRange{1,1} = get(handles.a11,'String');
default.appRange{1,2} = get(handles.a12,'String');
default.appRange{2,1} = get(handles.a21,'String');
default.appRange{2,2} = get(handles.a22,'String');
default.appRange{3,1} = get(handles.a31,'String');
default.appRange{3,2} = get(handles.a32,'String');

%default.full = get(handles.full,'Value');
%default.truncated = get(handles.truncated,'Value');
%default.linear = get(handles.linear,'Value');

default.coeffPath = get(handles.coeffPath,'String');

default.grbf = get(handles.grbf,'Value');
default.grbfcoeffPath = get(handles.grbfcoeffPath,'String');
default.grbfwPath = get(handles.grbfwPath,'String');
default.grbfcPath = get(handles.grbfcPath,'String');

%default.grbf_coeff = get(handles.grbfcoeff_FLAGcheck,'Value');
%default.loglog = get(handles.loglog_FLAGcheck,'Value');

%default.direct = get(handles.direct,'Value');
%default.indirect = get(handles.indirect,'Value');

[CurrentPath,~,~] = fileparts(mfilename('fullpath'));
fileName = [CurrentPath,filesep,'approx.ini'];
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
if length(alpha) == 2
    c = 26*(alpha(1)+1) + alpha(2);
else
    c = alpha;
end
a1(alpha_ind) = [];
r = str2num(a1)-1;
a = [r c];


% --- Executes on button press in grbftares_FLAGcheck.
function grbftares_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to grbftares_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of grbftares_FLAGcheck


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
if strcmp(ext,'.csv')
    %set(handles.appSave, 'Enable', 'on');
    set(handles.a11, 'Enable', 'on');
    set(handles.a12, 'Enable', 'on');
    set(handles.a21, 'Enable', 'on');
    set(handles.a22, 'Enable', 'on');
    set(handles.a31, 'Enable', 'on');
    set(handles.a32, 'Enable', 'on');
else
    %set(handles.appSave, 'Enable', 'off');
    set(handles.a11, 'Enable', 'off');
    set(handles.a12, 'Enable', 'off');
    set(handles.a21, 'Enable', 'off');
    set(handles.a22, 'Enable', 'off');
    set(handles.a31, 'Enable', 'off');
    set(handles.a32, 'Enable', 'off');
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
[FileName, PathName] = uigetfile('*.csv;*.approx');
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

function savePath = loadCSV(cva)
% This function pre-loads the csv's and saves the data as .mat for quicker
% reading.
% Input: type - String that changes depending on whether loading
% calibration, validation, or approximation data.
switch cva.type
    case 'calibrate'
        cal = cva;
        
        natzerosCalib  =    csvread(cal.Path,cal.CSV(1,1),cal.CSV(1,2),cal.Range{1});
        excessVecCalib =    csvread(cal.Path,cal.CSV(2,1),cal.CSV(2,2),cal.Range{2});
        
        [~,calName,~] = fileparts(cal.Path);
        fileName = [calName,'.appcal'];
        [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
        savePath = [CurrentPath,filesep,fileName];
        
        clear cva calName CurrentPath
        save(savePath);
    case 'approximate'
        app = cva;
        
        natzerosapprox  =    csvread(app.Path,app.CSV(1,1),app.CSV(1,2),app.Range{1});
        seriesapprox    =    csvread(app.Path,app.CSV(2,1),app.CSV(2,2),app.Range{2});
        excessVecapprox =    csvread(app.Path,app.CSV(3,1),app.CSV(3,2),app.Range{3});
        
        [~,appName,~] = fileparts(app.Path);
        fileNameapprox = [appName,'.approx'];
        [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
        savePathapprox = [CurrentPath,filesep,fileNameapprox];
        
        clear cva appName CurrentPath
        save(savePathapprox);
        savePath = savePathapprox;
end

% --- Executes on button press in indirect.
function indirect_Callback(hObject, eventdata, handles)
% hObject    handle to indirect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of indirect
if get(hObject,'Value') == 1
    set(handles.direct,'Value',0);
else
    set(handles.direct,'Value',1);
end

% --- Executes on button press in direct.
function direct_Callback(hObject, eventdata, handles)
% hObject    handle to direct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of direct
if get(hObject,'Value') == 1
    set(handles.indirect,'Value',0);
else
    set(handles.indirect,'Value',1);
end

% --- Executes on button press in res_FLAGcheck.
function res_FLAGcheck_Callback(hObject, eventdata, handles)
% hObject    handle to res_FLAGcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of res_FLAGcheck

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



function coeffPath_Callback(hObject, eventdata, handles)
% hObject    handle to coeffPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coeffPath as text
%        str2double(get(hObject,'String')) returns contents of coeffPath as a double


% --- Executes during object creation, after setting all properties.
function coeffPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coeffPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in coeffFind.
function coeffFind_Callback(hObject, eventdata, handles)
% hObject    handle to coeffFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.csv');
if FileName ~= 0
    [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
    CurrentPath = [CurrentPath,filesep];
    if strcmp(CurrentPath,PathName)
        FullPath = FileName;
    else
        FullPath = [PathName,FileName];
    end
    set(handles.coeffPath,'String',FullPath)
end



function grbfcoeffPath_Callback(hObject, eventdata, handles)
% hObject    handle to grbfcoeffPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grbfcoeffPath as text
%        str2double(get(hObject,'String')) returns contents of grbfcoeffPath as a double


% --- Executes during object creation, after setting all properties.
function grbfcoeffPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grbfcoeffPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in grbfcoeffFind.
function grbfcoeffFind_Callback(hObject, eventdata, handles)
% hObject    handle to grbfcoeffFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.csv');
if FileName ~= 0
    [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
    CurrentPath = [CurrentPath,filesep];
    if strcmp(CurrentPath,PathName)
        FullPath = FileName;
    else
        FullPath = [PathName,FileName];
    end
    set(handles.grbfcoeffPath,'String',FullPath)
end


function grbfwPath_Callback(hObject, eventdata, handles)
% hObject    handle to grbfwPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grbfwPath as text
%        str2double(get(hObject,'String')) returns contents of grbfwPath as a double


% --- Executes during object creation, after setting all properties.
function grbfwPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grbfwPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in grbfwFind.
function grbfwFind_Callback(hObject, eventdata, handles)
% hObject    handle to grbfwFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.csv');
if FileName ~= 0
    [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
    CurrentPath = [CurrentPath,filesep];
    if strcmp(CurrentPath,PathName)
        FullPath = FileName;
    else
        FullPath = [PathName,FileName];
    end
    set(handles.grbfwPath,'String',FullPath)
end


function grbfcPath_Callback(hObject, eventdata, handles)
% hObject    handle to grbfcPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grbfcPath as text
%        str2double(get(hObject,'String')) returns contents of grbfcPath as a double


% --- Executes during object creation, after setting all properties.
function grbfcPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grbfcPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in grbfcFind.
function grbfcFind_Callback(hObject, eventdata, handles)
% hObject    handle to grbfcFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.csv');
if FileName ~= 0
    [CurrentPath,~,~] = fileparts(mfilename('fullpath'));
    CurrentPath = [CurrentPath,filesep];
    if strcmp(CurrentPath,PathName)
        FullPath = FileName;
    else
        FullPath = [PathName,FileName];
    end
    set(handles.grbfcPath,'String',FullPath)
end
