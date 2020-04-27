%GUI created by Matlab 'GUIDE' program.  GUI generated is a popup window
%used for custom term selection.  Called by main AOX_GUI.m

function varargout = termSelect_GUI(varargin)
% TERMSELECT_GUI MATLAB code for termSelect_GUI.fig
%      TERMSELECT_GUI, by itself, creates a new TERMSELECT_GUI or raises the existing
%      singleton*.
%
%      H = TERMSELECT_GUI returns the handle to a new TERMSELECT_GUI or the handle to
%      the existing singleton*.
%
%      TERMSELECT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TERMSELECT_GUI.M with the given input arguments.
%
%      TERMSELECT_GUI('Property','Value',...) creates a new TERMSELECT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before termSelect_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to termSelect_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help termSelect_GUI

% Last Modified by GUIDE v2.5 04-Feb-2020 13:50:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @termSelect_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @termSelect_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    if contains(varargin{1}," ")==0
    gui_State.gui_Callback = str2func(varargin{1});
    end
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before termSelect_GUI is made visible.
function termSelect_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to termSelect_GUI (see VARARGIN)

% Choose default command line output for termSelect_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Terms are listed in following order:
%  F, |F|, F*F, F*|F|, F*G, |F*G|, F*|G|, |F|*G, F*F*F, |F*F*F|, F*G*G, F*G*H
termList={' F, ',' |F|, ', ' F*F, ', ' F*|F|, ', ' F*G, ', ' |F*G|, ', ' F*|G|, ', ' |F|*G, ', ' F*F*F, ', ' |F*F*F|, ',' F*G*G, ',' F*G*H '};
% Update terms with defaults
if isempty(varargin)==0
    defaultTerms=varargin{1};
    handles.T1.Value=contains(defaultTerms,termList{1});
    handles.T2.Value=contains(defaultTerms,termList{2});
    handles.T3.Value=contains(defaultTerms,termList{3});
    handles.T4.Value=contains(defaultTerms,termList{4});
    handles.T5.Value=contains(defaultTerms,termList{5});
    handles.T6.Value=contains(defaultTerms,termList{6});
    handles.T7.Value=contains(defaultTerms,termList{7});
    handles.T8.Value=contains(defaultTerms,termList{8});
    handles.T9.Value=contains(defaultTerms,termList{9});
    handles.T10.Value=contains(defaultTerms,termList{10});
    handles.T11.Value=contains(defaultTerms,termList{11});
    handles.T12.Value=contains(defaultTerms,termList{12});
end

% UIWAIT makes termSelect_GUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = termSelect_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

termList={" F, "," |F|, ", " F*F, ", " F*|F|, ", " F*G, ", " |F*G|, ", " F*|G|, ", " |F|*G, ", " F*F*F, ", " |F*F*F|, "," F*G*G, "," F*G*H "};
termInclude=zeros(10,1);
termInclude(1)=handles.T1.Value;
termInclude(2)=handles.T2.Value;
termInclude(3)=handles.T3.Value;
termInclude(4)=handles.T4.Value;
termInclude(5)=handles.T5.Value;
termInclude(6)=handles.T6.Value;
termInclude(7)=handles.T7.Value;
termInclude(8)=handles.T8.Value;
termInclude(9)=handles.T9.Value;
termInclude(10)=handles.T10.Value;
termInclude(11)=handles.T11.Value;
termInclude(12)=handles.T12.Value;

if any(termInclude)
    outStruct.termInclude = strcat(termList{logical(termInclude)});  
else
    outStruct.termInclude='Click to select terms';
end
uiresume(handles.figure1);
handles.output = outStruct;
guidata(hObject, handles);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes on button press in T1.
function T1_Callback(hObject, eventdata, handles)
% hObject    handle to T1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T1


% --- Executes on button press in T2.
function T2_Callback(hObject, eventdata, handles)
% hObject    handle to T2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T2


% --- Executes on button press in T3.
function T3_Callback(hObject, eventdata, handles)
% hObject    handle to T3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T3


% --- Executes on button press in T4.
function T4_Callback(hObject, eventdata, handles)
% hObject    handle to T4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T4


% --- Executes on button press in T5.
function T5_Callback(hObject, eventdata, handles)
% hObject    handle to T5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T5


% --- Executes on button press in T6.
function T6_Callback(hObject, eventdata, handles)
% hObject    handle to T6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T6


% --- Executes on button press in T7.
function T7_Callback(hObject, eventdata, handles)
% hObject    handle to T7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T7


% --- Executes on button press in T8.
function T8_Callback(hObject, eventdata, handles)
% hObject    handle to T8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T8


% --- Executes on button press in T9.
function T9_Callback(hObject, eventdata, handles)
% hObject    handle to T9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T9


% --- Executes on button press in T10.
function T10_Callback(hObject, eventdata, handles)
% hObject    handle to T10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T10


% --- Executes on button press in T11.
function T11_Callback(hObject, eventdata, handles)
% hObject    handle to T11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T11


% --- Executes on button press in T12.
function T12_Callback(hObject, eventdata, handles)
% hObject    handle to T12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T12
