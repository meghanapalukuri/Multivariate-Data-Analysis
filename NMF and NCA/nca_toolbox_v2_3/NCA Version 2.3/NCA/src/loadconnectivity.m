function varargout = loadconnectivity(varargin)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% LOADCONNECTIVITY M-file for loadconnectivity.fig
%      LOADCONNECTIVITY, by itself, creates a new LOADCONNECTIVITY or raises the existing
%      singleton*.
%
%      H = LOADCONNECTIVITY returns the handle to a new LOADCONNECTIVITY or the handle to
%      the existing singleton*.
%
%      LOADCONNECTIVITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADCONNECTIVITY.M with the given input arguments.
%
%      LOADCONNECTIVITY('Property','Value',...) creates a new LOADCONNECTIVITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before loadconnectivity_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to loadconnectivity_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help loadconnectivity

% Last Modified by GUIDE v2.5 12-Aug-2005 12:19:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @loadconnectivity_OpeningFcn, ...
                   'gui_OutputFcn',  @loadconnectivity_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before loadconnectivity is made visible.
function loadconnectivity_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to loadconnectivity (see VARARGIN)

% Choose default command line output for loadconnectivity
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes loadconnectivity wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = loadconnectivity_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

os=computer;
if strcmpi(os,'GLNX86'),
  a=dir('connectivity_db/*.mat');
elseif strcmpi(os,'MAC'),
  a=dir('connectivity_db/*.mat');      
else,
  a=dir('connectivity_db\*.mat');
end
m='';
for i=1:length(a),
    m=strvcat(m,(a(i).name));
%    j=[j;m];
end
m=cellstr(m);
%j=cellstr(j);
set(hObject,'String',m);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guihandles;
contents = get(handles.listbox1,'String');
filename = contents{get(handles.listbox1,'Value')};
os=computer;
%if strcmpi(os,'GLNX86'),
  h=msgbox(strcat('Please wait while the program loads file: ',filename));
  W=load(strcat('connectivity_db/',filename));
  %else,
 % W=load(strcat('connectivity_db\',filename));
  %end

assignin('base','Plabels',W.Plabels);
assignin('base','A',W.A);
assignin('base','Alabels',W.Alabels);
close(h); 
msgbox('Successfully loaded connectivity data. Please click OK to return to NCA Toolbox.');
h=evalin('base','h');
close(h);
evalin('base','clear h');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

status_ok=0;
disp('Select Connectivity Data File');
[filename,pathname]=uigetfile('*.txt','Open input data file');
if filename==0,
    return;
end
file = strcat(pathname,filename);
h=msgbox(strcat('Please wait while the program loads file: ',filename));
try,
    [A,Alabels,Plabels]=read_A(file);
catch,
    close(h);
    msgbox('Error Reading File');
    return;
end
    close(h);
    go=0;
    j=length(find(isnan(A)));
    k=length(find(isinf(A)));
    if j==0 & k==0,
        go=1;
    end

if length(Alabels)==0,
    go=0;
end

if go==1,
  assignin('base','A',A);
  assignin('base','Alabels',Alabels);
  assignin('base','Plabels',Plabels);  
  status_ok=1;
end
% check size of files to insure correctness of load.

if status_ok,
    msgbox('Successfully loaded Connectivity data.  Please click OK to return to NCA Toolbox.');
    h=evalin('base','h');
    close(h);
    evalin('base','clear h');
elseif go==0,
  msgbox('There is missing data present in the connectivity matrix. Please fix your file before loading into the NCA Toolbox.');
  uiwait;
else,
    warndlg(strcat('Error loading from data file',filename));
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guihandles;
contents = get(handles.listbox1,'String');
filename = contents{get(handles.listbox1,'Value')};
os=computer;
try, 
  W=load(strcat('connectivity_db/',filename),'AboutConnectivity'); 
  msgbox(W.AboutConnectivity);
catch,
  msgbox('Could not find a description of the connectivity file you selected');
end