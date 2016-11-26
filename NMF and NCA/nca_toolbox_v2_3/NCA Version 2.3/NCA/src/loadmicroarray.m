function varargout = loadmicroarray(varargin)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% LOADMICROARRAY M-file for loadmicroarray.fig
%      LOADMICROARRAY, by itself, creates a new LOADMICROARRAY or raises the existing
%      singleton*.
%
%      H = LOADMICROARRAY returns the handle to a new LOADMICROARRAY or the handle to
%      the existing singleton*.
%
%      LOADMICROARRAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADMICROARRAY.M with the given input arguments.
%
%      LOADMICROARRAY('Property','Value',...) creates a new LOADMICROARRAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before loadmicroarray_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to loadmicroarray_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help loadmicroarray

% Last Modified by GUIDE v2.5 11-Nov-2004 11:46:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @loadmicroarray_OpeningFcn, ...
                   'gui_OutputFcn',  @loadmicroarray_OutputFcn, ...
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


% --- Executes just before loadmicroarray is made visible.
function loadmicroarray_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to loadmicroarray (see VARARGIN)

% Choose default command line output for loadmicroarray
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes loadmicroarray wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = loadmicroarray_OutputFcn(hObject, eventdata, handles)
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

% search in db for the microarray .mat files 
% 

os=computer;
if strcmpi(os,'GLNX86'),
  a=dir('microarray_db/*.mat');
elseif strcmpi(os,'MAC'),
  a=dir('microarray_db/*.mat');      
else,
  a=dir('microarray_db\*.mat');
end

% format a for cellstr.
% construct a dynamic function to get what we want?
% this is a hack....

m='';
%j=cellstr;
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
filename = strcat('microarray_db/',filename);
h=msgbox(strcat('Please wait while the program loads file: ',filename));
W=load(filename);
assignin('base','E',W.E);
assignin('base','Elabels',W.Elabels);
assignin('base','Experiment_labels',W.Experiment_labels);
close(h);
msgbox('Successfully loaded Microarray Data. Please click OK to return to NCA Toolbox.');
h=evalin('base','h');
close(h);
evalin('base','clear h');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

status_ok=0;
disp('Select Microarray Data File');
[filename,pathname]=uigetfile({'*.txt'},'Open input data file');
if filename==0,
    return;
end
file = strcat(pathname,filename);
j=0;
k=0;
h=msgbox(strcat('Please wait while the program loads file: ',filename));
 try,
   [E,Elabels,Experiment_labels]=read_A(file);
      
 catch,
       close(h); 
       msgbox('Could not read file');
       return;
 end

go=0;
j=length(find(isnan(E)));
k=length(find(isinf(E)));

  
 if j==0 & k==0,
       go=1;
       if length(Elabels)==0,
           go=0;
       end
   end
   if 1==go,
       status_ok=1;
       assignin('base','E',E);
       assignin('base','Elabels',Elabels);
       assignin('base','Experiment_labels',Experiment_labels);
   else,
       warndlg('There are missing data points in the data matrix.  Please fix your file before loading into the NCA Toolbox.');
       status_ok=0;
       close(h);
       return;
   end

if status_ok,
    msgbox('Successfully loaded Microarray data from File. Please click OK to return to NCA Toolbox.');
    close(h);
    h=evalin('base','h');
    close(h);
    evalin('base','clear h');
    return;
elseif go==0,
    % skip
else,
    warndlg('Unspecified error loading data from file.  Please check that your data file is in the correct row/column format -- see README.TXT).');
end


function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guihandles;
contents = get(handles.listbox1,'String');
filename = contents{get(handles.listbox1,'Value')};
os=computer;
try,
  warning off;  
  W=load(strcat('microarray_db/',filename),'AboutMicroarray'); 
  msgbox(W.AboutMicroarray);
catch,
  msgbox('Could not find a description of the microarray file you selected');
end