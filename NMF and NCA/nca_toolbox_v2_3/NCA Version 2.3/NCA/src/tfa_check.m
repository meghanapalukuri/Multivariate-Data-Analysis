function varargout = tfa_check(varargin)
% TFA_CHECK M-file for tfa_check.fig
%      TFA_CHECK, by itself, creates a new TFA_CHECK or raises the existing
%      singleton*.
%
%      H = TFA_CHECK returns the handle to a new TFA_CHECK or the handle to
%      the existing singleton*.
%
%      TFA_CHECK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TFA_CHECK.M with the given input arguments.
%
%      TFA_CHECK('Property','Value',...) creates a new TFA_CHECK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tfa_check_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tfa_check_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help tfa_check

% Last Modified by GUIDE v2.5 05-Oct-2005 16:19:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tfa_check_OpeningFcn, ...
                   'gui_OutputFcn',  @tfa_check_OutputFcn, ...
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

% --- Executes just before tfa_check is made visible.
function tfa_check_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tfa_check (see VARARGIN)

% Choose default command line output for tfa_check
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using tfa_check.
if strcmp(get(hObject,'Visible'),'off')
    %%
    %% load up the Plot for the TFA matrix.
    PlotLabels=evalin('base','PlotLabels');
    P=evalin('base','P');
    set(handles.popupmenu1,'String',PlotLabels);
    
    % create a function to do this for us, and we will just call the
    % function...
    idx_selected=1;
    plot(P(idx_selected,:),'linewidth',2);
   
    gene_ids=evalin('base','gene_names_final');
    Af=evalin('base','Ap');
    
    gidx=find(Af(:,idx_selected)~=0);
    genes=gene_ids(gidx);
    cs=Af(gidx,idx_selected);
    for n=1:length(cs),
        sCS{n} = cs(n);
    end
    
   set(handles.listbox2,'String',genes);
   set(handles.listbox1,'String',cs);
end

% UIWAIT makes tfa_check wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tfa_check_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;




% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

P=evalin('base','P');
idx_selected=get(hObject,'Value');
plot(P(idx_selected,:),'linewidth',2);
gene_ids=evalin('base','gene_names_final');
Af=evalin('base','Ap');
gidx=find(Af(:,idx_selected)~=0);
genes=gene_ids(gidx);
cs=Af(gidx,idx_selected);
for n=1:length(cs),
    sCS{n} = cs(n);
end
   
set(handles.listbox2,'String',genes);
set(handles.listbox1,'String',cs);
set(handles.checkbox2,'Value',0);



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idx = get(handles.popupmenu1,'Value');
if idx>1,
    idx=idx-1;
end
set(handles.popupmenu1,'Value',idx);
P=evalin('base','P');
idx_selected=get(handles.popupmenu1,'Value');
plot(P(idx_selected,:),'linewidth',2);
gene_ids=evalin('base','gene_names_final');
Af=evalin('base','Ap');
gidx=find(Af(:,idx_selected)~=0);
genes=gene_ids(gidx);
cs=Af(gidx,idx_selected);
for n=1:length(cs),
    sCS{n} = cs(n);
end
   
set(handles.listbox2,'String',genes);
set(handles.listbox1,'String',cs);
set(handles.checkbox2,'Value',0);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = get(handles.popupmenu1,'Value');
totalTFAs = evalin('base','size(P,1);');
if idx~=totalTFAs,
    idx=idx+1;
end
set(handles.popupmenu1,'Value',idx);
P=evalin('base','P');
idx_selected=get(handles.popupmenu1,'Value');
plot(P(idx_selected,:),'linewidth',2);
gene_ids=evalin('base','gene_names_final');
Af=evalin('base','Ap');
gidx=find(Af(:,idx_selected)~=0);
genes=gene_ids(gidx);
cs=Af(gidx,idx_selected);
for n=1:length(cs),
    sCS{n} = cs(n);
end
   
set(handles.listbox2,'String',genes);
set(handles.listbox1,'String',cs);
set(handles.checkbox2,'Value',0);

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

idx=get(hObject,'Value');
set(handles.listbox2,'Value',idx);


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


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2

idx=get(hObject,'Value');
set(handles.listbox1,'Value',idx);

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox2

idx = get(handles.popupmenu1,'Value');
P = evalin('base','P');
A = evalin('base','Ap');

P(idx,:)=-1*P(idx,:);
A(:,idx)=-1*A(:,idx);

assignin('base','P',P);
assignin('base','Ap',A);

P=evalin('base','P');
idx_selected=get(handles.popupmenu1,'Value');
plot(P(idx_selected,:),'linewidth',2);
gene_ids=evalin('base','gene_names_final');
A=evalin('base','Ap');
gidx=find(A(:,idx_selected)~=0);
genes=gene_ids(gidx);
cs=A(gidx,idx_selected);
for n=1:length(cs),
    sCS{n} = cs(n);
end
   
set(handles.listbox2,'String',genes);
set(handles.listbox1,'String',cs);



% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% close ourself.
close;


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('The TFA Analysis Utility allows you to flip the TFA and CS values by multiplying by negative 1.  This is used in cases where you know the sign of the TFA experimentally or for a specific gene you know the sign of the control strength.  For example, activation/repression.'); 

