function varargout = tf_experiment(varargin)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tf_experiment_OpeningFcn, ...
                   'gui_OutputFcn',  @tf_experiment_OutputFcn, ...
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


% --- Executes just before tf_experiment is made visible.
function tf_experiment_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tf_experiment wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tf_experiment_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

try,
    exp_labels=evalin('base','exp_labels');    
    E0=evalin('base','E0');
    exp_labels=1:size(E0,2);    
    assignin('base','exp_labels',exp_labels);  
catch,
    % do nothing
end

set(hObject,'String',exp_labels);



% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

handles=guihandles;
contents = get(handles.listbox2,'String');
contents=cellstr(contents);
tf_selected = contents{get(handles.listbox2,'Value')};
contents = get(handles.listbox4,'String');
contents=cellstr(contents);
exp_selected = contents{get(handles.listbox4,'Value')};

% now pair the strings to together and record in a structure the
% constraints and save this to the base workspace.
%

b = strcat(tf_selected, '__');
b = strcat(b,exp_selected);
b = cellstr(b);
b=union(b,get(handles.listbox5,'String'));
set(handles.listbox5,'String',b);

%%


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guihandles;
contents = get(handles.listbox5,'String');
contents = cellstr(contents);
selected = contents{get(handles.listbox5,'Value')};
contents = setdiff(contents,selected);
set(handles.listbox5,'String',contents);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

handles=guihandles;
contents=get(handles.listbox5,'String');
c=cell(length(contents),2);
for i=1:length(contents),
    c(i,:)=strsplit('__',contents{i});
end
Exp_constraints = c(:,1);
Pconstraints = c(:,2);

assignin('base','Exp_constraints',Exp_constraints);
assignin('base','Pconstraints',Pconstraints);
h=evalin('base','h');
close(h);
msgbox('Saved Constraints to Global Workspace');

% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
try,
  tfa_ids_final = evalin('base','tfa_ids_final');
catch,
  warndlg('Error: Please Match Prune and Select TFs before setting P constraints');
end

try,
  tfa_ids = evalin('base','tfa_ids');
catch,
  warndlg('Could not find tfa_ids in base workspace.');    
end
set(hObject,'String',tfa_ids(tfa_ids_final));


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

set(hObject,'String',{});

% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
%empty

