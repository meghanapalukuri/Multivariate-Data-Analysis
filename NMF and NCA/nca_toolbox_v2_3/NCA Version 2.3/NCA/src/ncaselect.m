function varargout = ncaselect(varargin)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% NCASELECT M-file for ncaselect.fig
%      NCASELECT, by itself, creates a new NCASELECT or raises the existing
%      singleton*.
%
%      H = NCASELECT returns the handle to a new NCASELECT or the handle to
%      the existing singleton*.
%
%      NCASELECT('Property','Value',...) creates a new NCASELECT using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ncaselect_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      NCASELECT('CALLBACK') and NCASELECT('CALLBACK',hObject,...) call the
%      local function named CALLBACK in NCASELECT.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ncaselect

% Last Modified by GUIDE v2.5 23-Aug-2005 19:07:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ncaselect_OpeningFcn, ...
                   'gui_OutputFcn',  @ncaselect_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ncaselect is made visible.
function ncaselect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ncaselect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ncaselect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ncaselect_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


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


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles is empty.
  handles = guihandles;
  index_selected = get(handles.listbox1,'Value');
  data = get(handles.listbox1, 'String' );
  selected_tf = data{index_selected};

  % now add selected TFs to list, appending any already in listbox.
  me = char(get(handles.listbox2,'String'));
  me = union(selected_tf, me, 'rows'); 
  me=cellstr(me);
 
  
set(handles.listbox2,'String',me);
set(handles.text2,'String',length(me));
set(handles.pushbutton5,'Enable','off');

%update guidata to store latest changes of A_final
data =guidata(handles.listbox1);
data1 = guidata(handles.pushbutton4);
data1.A_final = data.A0;
%data1.tfa_ids_final = data.tfa_ids{data1.tfa_ids_final};
guidata(handles.pushbutton4,data1);

%
%--- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guihandles;
index_selected = get(handles.listbox2,'Value');
data = get(handles.listbox2,'String');
try,
    remove = data{index_selected};
    me = setdiff(char(data), remove, 'rows');
    if length(me)>0,
        data = cellstr(me);
    else,
        data = {};
    end
    set(handles.listbox2,'String',data);
    set(handles.listbox2,'Value',1);
    set(handles.text2,'String',length(data));
    set(handles.pushbutton5,'Enable','off');
end
    
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guihandles;
%M = length(EE);
% need to retrieve data from GUI model
data = guidata(handles.listbox1);
A0 = data.A0;
tfa_ids = data.tfa_ids;

val = get(handles.popupmenu1,'Value');
string_list = get(handles.popupmenu1,'String');
M = str2double(string_list{val}); % string_list{val};
gui_active(1);
f=progressbar([],0,'Searching for an NCA Compliant Network. Please wait.');
id=0;
u=[];
%M=M+1;
%for i=1:M,
%    j=M-i;
    progressbar(f,i/M);
    in_S.data=A0;
    in_S.dataE=evalin('base','E0');
    in_S.genename=evalin('base','gene_ids');
    in_S.tfname=tfa_ids;
    try,
        in_S.exptname=evaline('base','Experiment_labels');
    catch
        in_S.exptname=[];
    end
    out_S=Subnetwork(in_S,2);
    z=length(out_S.tfname);
    x=randperm(z);
    if M<z,
        out_S=pickNetwork(out_S,x(1:M));
    else
        out_S=pickNetwork(out_S,[1:z]);
    end
     
%end
progressbar(f,-1);
me = out_S.tfname;
set(handles.listbox2,'String',me);
set(handles.text2,'String',length(me));


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% this button checks NCA compliance, and saves the current matrix if
% it is compliant.


% just read from the listbox2.  Then read tfa_ids original and then store
% to tfa_ids_final .

handles = guihandles;
data = guidata(handles.listbox1);

selected_tfs = cellstr(get(handles.listbox2,'String'));

u=[];
for i=1:length(selected_tfs),
   for j=1:length(data.tfas),
    % get the index for each tf
       if strcmp(data.tfas(j),selected_tfs(i)),
         u = [u,j];
       end
   end
end


id = 0;
%u=data.tfa_ids_final;
[id,A,gene_index,add,remove,add_genes,remove_genes] = RankEval(data.A0,u);
% now check id and if identifiable update GUI and save experiment status.
% if not id then pop up a dialog box suggesting which genes to remove or
% add.
if (1==id),
        data.tfa_ids_final=u;
        data.A_final = A;
        gene_ids=evalin('base','gene_ids');
        data.gene_names = gene_ids(gene_index);        
        set(handles.pushbutton5,'Enable','on');
        guidata(handles.pushbutton4,data);    
        
else
        warndlg('Matrix is not identifiable with the TFs you have selected.  Auto Adjusting List...');
        u = union(u,add)
        u = setdiff(u,remove)
        set(handles.listbox2,'String',sort(data.tfas(u)))
        
        % incorporate suggestions? yes. 9/25/2004 SJG
        % show add list and remove list to get identifiable matrix.

end
    
    


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% this button saves the A matrix and the TFA list into SBML format
% including a source description

handles = guihandles;
data = guidata(handles.pushbutton4);
gene_names = data.gene_names;
guidata(handles.pushbutton4,data);
assignin('base','gene_names_final',gene_names);
assignin('base','tfa_ids_final',data.tfa_ids_final);
assignin('base','A_rank_final',data.A_final);
evalin('base','clear Pconstraints;clear Exp_constraints;');
msgbox('Saved the NCA Compliant Network to the Global Workspace.  You may now run NCA, or you may specify gNCA TF constraints and then run NCA.');
fprintf('Saved the NCA Compliant Network to the Global Workspace.\nYou may now run NCA, or you may specify gNCA TF constraints and then run NCA.\n');

h=evalin('base','h');
close(h);
evalin('base','clear h');

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guihandles;
set(handles.listbox2,'String',{});
set(handles.text2,'String',0);
data = guidata(handles.listbox1);
data1 = guidata(handles.pushbutton4);
data1.A_final = data.A0;
guidata(handles.pushbutton4,data1);
set(handles.pushbutton5,'Enable','off');


% --- Executes on button press in maxnetwork.
function maxnetwork_Callback(hObject, eventdata, handles)
% hObject    handle to maxnetwork (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guihandles;

% test if spectrum or Not.

selected_tfs = cellstr(get(handles.listbox2,'String'));
in_S.dataE=evalin('base','E0');
tfname=evalin('base','tfa_ids');
tfs=[];
for i=1:length(selected_tfs),
    j=find(strcmpi(selected_tfs{i},tfname));
        if j>0, 
            tfs = [tfs,j];
        end    
end

in_S.tfname = tfname;
in_S.data=evalin('base','A0');
in_S.genename=evalin('base','gene_ids');
in_S.exptname=evalin('base','Experiment_labels');
[new_S,Af] = pickNetwork(in_S,tfs);
out_S = Subnetwork(new_S,2);
my_tfas=[];
try,
    for i=1:length(selected_tfs),
    idx = find(strcmpi(selected_tfs{i},out_S.tfname));
    if idx>0,
        my_tfas=[my_tfas;idx];
    end
end
catch,
end

if length(my_tfas)<1,
    msgbox('Cannot find any NCA Compliant subnetwork from your selected regulators. Please adjust your TF selection.'); 
else,
    successStr=sprintf('Congratulations!  The network is NCA compliant for %d TFs', length(out_S.tfname));
    msgbox(successStr);
    
    set(handles.listbox2,'String',out_S.tfname);
    set(handles.text2,'String',length(out_S.tfname));
    data = guidata(handles.pushbutton4);
    data.A_final=out_S.data;
    data.gene_names=out_S.genename;
    assignin('base','Ef',out_S.dataE);
    tfs=[];
    for i=1:length(out_S.tfname),
        j=find(strcmpi(out_S.tfname{i},tfname));
        if j>0, 
            tfs = [tfs,j];
        end
    end
    data.tfa_ids_final=tfs;
    set(handles.pushbutton5,'Enable','on');
    guidata(handles.pushbutton4,data);
end




