function varargout = avg_contrib_by_tf(varargin)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% AVG_CONTRIB_BY_TF M-file for avg_contrib_by_tf.fig
%      AVG_CONTRIB_BY_TF, by itself, creates a new AVG_CONTRIB_BY_TF or raises the existing
%      singleton*.
%
%      H = AVG_CONTRIB_BY_TF returns the handle to a new AVG_CONTRIB_BY_TF or the handle to
%      the existing singleton*.
%
%      AVG_CONTRIB_BY_TF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AVG_CONTRIB_BY_TF.M with the given input arguments.
%
%      AVG_CONTRIB_BY_TF('Property','Value',...) creates a new AVG_CONTRIB_BY_TF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before avg_contrib_by_tf_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to avg_contrib_by_tf_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help avg_contrib_by_tf

% Last Modified by GUIDE v2.5 26-Jul-2005 13:42:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @avg_contrib_by_tf_OpeningFcn, ...
                   'gui_OutputFcn',  @avg_contrib_by_tf_OutputFcn, ...
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

% --- Executes just before avg_contrib_by_tf is made visible.
function avg_contrib_by_tf_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to avg_contrib_by_tf (see VARARGIN)

% Choose default command line output for avg_contrib_by_tf
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% This sets up the initial plot - only do when we are invisible
% so window can get raised using avg_contrib_by_tf.
if strcmp(get(hObject,'Visible'),'off')
end

% UIWAIT makes avg_contrib_by_tf wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = avg_contrib_by_tf_OutputFcn(hObject, eventdata, handles)
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



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

gene_names_final=evalin('base','gene_names_final');
set(hObject, 'String', gene_names_final);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guihandles;
popup_sel_index = get(handles.popupmenu2, 'Value');
% get selected string.  Then use to create our plot.
gene_names=evalin('base','gene_names_final');
selected_gene=gene_names(popup_sel_index);
A=evalin('base','Ap');
P=evalin('base','P');
E=evalin('base','Ep');
Plabels=evalin('base','PlotLabels');
idx = find(A(popup_sel_index,:)~=0);
figure;
Efit=A(popup_sel_index,idx)*P(idx,:);
plot(Efit,'b-+');
hold on;
legendstr{1} = strcat('Fitted Gene Expression: ',char(selected_gene));
x=1;
y=1;
z=1;
colorstr={'g';'c';'r';'y';'m';'k'};
for i=1:length(idx),
    if i<6,
        plot(A(popup_sel_index,idx(i))*P(idx(i),:),colorstr{i});
        legendstr{i+1} = strcat(char(Plabels(idx(i))),' Contribution');
        y=y-0.1;
    end
end
xlabel('Experiment Number');
ylabel('log mRNA ratio');
title(strcat('Contribution of Each TF for Gene: ',char(selected_gene)));
legend(legendstr);
