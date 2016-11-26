function varargout = multi_plot(varargin)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% MULTI_PLOT M-file for multi_plot.fig
%      MULTI_PLOT, by itself, creates a new MULTI_PLOT or raises the existing
%      singleton*.
%
%      H = MULTI_PLOT returns the handle to a new MULTI_PLOT or the handle to
%      the existing singleton*.
%
%      MULTI_PLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTI_PLOT.M with the given input arguments.
%
%      MULTI_PLOT('Property','Value',...) creates a new MULTI_PLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multi_plot_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multi_plot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help multi_plot

% Last Modified by GUIDE v2.5 09-Nov-2004 12:18:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multi_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @multi_plot_OutputFcn, ...
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


% --- Executes just before multi_plot is made visible.
function multi_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multi_plot (see VARARGIN)

% Choose default command line output for multi_plot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes multi_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = multi_plot_OutputFcn(hObject, eventdata, handles)
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

PlotLabels = evalin('base','PlotLabels');
set(hObject,'String',PlotLabels);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% get from the global workspace the Plabel matrix.




% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check that user has made a selection.
% update the plot to include any selected TFAs
% plot as different colors.

% get from the global workspace the Plabels and P matrix

selected_ind=get(handles.listbox1,'Value');
P = evalin('base','P');
P_dev = evalin('base','P_dev');
Plabels = evalin('base','PlotLabels');

plotSD=evalin('base','plotSD');

%plot(P(selected_ind,:));
%cla;
% make up to 4 subplots, or 2 or 1 or 3 if necessary.
 numplot = length(selected_ind);
% if numplot>1,

% create a new figure

organism=evalin('base','organism');
plot_fig = figure;
  for i=1:numplot,  
    subplot(numplot,1,i);
    try,
        z=[];
        Exp_labels=evalin('base','Experiment_labels');
        for j=1:length(Exp_labels),
            z=[z,str2double(Exp_labels(j))];
            if length(find(isnan(z)))>0,
                error('Cannot convert Experiment labels to time point');
            end
        end
        if plotSD,
            shade(z/60,P(selected_ind(i),:),P(selected_ind(i),:)-P_dev(selected_ind(i),:),P(selected_ind(i),:)+P_dev(selected_ind(i),:),1)
            xlabel('Time (Hour)')

        else,
            plot(z/60,P(selected_ind(i),:));
            xlabel('Time (Hour)')
        end
    catch, 
    if plotSD,
        shade([1:size(P,2)],P(selected_ind(i),:),P(selected_ind(i),:)-P_dev(selected_ind(i),:),P(selected_ind(i),:)+P_dev(selected_ind(i),:),1)
    else,
        plot(P(selected_ind(i),:));
    end    
    xlabel('Experiment');
    end
    if strcmpi(organism,'Spectra'),  
        ylabel('Absorption');
        xlabel('Wavelength');
        try,
            Exp_labels=evalin('base','Experiment_labels');
            set(gca,'XTickLabel',Exp_labels);
        end;
    else,
        
        ylabel('TFA Profile');
    end
    title(Plabels(selected_ind(i)));
    
  end  
%else,
%  handles.axes1=plot(P(selected_ind,:));
%end
%axes(handles.axes1);
