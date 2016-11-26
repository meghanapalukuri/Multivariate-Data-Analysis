function varargout = nca_menu(varargin)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% nca_menu is a graphical user interface (GUI) to set up
% an NCA analysis.  It allows the user to load in Connectivity
% and Expression (or Signal) data, match the two matrices by row
% and then select a set (or subset) of regulators, find the
% maximum NCA compliant network, and then run NCA.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nca_menu_OpeningFcn, ...
                   'gui_OutputFcn',  @nca_menu_OutputFcn, ...
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


% --- Executes just before nca_menu is made visible.
function nca_menu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nca_menu (see VARARGIN)

% Choose default command line output for nca_menu
handles.output = hObject;
nca_logo = imread('nca.jpg');
image(nca_logo);
axis equal;
axis off;
% Update handles structure

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = nca_menu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function NewWorkspace_2_Callback(hObject, eventdata, handles)
% hObject    handle to NewWorkspace_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open up a dialog box

x=questdlg('Do you want to Start a new NCA Session');
if strcmpi(x,'yes'),
    evalin('base','clear');
    handles=guihandles;
    contents=get(handles.popupmenu1,'String');
    organism=contents{get(handles.popupmenu1,'Value')};
    if strcmpi(organism,'Spectra'),
        set(handles.pushbutton4,'String','Load Spectra Data');
        set(handles.pushbutton7,'String','Match Spectra Data to Connectivity');
    else,
        set(handles.pushbutton4,'String','Load Microarray Data');
        set(handles.pushbutton7,'String','Match Microarray Data to Connectivity');
    end
    assignin('base','organism',organism);
    assignin('base','iterations',10);  % default settings
    evalin('base','clc;');
end


% --------------------------------------------------------------------
function LoadWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to LoadWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName, PathName] = uigetfile('saved_workspaces/*.mat','Select the Matlab Workspace');

if isequal(0,FileName),
    % do nothing
else, 
    file=strcat(PathName,FileName);
    evalin('base','clear;');
    assignin('base','file',file);
    evalin('base','load(file)');
end

% --------------------------------------------------------------------
function SaveWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to SaveWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname]=uiputfile('saved_workspaces/*.mat');
if isequal(0,filename),
    % do nothing
else,
    file=strcat(pathname,filename);
    assignin('base','file',file);
    evalin('base','save(file)');
end

% --------------------------------------------------------------------
function ExportTFA_Callback(hObject, eventdata, handles)
% hObject    handle to ExportTFA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try,
    exp_labels=evalin('base','Experiment_labels');
catch,
    exp_labels=[];
end
try,
    PlotLabels=evalin('base','PlotLabels');
catch
    PlotLabels=[];
end
try,
    P=evalin('base','P');
    go=1;
catch,
    msgbox('Cannot find TFA matrix (variable P). Please load a workspace or run NCA');
    go=0;
    return;
end

if go==1,
    [filename,path]=uiputfile('*.txt','Enter Filename');
    if ~isequal(filename,0),
        file=strcat(path,filename);
        write_tfa_matrix(P,PlotLabels,exp_labels,file);
        msgbox(strcat('Successfully wrote TFA Matrix to file ', file));
    end
end


% --------------------------------------------------------------------
function ExitAndClose_Callback(hObject, eventdata, handles)
% hObject    handle to ExitAndClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --------------------------------------------------------------------
function TFAplot_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try,
    Plabels=evalin('base','PlotLabels'); 
    evalin('base','plotSD=0;');
    multi_plot;
catch,
  warndlg('Please load a previously saved workspace, or run an NCA Analysis.');
  return;
end


% --------------------------------------------------------------------
function Untitled_23_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try,
    
    R=evalin('base','Ep-Ap*P;');
    figure;pcolor(R);colorbar;
    title Heatmap; 
    ylabel Genes; 
    xlabel Experiment;
catch,
    warndlg('Cannot create Heatmap.  Error most likely that NCA has not been run');
end

% --------------------------------------------------------------------
function AboutMenu_Callback(hObject, eventdata, handles)
% hObject    handle to AboutMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('The NCA Toolbox and NCA (c) 2003-2005. University of California Los Angeles.   Authors:  Simon Galbraith and Linh Tran. Please see: http://www.seas.ucla.edu/~liaoj for more information');

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=loadmicroarray;
assignin('base','h',h);    
% --- Executes on button press in load_conn_button.
function load_conn_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_conn_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=loadconnectivity;
assignin('base','h',h);    


% --- Executes on button press in selectnetwork.
function selectnetwork_Callback(hObject, eventdata, handles)
% hObject    handle to selectnetwork (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try,
A0 = evalin('base','A0');
gene_ids = evalin('base','gene_ids');
tfa_ids = evalin('base','tfa_ids');
TfaSelectGui(A0,gene_ids,tfa_ids);
catch,
    msgbox('Error: Please run Match Data Matrices to get an initial network');
end
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
%Match Data Matrix callback

organism = evalin('base','organism');
if strcmpi(organism,'Spectra'),
    go=2;
else,
    go=1;
end

prematched=get(handles.prematched_Checkbox,'Value');
if prematched,
    go=2;
end

try,
    A = evalin('base','A');
    Alabels = evalin('base','Alabels');
    Plabels=evalin('base','Plabels');
    E = evalin('base','E');
    Elabels = evalin('base','Elabels');
catch,
    msgbox('Please Load Connectivity and Microarray Data first.')
    return;
end

if 1==go & size(A,1)~=length(Alabels),
    warndlg('Error: Number of rows in A does not match number of genes in Alabels.'); 
    go=0;
    return;
end
if 1==go & length(Plabels)~=size(A,2),
    warndlg('Error: Number of TFs in Plabels does not match column span of A.');
    go=0;
    return;
end
if 1==go & length(Elabels)~=size(E,1),
    warndlg('Error: Number of rows in matrix E does not equal the length of Elabels');
    go=0;
    return;
end
if go==1,
  try,
  f=progressbar([],1,'Matching Connectivity Matrix to Microarray Matrix row label identifiers, please wait...');
  [E0,A0,gene_ids,tfa_ids]=MatchPrune(E,Elabels,A,Alabels,Plabels,organism);
 
  if size(E0,1)~=0,  
      assignin('base','E0',E0);
      assignin('base','A0',A0);
      assignin('base','tfa_ids',tfa_ids);
      assignin('base','gene_ids',gene_ids);
      progressbar(f,-1);
      msgbox('Successfully Matched Connectivity and Microarray Data Together.  Results are in variables E0 and A0 of the global workspace.');
      fprintf('Successfully Matched Connectivity and Microarray Data Together.\nResults are in variables E0 and A0 of the global workspace.\n');   
      
  else,
      msgbox('Error:  Cannot Match Connectivity and Data Matrices.  It is likely that one is loaded incorrectly, or that you have selected the wrong organism option. Please double-check your settings.'); 
      progressbar(f,-1);
      return;
  end
  catch,
      msgbox('Error:  Cannot match the connectivity and expression data matrices.  It is likely that one is loaded incorrectly, or that you have selected the wrong organism option. Please double-check your settings.');
      progressbar(f,-1);
      return;
  end
elseif go==2,
  % output the matrix directly to A0 as it is already matched  
try,
    A0=evalin('base','A');
    E0=evalin('base','E');
catch,
    warndlg('Please load a data matrix and a connectivity matrix first.');
    warn('Please load a data matrix and a connectivity matrix first.');
    return;
end
  if size(A0,1)~=size(E0,1),
      msgbox('You have chcked that the data is prematched by row, but the # of rows in E and A are not equal. Please prematch the connectivity and expression data in a separate program such as Microsoft Excel and reload into the NCA Toolbox.  If this is not the case, then you may need to transpose the E matrix.  To do this then please click on the transpose checkbox.');
  else,
     evalin('base','A0=A;E0=E;tfa_ids=Plabels;gene_ids=Elabels;');
     msgbox('You have checked that the data is prematched, this assumes that the matrices A and E are already matched by row. Results are in E0 and A0.');
  end
end
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
contents = get(hObject,'String');
%handles=guihandles;
organism=contents{get(hObject,'Value')};
assignin('base','organism',organism);
z=evalin('base','length(who)');

if strcmpi(organism,'Spectra'),
    set(handles.pushbutton4,'String','Load Spectra Data');
    set(handles.pushbutton7,'String','Match Spectra Data to Connectivity');
else,
    set(handles.pushbutton4,'String','Load Microarray Data');
    set(handles.pushbutton7,'String','Match Microarray Data to Connectivity');
end


if z>2,
    msgbox('You changed an important parameter, and it appears that you have an active workspace.  Please double check your settings.');   
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

contents = get(hObject,'String');
assignin('base','organism',contents{get(hObject,'Value')});
assignin('base','iterations',10);  % default settings




% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)

handles=guihandles;
organism = evalin('base','organism');
skip_norm=0;

try,
    iterations=evalin('base','iterations');
catch,
    iterations=10;       
end

if 0==strcmpi(organism,'Spectra'),
    try,   
        Ef=evalin('base','Ef');
        Af=evalin('base','A_rank_final');
        gene_names_final = evalin('base','gene_names_final'); 
        tfa_ids_final = evalin('base','tfa_ids_final');
        tfa_ids = evalin('base','tfa_ids');
        gene_ids=evalin('base','gene_ids');
        
        if length(Ef)==0,
            warndlg('Error cannot find matched Subnetwork.  This most likely means that you have selected the wrong organism to analyze in Step 1 then you have provided connectivity or microarray data for.   Please double-check your settings and fix any errors, or start a new NCA analysis.  If you recently loaded a saved workspace, then adjust the organism selection before running NCA.'); 
            return;
        end
        assignin('base','PlotLabels',tfa_ids(tfa_ids_final));
    catch,
        warndlg('Error: Could not run NCA analysis, please reload the data and connectivity matrices.');
        return;
    end
       
    tf_c=[];
    exp_c=[];
    Pf = rand(size(Af,2),size(Ef,2));  
    try,
        Pconstraints=evalin('base','Pconstraints');
        Exp_constraints=evalin('base','Exp_constraints');   
        tfa_ids = evalin('base','tfa_ids');
        for i=1:length(Pconstraints),    
            index=find(strcmpi(tfa_ids(tfa_ids_final),Pconstraints(i)));  
            Pf(index,str2double(char(Exp_constraints(i))))=0;
        end
    end

    % get status of regularization box...
    run_regularization=0;
    run_regularization=get(handles.checkbox2,'Value');
    assignin('base','regularized',run_regularization);
    gui_active(1);
    f=[];
    f = progressbar( f,0,'Running NCA. Please wait.' );
    for k=1:iterations,
        f=progressbar(f,1/iterations);        
        fprintf('Current Initial Guess Iteration %d of %d\n',k,iterations);
        if ~gui_active,
            progressbar(f,-1);
            return;
        else,    
            if run_regularization>0,
                %fprintf('Running GNCA with Regularization\n');    
                [a_iter(:,:,k), p_iter(:,:,k)] = gncar_cc(Ef,Af,Pf);
            else,
                %fprintf('Running GNCA\n');
                if get(handles.checkbox5,'Value'),
                 [a_iter(:,:,k), p_iter(:,:,k)] = gnca_cc_v2(Ef,Af,Pf,1);   
                else,
                 [a_iter(:,:,k), p_iter(:,:,k)] = gnca_cc_v2(Ef,Af,Pf,0);
                end
            end
        end
    end
    progressbar(f,-1);
    
else,
    Af = evalin('base','A0');
    Ef = evalin('base','E0');
    try,
        Plabels=evalin('base','tfa_ids');
        tfa_ids_final = evalin('base','tfa_ids_final');
        Ptfas=Plabels(tfa_ids_final);
    catch,
        Ptfas={};
        for i=1:size(Af,2),
            Ptfas{i} = int2str(i);
        end
    end
    assignin('base','PlotLabels',Ptfas);
    prompt = {'Please enter an index vector for an experiment to use for Wavlength normalization. Click CANCEL to use default normalization or if you are analyzing microarray data.','Enter a vector of absorptions to use for normalization of each regulator. Please separate by a space.'};
    response=inputdlg(prompt,'Spectra Normalization',1,{'46 46 46','1.3989 0.2019 1.2103'},'on');
    
    if length(response)>0,
        try,
            temp=strsplit(' ',response{1});
            wavelength=[];
            for i=1:size(temp,2)
                wavelength=[wavelength,str2double(temp{i})];
            end
            temp=strsplit(' ',response{2});
            absorp=[];
            for i=1:size(temp,2),
                absorp=[absorp,str2double(temp{i})];
            end       
        catch,
            msgbox('Could not parse the wavelength index or absorption vector.  These parameters must be space separated numbers.');
            return;
        end
        try,
            f=[];
            f = progressbar( f,0,'Running NCA. Please wait.' );
            gui_active(1);
            for k=1:iterations,
                f=progressbar(f,1/iterations);
                fprintf('Current Initial Guess Iteration %d of %d\n',k,iterations);
                if ~gui_active,
                    progressbar(f,-1);
                    return;
                else,                           
                    [a_iter(:,:,k),p_iter(:,:,k)]=nca_spectrum(Ef,Af,wavelength,absorp);
                    skip_norm=1;
                end
            end
            progressbar(f,-1);
            
        catch,
            warndlg('Could not use the User Input wavelength index or absorpotion parameters... using default normalization routine.');
            uiwait;
        
            f=[];
            gui_active(1);
            f = progressbar( f,0,'Running NCA. Please wait.' );
            for k=1:iterations,
            f=progressbar(f,1/iterations);        
            fprintf('Current Initial Guess Iteration %d of %d\n',k,iterations);
                if ~gui_active,
                    progressbar(f,-1);
                    return;       
                else,                                           
                    [a_iter(:,:,k), p_iter(:,:,k)] = gnca_cc_v2(Ef,Af);
                end
            end
            progressbar(f,-1);
        end
    else,
            answer = questdlg('You have selected to use the default normalization routine for microarray data to analyze Spectra data. If this is OK click Yes to continue.');
            if 0==strcmpi(answer,'Yes'),
                return;
            end;
        
            f=[];
            gui_active(1);
            f = progressbar( f,0,'Running NCA. Please wait.' );
            for k=1:iterations,
            f=progressbar(f,1/iterations);        
            fprintf('Current Initial Guess Iteration %d of %d\n',k,iterations);
                if ~gui_active,
                    progressbar(f,-1);
                    return;       
                else,                                           
                    [a_iter(:,:,k), p_iter(:,:,k)] = gnca_cc_v2(Ef,Af);
                end
            end
            progressbar(f,-1);   
    end
end
% now check if we need to compute average solution.
% if iteration is >1 then we do, else
% just report the first solution

Ap=[];
Pp=[];
if 1==iterations,
   Pp=p_iter(:,:,1);
   Ap=a_iter(:,:,1);
   Ps=zeros(size(Pp));
   As=zeros(size(Ap));
else,
    if 0==skip_norm,
        for i=1:iterations,
            [a_iter(:,:,i),p_iter(:,:,i)]=NormP(a_iter(:,:,i),p_iter(:,:,i));       
        end
        for i=1:iterations,
            [mA(:,:,i),mP(:,:,i)]=convertSignCor(squeeze(p_iter(:,:,1)),squeeze(a_iter(:,:,i)),squeeze(p_iter(:,:,i)));
        end    
    else,
        mP=p_iter;
        mA=a_iter;
    end
    
    % get average of each run..  
    % compute std. deviation over all runs.
    
    
    % now recalculate Ap from Pp
    
    Pp=mean(mP,3);
    if iterations>1,
        Ps=std(mP,0,3);
    else,
        Ps=zeros(size(Pp));
    end

    Ap = evalin('base','A_rank_final');

    %Remove the zero columns before whitening data

    
    for n=1:size(Ef,1),    
        index_A = find(Ap(n,:) ~= 0);
        P_red= Pp(index_A,:); 
        Ap(n,index_A)= Ef(n,:)/P_red;
    end
    
end  
    

if 0==skip_norm,
    [Ap,Pp,Ps]=Normalization(Ap,Pp,Ps);
end

assignin('base','P',Pp);
assignin('base','Ap',Ap);
assignin('base','Ep',Ef);
try,
    assignin('base','P_dev',Ps);
end

try,
    close(f);
end;

fprintf('\nSuccessfully ran NCA and the solution satisfies the 3rd criterion.\nYou may now plot the TFA profiles by selecting the plot menu option.\nYou may also cluster the TFA profiles by selecting the "Analyze TFAs" menu option, as well as export (File Menu)\nthe TFA profiles to a tab-delimited file.\n');
msgbox('Successfully ran NCA.  You may now plot the TFA profiles by selecting the plot menu option. You may also cluster the TFA profiles by selecting the "Analyze TFAs" menu option.');

% --- Executes on button press in Set_Constraints.
function Set_Constraints_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Constraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try,
    E0=evalin('base','E0');
    A0=evalin('base','A0');
    A_rank_final=evalin('base','A_rank_final');
    exp_labels=1:size(E0,2);    
    assignin('base','exp_labels',exp_labels);
catch,
    return;
end

h=open('tf_experiment.fig');
assignin('base','h',h);




% --------------------------------------------------------------------
function regularized_nca_Callback(hObject, eventdata, handles)
% hObject    handle to regularized_nca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=get(hObject,'Checked');
if strcmp(x,'on'),
    set(hObject,'Checked','off');
else,
    set(hObject,'Checked','on');
end


% --------------------------------------------------------------------
function normalized_TFAs_Callback(hObject, eventdata, handles)
% hObject    handle to normalized_TFAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Ask user to normalize by A or by P matrix.

x=get(hObject,'Checked');
if strcmp(x,'on'),
    set(hObject,'Checked','off');
else,
    set(hObject,'Checked','on');
end

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --------------------------------------------------------------------
function cs_select_Callback(hObject, eventdata, handles)
% hObject    handle to cs_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=get(hObject,'Checked');
if strcmp(x,'on'),
    set(hObject,'Checked','off');
else,
    set(hObject,'Checked','on');
end

% maybe I can have a global preferences structure that I can evaluate
% this would make things easier



% --------------------------------------------------------------------
function SelectTFsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SelectTFsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try, 
    A0 = evalin('base','A0');
    gene_ids = evalin('base','gene_ids');
    tfa_ids = evalin('base','tfa_ids');
    TfaSelectGui(A0,gene_ids,tfa_ids);
catch,
    warndlg('Please load Connectivity and Microarray Data before Selecting a Regulatory Network.');
end


% --------------------------------------------------------------------
function setConstraints_Callback(hObject, eventdata, handles)
% hObject    handle to setConstraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try,
    exp_labels=evalin('base','exp_labels');    
    E0=evalin('base','E0');
    exp_labels=1:size(E0,2);    
    assignin('base','exp_labels',exp_labels);
    fig=open('tf_experiment.fig');
catch,
    msgbox('Cannot find E0 or Experiment Labels, please Load data and run MatchPrune');
end


% --------------------------------------------------------------------
function PlotTFAMenu_21_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTFAMenu_21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try,
P=evalin('base','P');
tfa_ids = evalin('base','tfa_ids_final');
Plabels = evalin('base','tfa_ids');
assignin('base','Plabels',Plabels(tfa_ids));
multi_plot;
catch,
    msgbox('Cannot find TFA profiles.  Please load a saved workspace, or run a new NCA analysis');
end


% --------------------------------------------------------------------
function histogram_residuals_Callback(hObject, eventdata, handles)
% hObject    handle to histogram_residuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try,
R=evalin('base','Ep-Ap*P;');
r=[];
for i=1:size(R,1),
  r=[r,norm(R(i,:),2)];
end
figure;hist(r);
title 'Residual Error per Gene';
xlabel 'Residual Error';
ylabel 'Number of Genes';
catch,
    warndlg('Cannot create histogram.  Most likely NCA has not been run');
end

% --------------------------------------------------------------------
function hclust_Callback(hObject, eventdata, handles)
% hObject    handle to hclust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try, 
    P=evalin('base','P');
    tfname=evalin('base','PlotLabels');
catch,
    msgbox('Please run an NCA Analysis first.');
    return;
end
try, 
    expname=evalin('base','Experiment_labels');
catch,
    expname='';
end
    
try,
%    figure;
    clustergram(P,'RowLabels',tfname,'ColumnLabels',expname);
    title('Clustergram of TFA profiles');
    return;
catch,
    msgbox('Bioinformatics Toolbox not Found, using default cluster function.');
end
try,
    Z=pdist(P);
    T=linkage(Z);
    figure;	
    [a,b,perm]=dendrogram(T,'ORIENTATION','left','labels',tfname);
    title('Dendrogram of Hierarchical Clustering of TFAs');
    xlabel('Distance');
catch,
    msgbox('Cannot find Statistics toolbox.');
    return;
end



% --------------------------------------------------------------------
function avg_contrib_by_tf_Callback(hObject, eventdata, handles)
% hObject    handle to avg_contrib_by_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try,
    P=evalin('base','P');
    fig=open('avg_contrib_by_tf.fig');
catch,
   warndlg('Cannot find an NCA Analysis.  Please load a Workspace or run an NCA Analysis.');
end
% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try,
  evalin('base','E=transpose(E);temp=Elabels;Elabels=Experiment_labels;Experiment_labels=temp;');
  msgbox('You just transposed the data matrix E. This may not be the option you want unless you are analyzing Spectra (or other source signal) data.  Uncheck the box to revert to the original matrix');
catch,
  msgbox('Error:  Cannot Transpose Data Matrix E. Please check that signal data is loaded and variable E exists.')
  uiwait;
  set(hObject,'Value',0)
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Network_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Network_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NCARunTimeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to NCARunTimeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PlotMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PlotMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function TFAAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to TFAAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function HelpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to HelpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function AverageSolutions_select_Callback(hObject, eventdata, handles)
% hObject    handle to AverageSolutions_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

itercount=inputdlg('How many initial guess iterations do you want to run?','Set the Number Initial Guess Iterations',1,{'10'});

if length(itercount)>0,
    itercount=str2double(itercount);
    if (itercount<1 | isnan(itercount)),
        itercount=10;
        fprintf('You cannot run less than 1 iteration. Setting Number of Iterations to the default value of 10.\n');
    end
    assignin('base','iterations',itercount);
    
end

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2




% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NewWorkspace_2_Callback;



% --------------------------------------------------------------------
function Tutorial_Callback(hObject, eventdata, handles)
% hObject    handle to Tutorial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('tutorial.htm');



% --------------------------------------------------------------------
function ExportCS_Callback(hObject, eventdata, handles)
% hObject    handle to ExportCS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try,
    gene_ids=evalin('base','gene_names_final');
catch,
    gene_ids=[];
end
try,
    PlotLabels=evalin('base','PlotLabels');
catch
    PlotLabels=[];
end
try,
    Ap=evalin('base','Ap');
    go=1;
catch,
    msgbox('Cannot find CS matrix (variable Ap). Please load a workspace or run NCA');
    go=0;
    return;
end

if go==1,
    [filename,path]=uiputfile('*.txt','Enter Filename');
    if ~isequal(filename,0),
        file=strcat(path,filename);
        write_tfa_matrix(Ap,gene_ids,PlotLabels,file);
        msgbox(strcat('Successfully wrote A (CS values) Matrix to file ', file));
    end
end


% --------------------------------------------------------------------
function boxplot_tf_by_gene_Callback(hObject, eventdata, handles)
% hObject    handle to boxplot_tf_by_gene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try,
    PlotLabels=evalin('base','PlotLabels');
    evalin('base','plotSD=1;');    
    multi_plot
catch,
    warndlg('Please load a previously saved workspace, or start a new NCA analysis');
end


%
%




% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try,
    P=evalin('base','P');
    tfa_check;
catch,
    warndlg('Please load or run an NCA Analysis first');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try,
    P=evalin('base','P');
    tfa_check;
catch,
    warndlg('Please load or run an NCA Analysis first');
end



% --- Executes on button press in prematched_Checkbox.
function prematched_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to prematched_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of prematched_Checkbox




% --- Executes on button press in ztransform_data.
function ztransform_data_Callback(hObject, eventdata, handles)
% hObject    handle to ztransform_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ztransform_data

if get(hObject,'Value'),
    try,
        E=evalin('base','E');
        E=zscore(E')';
        evalin('base','Eorig=E;');
        assignin('base','E',E);
        msgbox('You have z-transformed the expression data.  You must rerun Step 3.');
    catch,
        warndlg('Cannot find Microarray Data. Please Reload');
        set(hObject,'Value',0);
        return;
    end

else,
    set(hObject,'Value',0);
    try,
        evalin('base','E=Eorig;');
    end
end


% --------------------------------------------------------------------
function export_tfa_devs_Callback(hObject, eventdata, handles)
% hObject    handle to export_tfa_devs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try,
    exp_labels=evalin('base','Experiment_labels');
catch,
    exp_labels=[];
end
try,
    PlotLabels=evalin('base','PlotLabels');
catch
    PlotLabels=[];
end
try,
    P_dev=evalin('base','P_dev');
    go=1;
catch,
    msgbox('Cannot find TFA Standard Deviation matrix (variable P_dev). Please load a workspace or run an NCA analysis.');
    go=0;
    return;
end

if go==1,
    [filename,path]=uiputfile('*.txt','Enter Filename');
    if ~isequal(filename,0),
        file=strcat(path,filename);
        write_tfa_matrix(P_dev,PlotLabels,exp_labels,file);
        msgbox(strcat('Successfully wrote TFA Matrix to file ', file));
    end
end




% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox5




% --------------------------------------------------------------------
function sigTestTFA_Callback(hObject, eventdata, handles)
% hObject    handle to sigTestTFA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% add code to check the statistical significance of each TFA.
% first check that NCA has been run.
% then compute random TFAs
% allow the user to choose the option they want... either 0-mean random
% variable -- or use NCA sample with bootstrap.
% also incorporate random periodicity test in menu.

iterations=evalin('base','iterations');
Af=evalin('base','A_rank_final');
P=evalin('base','P');
Pf=randn(size(P));
Ef=evalin('base','Ef');
b=inputdlg('How Many Sampling Iterations do you want to Run?','TFA Significance Testing');
try,
    b=str2double(b);
catch,
    b=20;
end
for x=1:b,
    Erand=Ef(randperm(size(Ef,1)),randperm(size(Ef,2)));
    [At(:,:,x),Pt(:,:,x)]=gnca_cc_v2(Erand,Af,Pf,true);
end
for i=1:size(P,1),
    for j=1:size(P,2),
        Px(i,j)=mean(Pt(i,j,:));
        Pd(i,j)=std(Pt(i,j,:));
    end
end

for i=1:size(P,1),
    for j=1:size(P,2),
        Z(i,j)=ztest(P(i,j),Px(i,j),Pd(i,j),0.05,'both');
    end
end


figure;
Peaks = axes('Position',[.1 .1 .75 .8],'Visible','off');
% Make the contour plot
pcolor(Z)
colorbar;
assignin('base','Significant_TFAs_by_Experiment',Z);
Plabels=evalin('base','Plabels');
set(Peaks,'YTickLabel',Plabels)
set(Peaks,'YMinorTick','on');
set(Peaks,'YTickLabelMode','manual');
set(Peaks,'YTickMode','manual');
set(Peaks,'YTick',[1:length(Plabels)]);
ylim([1 16]);
ylabel('TF Index');
xlabel('Experiment');
title('Heatmap of TFA Pertubation by Experiment (p-value < 0.05)'); 