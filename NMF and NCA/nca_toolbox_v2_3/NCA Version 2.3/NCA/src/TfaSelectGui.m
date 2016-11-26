% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% TfaSelectGUI.m
% (c) 2004 Simon J Galbraith
% Summary:
%     Creates GUI Wizard to modify GUI to work with Yeast and Ecoli Connectivity and Expression Data


function TfaSelectGUI(A0,gene_ids,tfas)

fig = open('ncaselect.fig');
     handles=guihandles(fig);
     tfas = cellstr(tfas);
     data=tfas;
assignin('base','h',fig);  
%global gtfa_ids_final;

% need to pass the A matrix and tfa_ids to the GUI.
handles.A0 = A0;
handles.tfas = tfas;
handles.tfa_ids=tfas;
handles.gene_ids = gene_ids;
guidata(handles.listbox1,handles);
set(handles.listbox1,'String',data);
set(handles.listbox2,'String','');
set(handles.pushbutton1,'String','Add TF'); 
set(handles.pushbutton2,'String','Delete TF');
set(handles.pushbutton3,'String','Auto Fill TF list');
set(handles.popupmenu1,'String',{[1:length(tfas)]});
