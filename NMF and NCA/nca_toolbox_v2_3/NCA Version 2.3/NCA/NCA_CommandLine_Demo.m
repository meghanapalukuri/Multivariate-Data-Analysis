% Script:   NCA_CommandLine_Demo.m 
% Author:   Simon J Galbraith.
% Description:   This script runs the NCA toolbox from the command line so that the 
% user can see how the NCA toolchain works.

%%%%%% Let's get started:
addpath('src');  % make sure we can see the NCA program files in the 'src' directory.
echo on;
organism='E. coli';

    % load in Connectivity Network (tab-delimited format)
[A,Alabels,Plabels]=read_A('Example Data/Ecoli_Connectivity_for_Kao_et_al_PNAS.txt');

    % load in Microarray dta (tab-delimited format)
[E,Elabels,Experiment_labels]=read_A('Example Data/Ecoli_Glucose_to_Acetate_Shift_Kao_et_al_PNAS.txt');

    % Match the data matrices and prune any disconnected TFs
[E0,A0,gene_ids,tfa_ids]=MatchPrune(E,Elabels,A,Alabels,Plabels,organism);

    % auto-select 5 random TFs that form an NCA compliant network
[id,u] = TfaSelect(A0,5,tfa_ids);  

    % Check the NCA compliance conditions
[id,A_rank_final,gene_id_final,add,remove]=RankEval(A0,u);

    % Now run NCA
if id,
    [Ap,P]=gnca_cc(E0(gene_id_final,:),A_rank_final);
    fprintf('Paused after running NCA analysis.  Now we will load the plotting function of the toolbox to visualize the results.\n');
    pause;
        % Normalize A and P
    [Ap,P]=Normalization(Ap,P);
    PlotLabels=tfa_ids(u);
    multi_plot;
else,
    fprintf('Cannot find an NCA compliant network\n');
end