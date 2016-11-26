% MATCHPRUNE - Match gene expression level matrix to connectivity 
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu)
%
% This code was developed by Riccardo Boscolo - University of California - Los Angeles
%  01/06/2004 - All Rights Reserved and published with the Biospice release
%  of NCA.
%
% Modified:  Simon J Galbraith

function [E0,A0,gene_ids,tfa_ids]= MatchPrune(Em,gene_names_E,Am,gene_names_A,tfa_names,organism)

% Check input arguments
if (nargin < 4)
    error('Usage: E0,A0,gene_ids,tfa_ids]= MatchPrune(Em,gene_names_E,Am,gene_names_A,tfa_names,organism)');
end

[Ne,M]= size(Em);
[Na,L]= size(Am);

if strcmpi(organism,'S. cerevisiae'),
  load('genebanks/yeast_genebank.mat');
elseif strcmpi(organism,'E. coli'),
  load('genebanks/ecoli_genebank.mat');
elseif strcmpi(organism,'H. sapiens'),
   load('genebanks/human_genebank.mat');  
end
   
if (length(gene_names_E)~=Ne)
    error('The vector ''gene_names_E'' should have as many elements as the # of rows in Em');
end

if (size(gene_names_A,1)~=Na)
    error('The vector ''gene_names_A'' should have as many elements as the # of rows in Am');
end

% Initialize variables
gene_ids= cell(min(Ne,Na),1);
E0= zeros(min(Ne,Na),M);
A0= zeros(min(Ne,Na),L);
gene_matches= zeros(10,1);
if ischar(gene_names_E),
  gene_names_E= cellstr(gene_names_E);
end

if ischar(gene_names_A),
    gene_names_A= cellstr(gene_names_A);
end
N=0;
% Find common genes
for n=1:Ne,
    gene_match_ind=0;
    % Find synonyms of gene_names_E(n)
    syn_ind= find(sum(strcmpi(gene_names_E(n),genebank),2));
if (length(syn_ind) > 1)
    syn_ind=syn_ind(1);  % choose the first one
    %       printf(['Found multiple matches for gene ' char(gene_names_E(n)) ' in the gene database.']);
end
    
    % Attempt to match each synonim to the name in the connectivity file
    for k=1:size(genebank,2),
        %genebank(syn_ind,k)       
        if (length(char(genebank(syn_ind,k))) > 0)
            v= find(sum(strcmpi(genebank(syn_ind,k),gene_names_A),2));
	    %warning([length(v)])
            if (length(v) > 1)
                warning(['Found multiple matches for gene ' char(genebank(syn_ind,k)) 'in the connectivity matrix.']);
            end
            if (length(v) > 0)
                gene_match_ind= gene_match_ind+1;
                gene_matches(gene_match_ind)= v(1);
            end
        end
    end
    
    % If the gene is matched, select the corresponding regulatory pattern 
    if (gene_match_ind >= 1)
        %fprintf('Found a Match\n');
        N=N+1;
        E0(N,:)= Em(n,:);
        A0(N,:)= Am(gene_matches(1),:);
        gene_ids(N)= gene_names_E(n);
    end
    % Check what happened if multiple matches occured
    if (gene_match_ind > 1)
        v= find(gene_matches==gene_matches(1));
        if (length(v) < gene_match_ind)
            warning(['Gene ' char(gene_names_E(n)) ' resulted in multiple matches, the first one was selected.']);
        end
    end    
    % No matches, too bad...
%     if (length(gene_matches)==0)
%         fprintf('\nCould not find gene %s in the connectivity master file.',char(gene_names_E(n)));
%     end
end

% Eliminate disconnected TFs
v= find(sum(abs(A0)));

A0= A0(:,v);
tfa_ids= tfa_names(v);

% remove any TFS not present in the Microarray and any genes only connected
% to it.


% Adjust vector sizes
E0= E0(1:N,:);
A0= A0(1:N,:);
gene_ids= gene_ids(1:N);