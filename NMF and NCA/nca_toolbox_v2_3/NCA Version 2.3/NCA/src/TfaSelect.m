% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% TFASELECT - Automatic TFs selection for NCAweb
%
% Riccardo Boscolo - University of California - Los Angeles
% and Young L Yang - University of California - Los Angeles
%  27/01/2004 - All Rights Reserved

function [id,u]= TfaSelect(A0,Mm,tfa_ids,verbose)

% Check input arguments
if (nargin < 4 | isempty(verbose))
    verbose=0;
end

M = min(Mm,length(tfa_ids));
% Initialize variables
id=0;
max_iter=5;
max_attempt=1;
dec_step= round(.2*M);
min_size= 5;

% Select the subnetwork size
%K= round(0.8*M);

K = M;

attempt=0;
iter=max_iter+1;
while (~id)
    % Check if we are iterating with no success
    if (iter > max_iter)
        % Check whether we have been trying this unsuccessfully too many times  
        attempt=attempt+1;
        if (attempt > max_attempt)
            % Try to reduce the network size
            K=K-dec_step;
            if (K < min_size) break; end
            attempt=1;
        end
        
        % Select an random set of 'K' TRs
        rand('state',sum(100*clock));
        p= randperm(size(A0,2));
        u= sort(p(1:K));
        
        % Display current selection of TRs
        if (verbose)
            fprintf('\nTranscriptional subnetwork size: %i',K);
            fprintf('\nUser selection: ');
            for m=1:K, fprintf('%s ',char(tfa_ids(u(m)))); end
            fprintf('\n');
            fprintf('\nIter# | Added  | Removed |  Rank\n');
        end
        iter=1;
    end    

    % Run main routine
    [id,A,genes,add,remove]= RankEval(A0,u);
    
    if (~id)
        % Adding and removing the TRs (just one at the time for now)
        if (verbose)
            for k=1:1, % length(remove)
                fprintf(' %2i   |  %4s  |  %4s  |  %i',iter,char(tfa_ids(add(k))),char(tfa_ids(remove(k))),R);
            end
            fprintf('\n');
        end

        % Identify the TFs to remove
        i= find(u==remove(1));
        u= u([1:i-1 i+1:K]);
        % Add TFs 
        u= sort([u add(1)]);
    end
    iter=iter+1;
end
if (id & verbose)
    fprintf('\nThe following set of TRs and genes was automatically identified:\n');
    for m=1:K, fprintf('%s ',char(tfa_ids(u(m)))); end
    fprintf('\n');
end
if (~id & verbose)
    fprintf('\nNo identifiable transcriptional subnetwork was discovered\n');
end