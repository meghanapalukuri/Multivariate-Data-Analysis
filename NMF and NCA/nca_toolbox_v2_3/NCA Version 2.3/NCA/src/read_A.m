% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
% read A matrix from a text file

function [A,Alabels,Plabels] = read_A( filename )

% Plabels is cellstr
% Alabels is cellstr
% A is double connectivity matrix with 1 and 0s


% 10/13/2005 :: SJG :: Rewrote the read program so that it does not use
% the statistics toolbox.  DONE 10/13/2005

    eol='\n';
    delim='\t';
    %fid=fopen(filename);
    %fseek(fid,0,-1);
    c = textread(filename,'%s','delimiter',eol);
    %fclose(fid);
    charC = char(c{:});
    
    % read header
    % now read each line
    Plabels=strread(charC(1,:),'%s','delimiter','\t');
    Plabels=Plabels(2:length(Plabels));
    
    Alabels={};
    q=1;
    for w=2:size(charC,1),
       line = strread(charC(w,:),'%s','delimiter','\t');
       Alabels{q} = line{1};
       q=q+1;
       for x=2:length(line),
               A(w-1,x-1) = str2double(line{x});
       end
    end
    Alabels=Alabels';
    % now we need to split each line of c by delim
