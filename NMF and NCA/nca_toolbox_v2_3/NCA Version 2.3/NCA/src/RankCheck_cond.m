% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
% RANKCHECK - Checks the rank of the augmented system matrix after pruning for
%       non-zero entries in A.
% 
%  R = RankCheck(A) - Returns the rank of Au= diag{A,A,...,A} after pruning the
%       rows corresponding to the non-zero entries in A.

function R= RankCheck_cond(alpha)

[N,L]= size(alpha);

R= 0;   % --> initialize rank of the augmented sparse system
nz= 0;  % --> initialize total # of zeros in the augmented system

for l= 1:L,
    v= find(alpha(:,l)==0);
    nz= nz+length(v);
    if l==1
        reduceA=alpha(v,2:L);
    else
        reduceA=alpha(v,[1:l-1 l+1:L]);
    end
    c=cond(reduceA);
    if c>1E7
        r=L-2;
    else
        r=L-1;
    end
    R = R + r;
end