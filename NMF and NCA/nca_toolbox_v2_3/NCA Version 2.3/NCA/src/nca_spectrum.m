function [A,P]=nca_spectrum(E,initialA,wavelength_index,absorption);
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
%This version is used to decompose spectra (i.e. E matrix) with random initial guesses
%for concentration matrix (i.e. initialA matrix).  The final solution A
%and P will have the minimum error (i.e. min ||E-AP||)
%Input:
%  E: is the spectra of all mixures, where E(ij) is the absorption of mixture i at wavelength j     
%  initial A: the concentration pattern of all mixture, where A(ik) represents the existence
%             of component k exists in mixture i. 
%  wavelength: is actually the index j in E matrix, which corresponds to the specific wavelength that we know
%              the absorption of pure components
%  absorption: absorption of pure components at the specified wavelength j
%Output:
%   A
%   P: pure spectra, where P(kj) is the absorption of component k at
%   wavelength j
%--------------------------------------------------------------------------
%Linh Tran, UCLA.
%Last revision August 20, 2005
%-------------------------------------------------------------------------

%Check matrix size
L=size(initialA,2);
[N,M]= size(E);
if (size(initialA,1) ~= N)
    error('The matrices E and initialA must have the same number of rows.');
end

%Check NCA crieria
if (RankCheck(initialA) < L*(L-1))
    error('The system specified by ''initilaA'' is NOT identifiable. Abort.');
end

%Run NCA
A=[];
P=[];
randA=rand(N,L).*sign(initialA);
[A,P]=nca(E,randA);

%Normalize solution
for j=1:L
    fraction=absorption(j)/P(j,wavelength_index(j));
    P(j,:)=P(j,:)*fraction;
    A(:,j)=A(:,j)/fraction;
end

%--------------------------------------------------------------------------

% RANKCHECK - Checks the rank of the augmented system matrix after pruning for
%       non-zero entries in A.
% 
%  R = RankCheck(A) - Returns the rank of Au= diag{A,A,...,A} after pruning the
%       rows corresponding to the non-zero entries in A.

function R= RankCheck(alpha)

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

%-------------------------------------------------------------
%new version of NCA program that check NCA criteria during iteration
function [A,P]= nca(E,A0)

% Variable initialization
[N,M]= size(E);
L= size(A0,2);
A= A0;
P= zeros(L,M);

%Check NCA crieria
if (RankCheck(A0) < L*(L-1))
    error('The system specified by ''initialA'' is NOT identifiable. Abort.');
end

% Parameters and constraints for the EM algorithm
epsl= 1.0e-9;           % --> convergence threshold

% Identify the sparsity pattern of 'alpha'
sp_pattern= sign(A0);

% EM algorithm - main loop
delta_R_rel= 1.0;
R= 1e6; 
iter=0;
while (delta_R_rel > epsl)
    iter=iter+1;
    % Step 1 - Re-estimation of 'P'
    if (RankCheck(A)<L*(L-1))   %Check the 2nd criterion
        A=rand(N,L).*sign(A0);   %re-start iteration
    end
    P= A\E;
    
    % Step 2 - Re-estimation of 'A' solving a sequence of least squares problems
    f=0;  %Pseudo variable to check the third criterion
    n=1;
    while (f==0 & n<=N)
        % Build the low-dimension least square problem for the n_th gene
        index= find(sp_pattern(n,:) ~= 0);
        P_red= P(index,:);  
        if (cond(P_red)>100)  %Check the 3rd criterion, if it is violated,
            f=1;
            A=rand(N,L).*sign(A0);  %re-start iteration
        else
            A(n,index)= E(n,:)/P_red;
            n=n+1;
        end
    end
    
    % (c) - Evaluate the current mismatch
    err= E-A*P;
    R_old= R;
    R= sum(err(:).^2);
    delta_R_rel= abs(R_old-R)/R_old;
    if 0==mod(iter,1000),
        fprintf(1,'ssr= %.8f\n',R);
    end
end