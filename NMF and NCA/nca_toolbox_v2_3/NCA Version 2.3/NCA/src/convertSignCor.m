function [newA,newP]=convertSignCor(orgP,A,P)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
%[newA,newP]=convertSignCor(orgP,A,P)
%Change sign of P base on reference matrix, orgP.
%If the correlation coeficient between TFA k in P and corresponding TFA in
%orgP is less than 0, change the sign of that TFA.
%Linh Tran, UCLA, 4/26/2005

if (nargin <3)
    error('You must input reference P along with A,P set');
end

L=size(orgP,1);

newA=A;
newP=P;

for k=1:L
    cc=corrcoef(orgP(k,:)',P(k,:)');
    if cc(2,1)<0
        newA(:,k)=-newA(:,k);
        newP(k,:)=-newP(k,:);
    end
end

function [newA,newP]= NormalizationP(oldA,oldP)

L=size(oldP,1);

for k=1:L
    [max_v,max_i]=max(abs(oldP(k,:)));
    newP(k,:)=oldP(k,:)/max_v;
    newA(:,k)=oldA(:,k)*max_v;
end