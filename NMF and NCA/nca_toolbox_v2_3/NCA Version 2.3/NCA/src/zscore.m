% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
% Written by: Simon J Galbraith
function Z=zscore(M),
   Z=zeros(size(M));
   for j=1:size(M,2),
       u_j=mean(M(:,j));
       o_j=std(M(:,j));
       Z(:,j)=(M(:,j)-u_j)/o_j;
   end
end