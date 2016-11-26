
function [gene_set] = diffexp(E,Elabels,fold)

% (c) 2004 Simon J Galbraith.
% select genes from E that are only differentially expressed in one or more
% columns.

Eabs = abs(E);
gene_set=[];
for i=1:size(E,2),
  gene_set=union(gene_set,find(Eabs(:,i)>=2.0));
end

