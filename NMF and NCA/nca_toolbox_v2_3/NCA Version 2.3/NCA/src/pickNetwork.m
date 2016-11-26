function [newS,A]=pickNetwork(oldS,selectedTF)
% (c) 2005 University of California Los Angeles, All Rights Reserved.
% The NCA Toolbox was written by Simon J Galbraith (sgalbrai@cs.ucla.edu)
% and Linh Tran (ltran@seas.ucla.edu).
%
%selectedTF vector stores the indices of interested TFs which are a subset
%of TFs stored in oldS


if nargin<2
    error('You must select interested TFs');
end

%initial network
Ao=oldS.data;
Eo=oldS.dataE;
geneName=oldS.genename;
tfName=oldS.tfname;
newS.exptname=oldS.exptname;

[N,L]=size(Ao);
Ls=length(selectedTF);

%triming A base on selectedTF
A=Ao;
startTF=1;
selectedTF=sort(selectedTF);
for k=1:L
    if startTF<=Ls
        id=selectedTF(startTF);
    else
        id=0;
    end
    if k~=id
        reggene=find(A(:,k)~=0);
        A(reggene,:)=zeros(length(reggene),L);
        startTF=startTF-1;
    end
    startTF=startTF+1;
end

ngenes=sum(abs(sign(A)),2);
selectedGenes=find(ngenes);
newA=A(selectedGenes,selectedTF);
newS.data=newA;
if ~isempty(Eo);
    temp=Eo(selectedGenes,:);
    newS.dataE=temp;
    %newS.dataE=Eo(selectedGenes,:);
else
    newS.dataE=[];
end
newS.genename=geneName(selectedGenes);
newS.tfname=tfName(selectedTF);

