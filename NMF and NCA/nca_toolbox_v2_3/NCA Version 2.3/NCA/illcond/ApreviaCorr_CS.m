[N,M]=size(E);
L=size(A,2);

[a,p]=gncar(E,A,P1);
[an,pn]=normalization(a,p);

Ptemp=[];    
for k=1:L
    [temp1,temp2]=max(abs(pn(k,:)));
    fact=1/temp1;
    Ptemp(k,:)=pn(k,:)*fact;
end

Enorm=[];
for k=1:N
    [temp1,temp2]=max(abs(E(k,:)));
    fact=1/temp1;
    Enorm(k,:)=E(k,:)*fact;
end

Acs=[];
Acor=[];
for k=1:N
    for j=1:L
        cs=Ptemp(j,:)'\Enorm(k,:)';
        Acs(k,j)=cs;
        cor=corrcoef(Ptemp(j,:)',Enorm(k,:)');
        Acor(k,j)=cor(1,2);
    end
end

Ascore=abs(Acs)+abs(Acor);
thres=prctile(Ascore(:),75);

Apre=[];
for k=1:N
    for j=1:L
        score=Ascore(k,j);
        if score<thres
            Apre(k,j)=0;
        else
            Apre(k,j)=1;
        end
    end
end
         
Arate=A-Apre;
fp=length(find(Arate==-1))/M/N;
fn=length(find(Arate==1))/M/N;



iter=10;
for k=1:10
    

s=sum(abs(sign(A)),2);
L=size(A,2);
comp=find(s>1);
ncomp=length(comp);
A1g_s=A;
A1g_s(comp,:)=zeros(ncomp,L);

A1g2_c=A;
[row,col]=find(A);
nnz=length(row);
for k=1:nnz
    cs=P1gn_n(col(k),:)'\Eg1(row(k),:)';
    if cs<0.3 & cs>-0.3
        A1g2_c(row(k),col(k))=0;
    end
end

A1g_new2=A1g_s+A1g2_c;
