

% AICfit :: Simon J Galbraith
%
% randomly perturb an edge (rewire, add, delete) 
%

function [AS,A0] = AICfit(E0,A0),
AS=[];
Aorig=A0;
AICscore_last=inf;


for t=1:500,
    T=501-t;
    while 1,
        z=randperm(size(E0,1));
        z=z(1);
        if length(find(A0(z,:)>0))>1,
            x=randn();
            if x>0.5,  % find an edge to remove
                g=find(A0(z,:))>0;
                y=randperm(length(g));
                A0(z,y(1))=0;
                u=[1:size(A0,2)];
                if RankEval(A0,u),
                    break;
                else,
                    A0(z,y(1))=1;
                end
            end
        end
    end

    [Az,Pz]=gncar_cc(E0,A0);
    R=norm(E0-Az*Pz);
    AICscore = -2*log(R)+2*length(find(A0>0));
    x=randn();

    if AICscore < AICscore_last,
        AICscore_last = AICscore;
        Aorig=A0;
    elseif x > exp(abs(AICscore_last-AICscore)/T),
        AICscore_last = AICscore;
        Aorig=A0;
    else,
        A0=Aorig;
    end
    AS = [AS,AICscore_last];
end