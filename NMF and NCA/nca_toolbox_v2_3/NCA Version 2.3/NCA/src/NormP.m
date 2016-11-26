function [norA,norP]=NormP(A,P)
%(c) 2005. University of California, Los Angeles. All Rights Reserved.
% NormP normalizes each k TF by the absolute highest norm of P(k,:) 

    norA=A;
    norP=P;
    for k=1:size(P,1),
        norP(k,:)=P(k,:)/abs(norm(P(k,:),2));
        norA(:,k)=A(:,k)*abs(norm(P(k,:),2));
    end
end