


M=zeros(size(A0,2),(size(E0,2)-1));

K=zeros(18,size(E0,2));

for i=1:size(A0,2),
    for j=2:(size(E0,2)),
      h=hist(E0(find(abs(E0(find(A0(:,i)>0),j))>2),1),17);
      if length(h)>0,
        if max(h)>0,
          z=h(find(h>0));
          M(i,j)=z(1);
          for a=1:length(z),
              K(z(a),j)=K(z(a),j)+1;
              
          end
        end
      end
    end
end