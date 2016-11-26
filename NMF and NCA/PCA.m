function DenoisedModelData = PCA(ModelData,k)
[M,N] = size(ModelData);
% subtract off the mean for each dimension
mn = mean(ModelData,2);
ModelData = ModelData - repmat(mn,1,N);
% construct the matrix Y
Y = ModelData' / sqrt(N-1);
% SVD does it all
[u,S,PC] = svd(Y,'econ');
% % calculate the variances
% S = diag(S);
% V = S .* S;
% % project the original data
% signals = PC' * data;
DenoisedModelData=u(:,1:k)*S(1:k,1:k)*PC(:,1:k)';
end