clc
clear all
close all

%% a)PCA
Data=load('ncadata.mat');

ModelData=Data.measabs;
[nsamples,nvar]=size(ModelData);
PureComp=Data.pureabs;
%Denoising the data using PCA
k=3;
[Umatrix,DenoisedModelData] = PCApar(ModelData,k);
Ahat=Umatrix(:,1:k)';
% Rotation Matrix based on connectivity 7*3 matrix in the order of Co, Cr,
% Ni
ConnectivityMatrix=[ 1 1 0 ; 1 0 1; 0 1 1; 1 0 1; 1 1 0; 1 0 1; 0 1 1];
RotMat=ConnectivityMatrix
PureCompPred=(DenoisedModelData*RotMat)';


% Calculating RMSE
RMSE=0;
for i=1:3
    for j=1:176
RMSE=RMSE+(PureComp(i,j)-PureCompPred(i,j))^2;
    end
end
RMSE_PCA=sqrt(RMSE)

% Correlation Coefficients

for i = 1:3
    for j = 1:3
        correlation_coeff(i,j) = corr2(PureComp(i,:),PureCompPred(j,:));
    end
end

correlation_coeff_PCA=correlation_coeff

% Comp1=corrcoef(PureComp(1,:),PureCompPred(1,:))
% Comp2=corrcoef(PureComp(2,:),PureCompPred(2,:))
% Comp3=corrcoef(PureComp(3,:),PureCompPred(3,:))

% figure(1)
% hold on
% plot(300:2:650,PureComp(1,:),'b-')
% plot(300:2:650,PureCompPred(1,:),'g-')
% hold off
% 
% figure(2)
% hold on
% plot(300:2:650,PureComp(2,:),'b-')
% plot(300:2:650,PureCompPred(2,:),'g-')
% hold off
% 
% figure(3)
% hold on
% plot(300:2:650,PureComp(3,:),'b-')
% plot(300:2:650,PureCompPred(3,:),'g-')
% hold off

figure(1)
plot(300:2:650,PureComp(1,:),'b-')
title('Pure Co Spectrum')

figure(2)
plot(300:2:650,PureCompPred(1,:),'g-')
title('Predicted Pure Co Spectrum')

figure(3)
plot(300:2:650,PureComp(2,:),'b-')
title('Pure Cr spectrum')

figure(4)
plot(300:2:650,PureCompPred(2,:),'g-')
title('Predicted Pure Cr spectrum')

figure(5)
plot(300:2:650,PureComp(3,:),'b-')
title('Pure Ni spectrum')

figure(6)
plot(300:2:650,PureCompPred(3,:),'g-')
title('Predicted Pure Ni spectrum')
%% b) NCA
% Formatting and saving for NCA solver
Name=[1 , 2, 3];
Number=0:7;
ConMat=cat(1,Name,ConnectivityMatrix);
ConMat=cat(2,Number',ConMat);
save('Connectivity.txt','ConMat','-ascii','-tabs');

Name=1:176;
Number=0:7;
ModelTab=cat(1,Name,ModelData);
ModelTab=cat(2,Number',ModelTab);
save('Model.txt','ModelTab','-ascii','-tabs');

% Applying NCA; running NCA and saving the A and E matrices
Amatrix=load('A.mat');
Ematrix=load('E.mat');

A=Amatrix.A0;
E=Ematrix.E0;
P=pinv(A)*E;

PureCompPred=P;
% Calculating RMSE
RMSE=0;
for i=1:3
    for j=1:176
RMSE=RMSE+(PureComp(i,j)-PureCompPred(i,j))^2;
    end
end
RMSE_NCA=sqrt(RMSE)

% Correlation Coefficients

for i = 1:3
    for j = 1:3
        correlation_coeff(i,j) = corr2(PureComp(i,:),PureCompPred(j,:));
    end
end

correlation_coeff_NCA=correlation_coeff
% Comp1=corrcoef(PureComp(1,:),PureCompPred(1,:))
% Comp2=corrcoef(PureComp(2,:),PureCompPred(2,:))
% Comp3=corrcoef(PureComp(3,:),PureCompPred(3,:))

%  figure(1)
% hold on
% plot(300:2:650,PureComp(1,:),'b-')
% plot(300:2:650,PureCompPred(1,:),'g-')
% hold off
% 
% figure(2)
% hold on
% plot(300:2:650,PureComp(2,:),'b-')
% plot(300:2:650,PureCompPred(2,:),'g-')
% hold off
% 
% figure(3)
% hold on
% plot(300:2:650,PureComp(3,:),'b-')
% plot(300:2:650,PureCompPred(3,:),'g-')
% hold off
% % 
figure(1)
plot(300:2:650,PureComp(1,:),'b-')
title('Pure Co Spectrum')

figure(2)
plot(300:2:650,PureCompPred(1,:),'g-')
title('Predicted Pure Co Spectrum')

figure(3)
plot(300:2:650,PureComp(2,:),'b-')
title('Pure Cr spectrum')

figure(4)
plot(300:2:650,PureCompPred(2,:),'g-')
title('Predicted Pure Cr spectrum')

figure(5)
plot(300:2:650,PureComp(3,:),'b-')
title('Pure Ni spectrum')

figure(6)
plot(300:2:650,PureCompPred(3,:),'g-')
title('Predicted Pure Ni spectrum')