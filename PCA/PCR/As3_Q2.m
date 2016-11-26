clc
clear all
close all

F1=load('flowdata1.mat')
F2=load('flowdata2.mat')

[N,n]=size(F1.Fmeas); % N - no. of samples, n- no. of variables
X=(F1.Fmeas)';
Atrue=F1.Atrue;
STD=F1.std; % Standard deviation of errors in variables
%% Data Set2 - Uncomment and Run to get solutions for 2nd dataset
% 
% X=(F2.Fmeas)';
% Atrue=F2.Atrue;

%% Q2a. 
% True Regression Matrix
Ad=[Atrue(:,3) Atrue(:,4) Atrue(:,5)];
Ai=[Atrue(:,1) Atrue(:,2)];

RegressionMatrix=-inv(Ad)*Ai

% Applying PCA to estimate Regression Matrix
k=2;
avg=mean(X,2);
Xs=X-repmat(avg,1,N);
[U S V]=svd(Xs,'econ');
Ahat=(U(:,k+1:n))';
Adhat=[Ahat(:,3) Ahat(:,4) Ahat(:,5)];
Aihat=[Ahat(:,1) Ahat(:,2) ];

RegressionMatrixEst=-inv(Adhat)*Aihat

RegError=RegressionMatrixEst-RegressionMatrix;
MaxAbsErrorPCA=max(max(abs(RegError)))

%% Q2b.

% figure(1)
% plot([1:5],diag(S))
% xlabel('Variable number')
% ylabel('Singular Value')


%% Q2c : PCA with Scaling

%Cov=diag(STD.*STD)
%L=chol(Cov,'lower')
Linv=inv(diag(STD));
Xs=Linv*Xs;

[U S V]=svd(Xs,'econ');
Ahat=(U(:,k+1:n))';
Adhat=[Ahat(:,3) Ahat(:,4) Ahat(:,5)];
Aihat=[Ahat(:,1) Ahat(:,2) ];

RegressionMatrixEst=-inv(Adhat)*Aihat

RegError=RegressionMatrixEst-RegressionMatrix;
MaxAbsErrorPCA=max(max(abs(RegError)))

%% Q2d : Performing PCA and Estimating different regression matrices with different choice of independent variables

% 
% X=(F2.Fmeas)';
% Atrue=F2.Atrue;

ind=[4,5] % Choice of independent variables
dep=[1,2,3 ];
% True Regression Matrix
Ad=[Atrue(:,dep(1)) Atrue(:,dep(2)) Atrue(:,dep(3))];
Ai=[Atrue(:,ind(1)) Atrue(:,ind(2))];

RegressionMatrix=-inv(Ad)*Ai

% Applying PCA to estimate Regression Matrix
k=2;
avg=mean(X,2);
Xs=X-repmat(avg,1,N);
[U S V]=svd(Xs,'econ');
Ahat=(U(:,k+1:n))';
Adhat=[Ahat(:,dep(1)) Ahat(:,dep(2)) Ahat(:,dep(3))];
Ad_determinant=det(Adhat)
Aihat=[Ahat(:,ind(1)) Ahat(:,ind(2)) ];

Ad=[Atrue(:,dep(1)) Atrue(:,dep(2)) Atrue(:,dep(3))];
Ai=[Atrue(:,ind(1)) Atrue(:,ind(2))];

RegressionMatrixEst=-inv(Adhat)*Aihat

RegError=RegressionMatrixEst-RegressionMatrix;
MaxAbsErrorPCA=max(max(abs(RegError)))
