clc
clear all
close all

Data=load('Inorfull.mat');
% Data structure consists of the follow attributes:
%           CONC: [130x3 double]
%           DATA: [130x176 double]
%         PureCo: [1x176 double]
%         PureCr: [1x176 double]
%         PureNi: [1x176 double]
%            WAV: [1x176 double]
%        stdDATA: [130x176 double]
%     PureCoCONC: 0.1720
%     PureCrCONC: 0.0764
%     PureNiCONC: 0.1965 
    
PureComp=[Data.PureCo; Data.PureCr; Data.PureNi];
% Final Data: Absorbance spectra (176 pts.) measured for each sample
% (Total: 26 samples). , each consisting of 3 components with their
% concentrations noted. 

% Model Data : 25 samples. Validation data: Left out sample. Calculate RMSE
% in predicting this conc.

% Goal is to predict concentrations given absorbance of mixture

%% Picking one replicate measurement per sample- Random measurement
NewData=[];
NewStdData=[];
NewConc=[];

rand('seed',95) % So as to provide consistent measurements
for i=1:26
    dt_num = randi(5);
    istart = 5*(i-1)+dt_num;
    NewStdData = [NewStdData; Data.stdDATA(istart,:)];
     NewData = [NewData; Data.DATA(istart,:)];
    NewConc = [NewConc;Data.CONC(istart,:)];
end

X=NewData;
Y=NewConc;

 
%% Q1.e Average measurement - Uncomment this to solve for average case

NewData=ones(26,176);
NewStdData=ones(26,176);
NewConc=ones(26,3);

count=1;
for i=1:5:130
      NewConc(count,:)=(Data.CONC(i,:)+Data.CONC(i+1,:)+Data.CONC(i+2,:)+Data.CONC(i+3,:)+Data.CONC(i+4,:))/5;
    NewData(count,:)=(Data.DATA(i,:)+Data.DATA(i+1,:)+Data.DATA(i+2,:)+Data.DATA(i+3,:)+Data.DATA(i+4,:))/5;
     NewStdData(count,:)=(Data.stdDATA(i,:)+Data.stdDATA(i+1,:)+Data.stdDATA(i+2,:)+Data.stdDATA(i+3,:)+Data.stdDATA(i+4,:))/5;
    count=count+1;
end

X=NewData;
Y=NewConc;
%% Q1. f: With end wavelengths removed- Uncomment and get results for Q. b and c

 X=X(:,26:151);
%% Q1a. OLS 
NewData = X-(ones(26,1)*mean(X));
NewConc = Y - (ones(26,1)*mean(Y));

[nsamples , nvar] = size(NewData);
sumsqrerr = zeros(3,1);
RMSE=0;
for k=1:26 
 ModelConc =[NewConc(1:k-1,:);NewConc(k+1:nsamples,:)];
    ModelData=[NewData(1:k-1,:);NewData(k+1:nsamples,:)];
     ModelStdData=[NewData(1:k-1,:);NewData(k+1:nsamples,:)];
     
     ValidateConc=NewConc(k,:);
    ValidateData=NewData(k,:);
     ValidateStdData=NewStdData(k,:);

 ModelX=ModelData;
 ModelY=ModelConc;
 %OLS - assuming errors in absorbance(ModelX)- Inverse OLS
 alphainv=pinv(ModelY'*ModelX)*ModelY'*ModelY;
 
PredictConc=ValidateData*alphainv;
err=0;
for i = 1:3
error = ValidateConc(i)-PredictConc(i);
sumsqrerr(i) = sumsqrerr(i) + error*error;
end

end

RMSE_1a = sqrt(sumsqrerr/nsamples);
RMSETotal_1a = sqrt(sum(sumsqrerr)/(nsamples*3));


%% Q1.d: Q1b and 1c with different no. of PCs - Change here from 1 to 6
nPC=6;

%% Q1.b: PCR with error free scores
NewData = X-(ones(26,1)*mean(X));
NewConc = Y - (ones(26,1)*mean(Y));
%Denoising the data using PCA
[DenoisedModelData,ScoresMatrix,SingularValues] = PCA(NewData',nPC);

NewData=ScoresMatrix;

[nsamples , nvar] = size(NewData);
sumsqrerr = zeros(3,1);
RMSE=0;
for k=1:26 
 ModelConc =[NewConc(1:k-1,:);NewConc(k+1:nsamples,:)];
    ModelData=[NewData(1:k-1,:);NewData(k+1:nsamples,:)];
     ModelStdData=[NewData(1:k-1,:);NewData(k+1:nsamples,:)];
     
     ValidateConc=NewConc(k,:);
    ValidateData=NewData(k,:);
     ValidateStdData=NewStdData(k,:);

 ModelX=ModelData;
 ModelY=ModelConc;
 %OLS - assuming errors in concentrations(ModelY)-  OLS
 alphainv=pinv(ModelX'*ModelX)*(ModelY'*ModelX)';
 
PredictConc=ValidateData*alphainv;
err=0;
for i = 1:3
error = ValidateConc(i)-PredictConc(i);
sumsqrerr(i) = sumsqrerr(i) + error*error;
end

end

RMSE_1b = sqrt(sumsqrerr/nsamples)
RMSETotal_1b = sqrt(sum(sumsqrerr)/(nsamples*3))

%% Q1.c: PCR with error free concentrations
NewData = X-(ones(26,1)*mean(X));
NewConc = Y - (ones(26,1)*mean(Y));
%Denoising the data using PCA
[DenoisedModelData,ScoresMatrix,SingularValues] = PCA(NewData',nPC);

NewData=ScoresMatrix;

[nsamples , nvar] = size(NewData);
sumsqrerr = zeros(3,1);
RMSE=0;
for k=1:26 
 ModelConc =[NewConc(1:k-1,:);NewConc(k+1:nsamples,:)];
    ModelData=[NewData(1:k-1,:);NewData(k+1:nsamples,:)];
     ModelStdData=[NewData(1:k-1,:);NewData(k+1:nsamples,:)];
     
     ValidateConc=NewConc(k,:);
    ValidateData=NewData(k,:);
     ValidateStdData=NewStdData(k,:);

 ModelX=ModelData;
 ModelY=ModelConc;
 %OLS - assuming errors in scores(ModelX) i.e absorbance- Inverse OLS
 alphainv=pinv(ModelY'*ModelX)*ModelY'*ModelY;
 
PredictConc=ValidateData*alphainv;
err=0;
for i = 1:3
error = ValidateConc(i)-PredictConc(i);
sumsqrerr(i) = sumsqrerr(i) + error*error;
end

end

RMSE_1c = sqrt(sumsqrerr/nsamples)
RMSETotal_1c = sqrt(sum(sumsqrerr)/(nsamples*3))