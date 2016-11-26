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
    
%%  a) Random replicate
NewData=[];
NewStdData=[];
NewConc=[];

rand('seed',95)
% count=1;
% for i=1:5:130
%     k=randi([i,i+4]);
%       NewConc(count,:)=Data.CONC(k,:);
%     NewData(count,:)=Data.DATA(k,:);
%      NewStdData(count,:)=Data.stdDATA(k,:);
%     count=count+1;
% end
for i=1:26
    dt_num = randi(5);
    istart = 5*(i-1)+dt_num;
    NewStdData = [NewStdData; Data.stdDATA(istart,:)];
     NewData = [NewData; Data.DATA(istart,:)];
    NewConc = [NewConc;Data.CONC(istart,:)];
end

% Final Data: Absorbance spectra (176 pts.) measured for each sample
% (Total: 26 samples). , each consisting of 3 components with their
% concentrations noted. 

PureComp=[Data.PureCo; Data.PureCr; Data.PureNi];

% Goal is to predict concentrations given absorbance of mixture

% Model: Absorbance spectra of a mixture = Linear Combination of absorbance
% spectras of individual pure component spectra with weights corresponding
% to concentrations of the components . Fit a model between mixture spectra and the weighted pure component spectra. Ideally coefficients should be 1 1
% 1 i.e , ideally : ModelData=m*ModelConc*PureComp
%  ModelY= NewConc;
%  ModelX=NewData;
 % ModelX=NewData*pinv(PureComp); % More accurate than with just ModelData, i.e here we are considering Beer Lambert's Laws for a better fit
 
 ModelY= NewConc;
%  Linv = inv(diag(mean(NewStdData)));
%  ModelX=NewData*Linv;
ModelX=NewData;


%Denoising the data using PCA
DenoisedModelData = PCA(ModelX,3);
DenoisedModelData =DenoisedModelData';
%DenoisedModelData = (DenoisedModelData')*(diag(mean(NewStdData)));
k=3;
% Applying NMF to extract the 3 sources
[W,H,iter,HIS]=nmf(DenoisedModelData,k);

% Finding correlation coefficients to determine which component it
% corresponds to


correlation_coeff = zeros(3,3);
for i = 1:3
    for j = 1:3
        correlation_coeff(i,j) = corr2(PureComp(i,:),H(j,:));
    end
end

correlation_coeff_rand=correlation_coeff
PureCompPred=[H(3,:) ; H(2,:); H(1,:)];
% count=1;
% it=size(W);
% Szer=sum(W,1);
% for i=1:it(2)
%     if(Szer(i)~=0)
%     Wnew(:,count)=W(:,i);
%     count=count+1;
%     end
% end
% 
%  
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

% it=size(Wnew);
% for i=1:it(2)
%     figure(i+3)
%     plot(300:2:650,Wnew(:,i),'g-')
% end

% b)With Average of the replicates
clear all
Data=load('Inorfull.mat');
% Average measurement
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

% Final Data: Absorbance spectra (176 pts.) measured for each sample
% (Total: 26 samples). , each consisting of 3 components with their
% concentrations noted. 

PureComp=[Data.PureCo; Data.PureCr; Data.PureNi];

% Goal is to predict concentrations given absorbance of mixture

% Model: Absorbance spectra of a mixture = Linear Combination of absorbance
% spectras of individual pure component spectra with weights corresponding
% to concentrations of the components . Fit a model between mixture spectra and the weighted pure component spectra. Ideally coefficients should be 1 1
% 1 i.e , ideally : ModelData=m*ModelConc*PureComp
%  ModelY= NewConc;
%  ModelX=NewData;
 % ModelX=NewData*pinv(PureComp); % More accurate than with just ModelData, i.e here we are considering Beer Lambert's Laws for a better fit
 
 ModelY= NewConc;
%  Linv = inv(diag(mean(NewStdData)));
%  ModelX=NewData*Linv;
ModelX=NewData;


%Denoising the data using PCA
DenoisedModelData = PCA(ModelX,3);
DenoisedModelData=DenoisedModelData';
% DenoisedModelData = (DenoisedModelData')*(diag(mean(NewStdData)));
k=3;
% Applying NMF to extract the 3 sources
[W,H,iter,HIS]=nmf(DenoisedModelData,k);

% Finding correlation coefficients to determine which component it
% corresponds to


correlation_coeff = zeros(3,3);
for i = 1:3
    for j = 1:3
        correlation_coeff(i,j) = corr2(PureComp(i,:),H(j,:));
    end
end

correlation_coeff_avg=correlation_coeff
PureCompPred=[H(3,:) ; H(2,:); H(1,:)];
% count=1;
% it=size(W);
% Szer=sum(W,1);
% for i=1:it(2)
%     if(Szer(i)~=0)
%     Wnew(:,count)=W(:,i);
%     count=count+1;
%     end
% end
% 
 
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

% it=size(Wnew);
% for i=1:it(2)
%     figure(i+3)
%     plot(300:2:650,Wnew(:,i),'g-')
% end

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