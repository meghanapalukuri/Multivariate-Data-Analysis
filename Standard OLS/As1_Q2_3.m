clc
clear all
close all
%% Question 2

%Reading numerical data into num and column names into txt
[num,txt,raw] = xlsread('prostate.xlsx');
Nvar=length(txt); % No. of variables

%Determining clinical measure has the highest correlation with LPSA 
  Correlation= corrcoef(num); % Returns correlation coefficient matrix 
  
%Correlation coefficient b/w each clinical measure(X)and LPSA(Y)-lastcolumn)
  Corr_XY = Correlation(1:Nvar-1,Nvar)
   
% Determining measure with highest correlation to LPSA
 [ HighestCorr, Index]=max(Corr_XY) 
 ClinMeas=cell2mat(txt(Index)) % Clinical Measure with Highest Correlation to LPSA
 
 % Fitting linear model using standard ols (without weights)
 mdl = fitlm(num(:,Index),num(:,Nvar),'linear','VarNames',{ClinMeas,'LPSA'})
 
% Plotting the model
figure(1)
plot(mdl)
 CoeffOfDetermination = mdl.Rsquared.Ordinary % R squared Value
 %% Question 3
 
 %Reading numerical data into num2 
num2 = xlsread('anscombe.xls');

data1=num2(1:11,1:2); % [X Y] for 1st data set

 % Fitting linear model using standard ols (without weights)
 mdl1 = fitlm(data1(:,1),data1(:,2),'linear','VarNames',{'X1','Y1'})
 Rsq1 = mdl1.Rsquared.Ordinary % Goodness of fit
 figure(2) % Plotting the model
 plotResiduals(mdl1,'fitted')
 title('Dataset1')
 
data2=num2(1:11,3:4); % [X Y] for 2nd data set
 mdl2 = fitlm(data2(:,1),data2(:,2),'linear','VarNames',{'X2','Y2'})
 Rsq2 = mdl2.Rsquared.Ordinary % Goodness of fit
  figure(3) % Plotting the model
  plotResiduals(mdl2,'fitted')
  title('Dataset2')
 
data3=num2(1:11,5:6); % [X Y] for 3rd data set
 mdl3 = fitlm(data3(:,1),data3(:,2),'linear','VarNames',{'X3','Y3'})
 Rsq3 = mdl3.Rsquared.Ordinary % Goodness of fit
  figure(4) % Plotting the model
  plotResiduals(mdl3,'fitted')
 title('Dataset3')
 
data4=num2(1:11,7:8); % [X Y] for 4th data set
 mdl4 = fitlm(data4(:,1),data4(:,2),'linear','VarNames',{'X4','Y4'})
 Rsq4 = mdl4.Rsquared.Ordinary % Goodness of fit
  figure(5) % Plotting the model
 plotResiduals(mdl4,'fitted')
 title('Dataset4')