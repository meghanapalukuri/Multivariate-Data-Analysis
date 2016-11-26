clc
clear all

Data=load('foundry.txt');
Estimation_Data=Data(1:300,:);
Est_X=Estimation_Data(:,1:(size(Estimation_Data,2)-1));
Est_y=Estimation_Data(:,size(Estimation_Data,2));
Validation_Data=Data(301:size(Data,1),:);
Validation_X=Validation_Data(:,1:(size(Estimation_Data,2)-1));
Validation_y=Validation_Data(:,size(Estimation_Data,2));

% Pre-processing the data
% Scaling using range for each variable - The same line ! We're just
% shifting the data.
RangeEst=range(Est_X,1);
RangeVal=range(Validation_X,1);
for i=1:size(Est_X,1)
Pre1_Est_X(i,:)=(Est_X(i,:)-min(Est_X,[],1))./RangeEst;
end

for i=1:size(Validation_X,1)
Pre1_Validation_X(i,:)=(Validation_X(i,:)-min(Validation_X,[],1))./RangeVal;
end


%% Standard OLS Regression : model - using 300 data points and testing the model with the remaining

% Fitting linear model using standard ols (without weights)
mdl = fitlm(Est_X,Est_y)
% Plotting the model
figure(1)
plot(mdl)
CoeffOfDetermination = mdl.Rsquared.Ordinary % R squared Value

%Validating the model using TestData
Predicted_y=feval(mdl,Validation_X);
Pred_Error=Predicted_y-Validation_y;
MaxAbs_Pred_Error=max(abs(Pred_Error))
Standard_Deviation = std(Pred_Error)
figure(2)
hist(Pred_Error,15)
title('Validation Data Error Plot')
xlabel('Predicted value error')
ylabel('Frequency')

% Fitting linear model using standard ols (without weights)
mdl2 = fitlm(Pre1_Est_X,Est_y)
% Plotting the model
figure(3)
plot(mdl2)
CoeffOfDetermination = mdl2.Rsquared.Ordinary % R squared Value

%Validating the model using TestData
Predicted_y=feval(mdl2,Pre1_Validation_X);
Pred_Error=Predicted_y-Validation_y;
MaxAbs_Pred_Error=max(abs(Pred_Error))
Standard_Deviation = std(Pred_Error)
figure(4)
hist(Pred_Error,15)
title('Validation Data Error Plot')
xlabel('Predicted value error')
ylabel('Frequency')

