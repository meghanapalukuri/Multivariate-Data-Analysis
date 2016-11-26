%% Question 1.a
clear
clc
load('Inorfull.mat');

stdavg = [];
Ydata = [];
Conc = [];
rand('seed',45)
for i=1:26
    dt_num = randi(5);
    istart = 5*(i-1)+dt_num;
    stdavg = [stdavg; stdDATA(istart,:)];
    Ydata = [Ydata; DATA(istart,:)];
    Conc = [Conc;CONC(istart,:)];
end

Linv = inv(diag(mean(stdavg)));
Ydata =Ydata*Linv;

nfact = 3;
conc = Conc;
[nsamples, nwave] = size(Ydata);


% Perform PCA to obtain scores matrix
[U S V] = svd(Ydata,'econ');
T = U(:,1:nfact)*S(1:nfact,1:nfact);

A = T*V(:,1:nfact)';
A = A*(diag(mean(stdavg)));

k = 3;

[W,H,iter,HIS]=nmf(A,k);
pure = [PureCo ;PureCr; PureNi];
correlation_mat = zeros(3,3);
for i = 1:3
    for j = 1:3
        correlation_mat(i,j) = corr2(pure(i,:),H(j,:));
    end
end
correlation_mat