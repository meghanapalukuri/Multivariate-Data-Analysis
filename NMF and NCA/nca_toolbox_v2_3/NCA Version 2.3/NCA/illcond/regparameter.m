% Matlab script that produces the figures in the handout on
% choosing regularization parmeters.

% The dataset used is that of homework 5.

g          = load('gravity.dat');
X          = load('integral.dat');
[n, p]     = size(X);    
q          = min(n, p);   % maximum number of nonzero singular values
noise_std  = .1;          % standard deviation of "measurement" noise

% exact solution
t          = linspace(0, 1, n);
f          = sin(pi*t) + 0.5 * sin(2*pi*t);

% singular value decomposition of integral operator
[U,s,V]    = svd(X,0);
s          = diag(s);     % vector of singular values (SVD returns matrix)

% assemble vector of regularization parameters
nh         = 200;         % number of regularization parameters
h_max      = s(1);        % maximum reg. par. = max(singular value)
h_min      = (eps)^(1/3); % minimum reg. par. (dependent on machine accuracy)
h          = h_max * (h_min/h_max).^( [0:nh-1] / (nh-1) );

% compute ridge regression estimates associated with regularization
% parameters in vector h
[f_h, rss_h f_h_ss, dof_h] = ridge(g, X, h);

% get regularization parameter from discrepancy principle
[d, idis]  = min(abs( rss_h - n* noise_std^2 ));
h_dis      = h(idis);

disp(['From discrepancy principle: h = ', num2str(h_dis)])

% plot corresponding estimate of density variations (along with
% exact solution)
set(gcf, 'DefaultAxesTickDir', 'out') %  property of all plots
figure(1)
clf
plot(t, f_h(:, idis), 'kx', ...
     t, f, 'k-')
set(gca, 'Box', 'off', ... 
         'YTick', [0:.5:1.5], ...
	 'XTick', [0:.5:1])
axis([0 1 0 1.5])
xlabel('t')
ylabel('Density variation')

% plot the L-curve (mark the points where the regularization parameter
% equals one of the singular values, and mark the corner, which was
% found by inspection of the curve)

h_l     = 0.14;           % corner of L-curve

% get corresponding estimate of density variations
[f_hl, rss_hl, f_hl_ss] = ridge(g, X, h_l);

% find indices of vector of regularization parameters where h=s(k)
ihs     = zeros(q, 1);
for k = 1:q
  [dum, i] = min(abs( s(k)-h ));
  ihs(k)   = i;           % index where h = singular value
end

figure(2)
clf

loglog(rss_h, f_h_ss, 'k-', ...
       rss_h(ihs), f_h_ss(ihs), 'kx')
hold on

% label the singular values
for k = 1:q   
  hd = text(rss_h(ihs(k)), f_h_ss(ihs(k)), num2str(h(ihs(k))));
  set(hd, 'VerticalAlignment', 'bottom')
end

% mark corner of L-curve
loglog([1e-3 rss_hl], [f_hl_ss f_hl_ss], 'k--', ...
       [rss_hl rss_hl], [1 f_hl_ss], 'k--')
set(gca, 'Box', 'off')
axis([1e-3  2e2  1 2e4])
xlabel('Residual SS')
ylabel('Estimate SS')

disp(['From L-curve: h = ', num2str(h_l)]);

% Plot estimate of density variations corresponding to corner of L-curve
figure(3)
clf
plot(t, f_hl, 'kx', ...
     t, f, 'k-')
set(gca, 'Box', 'off', ...
	 'YTick', [0:.5:1.5], ...
	 'XTick', [0:.5:1])
axis([0 1 0 1.5])
xlabel('t')
ylabel('Density variation')

% plot the generalized cross-validation function for the
% regularization parameters in vector h
G       = rss_h ./ dof_h.^2;   % GCV function

% find minimum of GCV function and corresponding estimate of
% density variations
[f_hgcv, rss_hgcv, f_hgcv_ss, dof_hgcv, h_gcv, G_hgcv] = ...
    ridge(g, X);

disp(['From GCV: h = ', num2str(h_gcv)]);

figure(4)
clf
loglog(h, G, 'k-', ...
       [h_gcv h_gcv], [10^(-3) G_hgcv], 'k--')
set(gca, 'Box', 'off')
xlabel('Regularization parameter h')
ylabel('GCV function')
axis([1e-3 5 1e-3 1])

% plot estimate of density variations corresponding to GCV
% regularization parameter
figure(5)
clf

plot(t, f_hgcv, 'kx', ...
     t, f, 'k-')
set(gca, 'Box', 'off', ...
	 'YTick', [0:.5:1.5], ...
	 'XTick', [0:.5:1])
axis([0 1 0 1.5])
xlabel('t')
ylabel('Density variation')






