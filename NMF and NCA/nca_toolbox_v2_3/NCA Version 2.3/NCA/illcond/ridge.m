function [f_h, rss, f_h_ss, dof, varargout] = ridge(g, X, h)
% RIDGE  Ridge regression estimates.
%
%    Given a vector g, a design matrix X, and a regularization parameter
%    h, 
%           
%           [f_h, rss, f_h_ss, dof] = ridge(g, X, h) 
% 
%    returns the ridge regression estimate of the vector f in the
%    linear regression model
% 
%           g = X*f + noise.
%
%    Also returned are the residual sum of squares rss, the sum of
%    squares f_h_ss of the elements of the ridge regression estimate
%    f_h (the squared norm of f_h), and the effective number of
%    residual degrees of freedom dof.
%
%    If h is a vector of regularization parameters, the i-th column
%    f_h(:,i) is the ridge regression estimate for the regularization
%    parameter h(i); the i-th elements of rss and f_h_ss are the
%    associated residual sum of squares and estimate sum of squares.
%
%    If no regularization parameter h is given, generalized
%    cross-validation is used to determine the regularization
%    parameter. The chosen regularization parameter h and the value of
%    the GCV function are then returned as the fifth and sixth
%    output arguments 
%
%            [f_h, rss, f_h_ss, dof, h, G] = ridge(g, X);
% 
  
%    Adapted from various routines in Per Christian Hansen's
%    Regularization Toolbox.

  % Size of inputs
  [n, p]      = size(X); 
  q           = min(n, p);
  
  if nargin < 3
    nh        =  1;
  else
    nh        = length(h);
    % Possible choice of regularization parameter?
    if (min(h) < 0)
      error('Impossible regularization parameter h.')
    end
  end
    
  % Initialize outputs
  f_h         = zeros(p, nh);
  rss         = zeros(nh, 1); 
  f_h_ss      = zeros(nh, 1);
  dof         = zeros(nh, 1);
  
  % Compute SVD of X
  [U, S, V]   = svd(X, 0);  
  s           = diag(S);      % vector of singular values
  s2          = s.^2;
  
  % Coefficients in expansion of solution in terms of right singular
  % vectors
  fc          = U(:, 1:q)'*g;
  zeta        = s .* fc;

  % Determine regularization parameter by GCV if none is given
  if nargin < 3
    [h, G]       = gcv(U, s, g, 'ridge');
    varargout(1) = {h};
    varargout(2) = {G};
  end
    
  % Treat each regularization parameter separately.
  for j = 1:nh    
    f_h(:, j) = V(:, 1:q) * (zeta ./ (s2 + h(j)^2));
    f_h_ss(j) = sum(f_h(:, j).^2);
    rss(j)    = h(j)^4 * sum(fc.^2 ./ (s2 + h(j)^2).^2);
    dof(j)    = n - sum(s2./(s2 + h(j)^2));
  end
  
  % In overdetermined case, add rss of least-squares problem
  if (n > p)
    rss       = rss + sum((g - U(:, 1:q)*fc).^2);
  end
   

