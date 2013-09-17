function [x, resnorm, evalf, output, jacob] = newtonraphson(fun, x0, options)
% NEWTONRAPHSON Solve set of non-linear equations using Newton-Raphson method.
%
% FUN is a function that returns a vector of residuals of the governing
% equations and takes a vector, x, as its only argument. When the equations are
% solved by x, then F(x) == zeros(size(x(:), 1)). If FUN is an anonymous
% function it must not use DEAL, because there's no way to check NARGOUT.
%
% Optionally FUN may return the Jacobian, Jij = dFi/dxj, as an additional
% output. If FUN only returns one output, then J is estimated using a center
% difference approximation,
%   Jij = dFi/dxj = (Fi(xj + dx) - Fi(xj - dx))/2/dx.
% The Jacobian should be a square matrix whose columns correspond to d/dxj and
% whose rows correspond to dFi/d.
% EG:  J23 = dF2/dx3 is the 2nd row ad 3rd column.
%
% X0 is a vector of initial guesses.
%
% OPTIONS is a structure of solver options created using OPTIMSET.
% EG: options = optimset('TolX', 0.001).
%
% The following options can be set:
% * OPTIONS.TOLFUN is the maximum tolerance of the norm of the residuals
%   [1e-6].
% * OPTIONS.MAXITER is the maximum number of iterations before giving up
%   [100]
% * OPTIONS.DISPLAY sets the level of display: {'off', 'iter'}
%   ['iter']
%
% X is the solution that solves the set of equations within the given tolerance.
% RESNORM is norm(FEVAL) and FEVAL is F(X). OUTPUT is a structure
% containing the number of iterations and JACOB is the J(X).
%
% See also OPTIMSET, OPTIMGET, FMINSEARCH, FZERO, FMINBND, FSOLVE, LSQNONLIN
%
% References:
% http://en.wikipedia.org/wiki/Newton's_method
% http://en.wikipedia.org/wiki/Newton's_method_in_optimization
%
% Version 0.2
% * Remove `options.FinDiffRelStep` and `options.TypicalX` since not in MATLAB.
% * Set `dx = eps^(1/3)` in `jacobian` function.
% * Remove `options` argument from `funwrapper` & `jacobian` functions
%   since no longer needed.
% * Set typx = x0; typx(x0==0) = 1; % use initial guess as typx, if 0 use 1.
% * Replace `feval` with `evalf` since `feval` is builtin.
%% initialize
% There are no argument checks!
x0 = x0(:); % needs to be a column vector
% set default options
oldopts = optimset( ...
    'TolFun', 1e-6, 'MaxIter', 100, 'Display', 'iter');
if nargin<3
    options = oldopts; % use defaults
else
    options = optimset(oldopts, options); % update default with user options
end
F = @(x)funwrapper(fun, x); % wrap FUN so it always returns J
%% get options
tol = optimget(options, 'TolFun'); % tolerance
maxiter = optimget(options, 'MaxIter'); % max number of iterations
display = strcmpi('iter', optimget(options, 'Display')); % display iterations
typx = x0; % x scaling value is the initial guess
typx(x0==0) = 1; % remove zeros from typical x
%% set scaling values
Ftyp = F(typx); % typical f
J0 = Ftyp*(1./typx'); % Jacobian scaling matrix
%% set display
if display
    fprintf('\n%10s %10s %10s %12s\n', 'Niter', 'resnorm', 'stepnorm', ...
        'convergence')
    for n = 1:45,fprintf('-'),end,fprintf('\n')
    printout = @(n, r, s, c)fprintf('%10d %10.4g %10.4g %12.4g\n', n, r, s, c);
end
%% check initial guess
x = x0; % initial guess
[evalf, J] = F(x); % evaluate initial guess
resnorm = norm(evalf); % calculate norm of the residuals
%% solver
Niter = 0; % start counter
if display,printout(Niter, resnorm, 0, 0);end
while resnorm>tol && Niter<maxiter
    Niter = Niter+1; % increment counter
    Jstar = J./J0; % scale Jacobian
    Fstar = evalf./Ftyp; % scale F
    dx_star = -Jstar\Fstar; % Newton-Raphson solver
    dx = dx_star.*typx; % rescale x
    x = x+dx; % next guess
    [evalf, J] = F(x); % evaluate residuals
    resnorm0 = resnorm; % old resnorm
    resnorm = norm(evalf); % calculate new resnorm
    convergence = log(resnorm0/resnorm); % calculate convergence rate
    stepnorm = norm(dx); % norm of the step
    if display,printout(Niter, resnorm, stepnorm, convergence);end
end
%% output
output.iterations = Niter;
jacob = J;
end

function [evalf, J] = funwrapper(fun, x)
% if nargout<2 use finite differences to estimate J
try
    [evalf, J] = fun(x);
catch
    evalf = fun(x);
    J = jacobian(fun, x); % evaluate center diff if no Jacobian
end
evalf = evalf(:); % needs to be a column vector
end

function J = jacobian(fun, x)
% estimate J
dx = eps^(1/3); % finite difference delta
dof = size(x, 1); % degrees of freedom
J = zeros(dof); % a square matrix
for n = 1:dof
    % create a vector of deltas, change delta_n by dx
    delta = zeros(dof, 1); delta(n) = delta(n)+dx;
    dF = fun(x+delta)-fun(x-delta); % delta F
    J(:, n) = dF/dx/2; % derivatives dF/d_n
end
end
