function [x, resnorm, F, output, jacob] = newtonraphson(fun, x0, options)
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
% * http://en.wikipedia.org/wiki/Newton's_method
% * http://en.wikipedia.org/wiki/Newton's_method_in_optimization
% * 9.7 Globally Convergent Methods for Nonlinear Systems of Equations 383,
%   Numerical Recipes in C, Second Edition (1992),
%   http://www.nrbook.com/a/bookcpdf.php

% Version 0.3
% * Display RCOND each step.
% * Replace nargout checking in funwrapper with ducktypin.
% * Remove Ftyp and F scaling b/c F(typx)->0 & F/Ftyp->Inf!
% * User Numerical Recipies minimum Newton step, backtracking line search
%   with alpha = 1e-4.
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
FUN = @(x)funwrapper(fun, x); % wrap FUN so it always returns J
%% get options
TOL = optimget(options, 'TolFun'); % tolerance
MAXITER = optimget(options, 'MaxIter'); % max number of iterations
DISPLAY = strcmpi('iter', optimget(options, 'Display')); % display iterations
TYPX = x0; % x scaling value is the initial guess
TYPX(x0==0) = 1; % remove zeros from typical x
ALPHA = 1e-4; % criteria for decrease
MIN_LAMBDA = 0.1; % min lambda
MAX_LAMBDA = 0.5; % max lambda
%% set scaling values
% TODO: let user set weights
weight = ones(size(x0));
J0 = weight*(1./TYPX'); % Jacobian scaling matrix
%% set display
if DISPLAY
    fprintf('\n%10s %10s %10s %10s %10s %12s\n', 'Niter', 'resnorm', 'stepnorm', ...
        'lambda', 'rcond', 'convergence')
    for n = 1:67,fprintf('-'),end,fprintf('\n')
    fmtstr = '%10d %10.4g %10.4g %10.4g %10.4g %12.4g\n';
    printout = @(n, r, s, l, rc, c)fprintf(fmtstr, n, r, s, l, rc, c);
end
%% check initial guess
x = x0; % initial guess
[F, J] = FUN(x); % evaluate initial guess
Jstar = J./J0; % scale Jacobian
rc = rcond(Jstar); % reciprocal condition
resnorm = norm(F); % calculate norm of the residuals
%% solver
Niter = 0; % start counter
lambda = 1; % backtracking
if DISPLAY,printout(Niter, resnorm, 0, 1, rc, 0);end
while (resnorm>TOL && Niter<MAXITER) || lambda<1
    if lambda==1
        %% Newton-Raphson solver
        Niter = Niter+1; % increment counter
        dx_star = -Jstar\F; % calculate Newton step
        % NOTE: use isnan(f) instead of STPMAX
        dx = dx_star.*TYPX; % rescale x
        g = F'*Jstar; % gradient of resnorm
        slope = g*dx_star; % slope of gradient
        fold = F'*F; % objective function
        xold = x; % initial value
    end
    x = xold+dx*lambda; % next guess
    [F, J] = FUN(x); % evaluate next residuals
    Jstar = J./J0; % scale next Jacobian
    f = F'*F; % new objective function
    %% check for convergence
    lambda1 = lambda; % save previous lambda
    if f>fold+ALPHA*lambda*slope
        if lambda==1
            lambda = -slope/2/(f-fold-slope); % calculate lambda
        else
            A = 1/(lambda1 - lambda2);
            B = [1/lambda1^2,-1/lambda2^2;-lambda2/lambda1^2,lambda1/lambda2^2];
            C = [f-fold-lambda1*slope;f2-fold-lambda2*slope];
            coeff = num2cell(A*B*C);
            [a,b] = coeff{:};
            if a == 0
                lambda = -slope/2/b;
            else
                discriminant = b^2 - 3*a*slope;
                if discriminant<0
                    lambda = MAX_LAMBDA*lambda1;
                elseif b<=0
                    lambda = (-b+sqrt(discriminant))/3/a;
                else
                    lambda = -slope/(b+sqrt(discriminant));
                end
            end
            lambda = min(lambda,MAX_LAMBDA*lambda1); % minimum step length
        end
        lambda2 = lambda1;f2 = f; % save 2nd most previous value
        lambda = max(lambda,MIN_LAMBDA*lambda1); % minimum step length
        continue
    elseif isnan(f)
        lambda = MAX_LAMBDA*lambda1;
        continue
    else
        lambda = 1; % fraction of Newton step
    end
    %% display
    resnorm0 = resnorm; % old resnorm
    resnorm = norm(F); % calculate new resnorm
    convergence = log(resnorm0/resnorm); % calculate convergence rate
    stepnorm = norm(dx); % norm of the step
    rc = rcond(Jstar); % reciprocal condition
    if DISPLAY,printout(Niter, resnorm, stepnorm, lambda1, rc, convergence);end
end
%% output
output.iterations = Niter;
jacob = J;
end

function [F, J] = funwrapper(fun, x)
% if nargout<2 use finite differences to estimate J
try
    [F, J] = fun(x);
catch
    F = fun(x);
    J = jacobian(fun, x); % evaluate center diff if no Jacobian
end
F = F(:); % needs to be a column vector
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
