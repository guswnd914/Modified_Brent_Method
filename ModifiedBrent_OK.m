function [root,info] = ModifiedBrent_OK(func, Int, params)
%
% Modified Brent Method Algorithm
%
% Input:
%   [func]: Function for approximation.
%   [Int] : Interval that sets limit of approximation.
%   [params] : Conditions for tolerances and max iterations.
%       params.func_tol : if f(x) < func_tol then algorithm stops.
%       params.root_tol : if x < root_tol then algorithm stops.
%       params.maxit : if no. iteration > maxit, then algorithm fails. 
%
% Output:
%   [root] : Approximated point s.t. func(root) = 0.
%   [info] : Number of iterations and successibility.
%       info.flag <- 0 if succeed, 1 if failed.
%       info.it <- count of iterations
%
%------------------------------------------------------------------------
% >> func = @(x) 10e-12*(exp(x) - x^2 -2)
% >> Int = [-1, 3]
% >> params.root_tol = 1e-15; params.func_tol = 1e-15; params.maxit = 100;
% >> [root, info] = ModifiedBrent_OK(func, Int, params)
%
% root = 
%
%     1.3191
% 
% 
% info = 
% 
%   struct with fields:
% 
%       it: 6
%     flag: 0
%-------------------------------------------------------------------------
%
x0 = Int(1);                                        % Initials
x1 = Int(2);
fx0 = func(x0); 
fx1 = func(x1);
info.it = 0;
if (fx0*fx1) > 0 || params.maxit <= 0               % Validate root existence and iteration feasibility
    info.flag = 1; 
    root = false;
    return 
end                                                 
if abs(fx0) < abs(fx1)                              % swap(x0, x1) s.t. f(x0) > f(x1) always
    [x1, x0] = deal(x0, x1);
    fx1 = fx0;
end
x2 = x0;                                            % s.t. f(x0) > f(x1) < f(x2)
prev_x0 = x0;                                       % Previous points for bisection criteria
prev_x1 = x1; 
prev_fx = fx1;
step = 0;                                           % Successive interpolation steps                 
while abs(x1-x0) > params.root_tol && info.it <= params.maxit
    info.it = info.it + 1;                          % Repeat until convergence or above max iteration
    fx0 = func(x0); 
    fx1 = func(x1); 
    fx2 = func(x2);
    if fx0 ~= fx2 && fx1 ~= fx2                     
        L0 = (x0*fx1*fx2) / ((fx0-fx1) * (fx0-fx2));% Inverse quadtratic interporlation (IQI)
        L1 = (x1*fx0*fx2) / ((fx1-fx0) * (fx1-fx2));
        L2 = (x2*fx0*fx1) / ((fx2-fx0) * (fx2-fx1));
        new = L0 + L1 + L2;
    else                                            
        new = x1 - (fx1 * ((x1-x0) / (fx1-fx0)));   % Linear Interpolation (Secant Method)
    end
    fnew = func(new);
    if (fx0*fnew) < 0                               % Get opposite direction point
        tmp = x0;
    else
        tmp = x1;
    end
    step = step + (abs(new-tmp) > abs(prev_x1-prev_x0)/2);
    if step == 5 || abs(fnew) > abs(prev_fx)/2      % Bisection
        new = (x0+x1) / 2;
        fnew = func(new);
        step = 0;
    end
    if abs(fnew) < params.func_tol                  % Stop if value is converged
        break
    end
    if abs(fnew) < abs(fx1)                         % s.t. previous fx < current fx 
        prev_fx = fnew;
    end
    x2 = x1;
    if (fx0*fnew) < 0                               % s.t. f(x0) is always at opposite side of f(x1)
        x1 = new;
    else
        x0 = new;
    end
    if abs(fx0) < abs(fx1)                          % swap(x0, x1) s.t. f(x0) > f(x1) always
        [x1, x0] = deal(x0 , x1);
    end
    if step == 0                                    % Update previous points right after bisection
        prev_x0 = x0; prev_x1 = x1;
    end
end
info.flag = info.it > params.maxit;                 % Exceeded max iteration? Yes = 1 = Failed, No = 0 = Succeed
root = new;
