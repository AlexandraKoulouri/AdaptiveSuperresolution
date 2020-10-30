function [x, iter] = STORM_Homotopy(A, y, varargin)
% PURPOSE:
% Solve optimization problem:
% minimize_x ||x||_1
% subject to ||Ax - y||_2 < tol
%---------------------------------------------------
% USAGE:
% x = STORM_Homotopy(A, im, 'tol', eps, 'isnonnegative', isnonnegative,...
%   'maxiteration', maxiteration);
%---------------------------------------------------
% INPUTS:
% A:                measurement matrix
% b:                measurements
% varargin:         parameters
%---------------------------------------------------
% OUTPUTS:
% x:                reconstructed x
% iter:             number of iterations
%---------------------------------------------------
% REFERENCE:
% Donoho, David L., and Yaakov Tsaig. "Fast solution of l1-norm
% minimization problems when the solution may be sparse."
% Information Theory, IEEE Transactions on 54, no. 11 (2008): 4789-4812.
%---------------------------------------------------

y = double(y(:));

% default parameter
maxiter = 5000;
isNonnegative = false;
tol = 100;

% parse the optional inputs
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin) / 2;

for parameterIndex = 1 : parameterCount
    parameterName = varargin{parameterIndex * 2 - 1};
    parameterValue = varargin{parameterIndex * 2};
    switch lower(parameterName)
        case 'maxiteration'
            maxiter = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        case 'tol'
            tol = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''...
                ' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin;

%% Initial setting of l1 homotopy
[~, N] = size(A);
% initial correlation
C = A' * y;
% I is active set
if isNonnegative
    [lambda, I]  = max(C);
    lambda = max(lambda, 0);
else
    [lambda, I] = max(abs(C));
end

x = zeros(N, 1);
iter = 0;

ATA = A(:, I)' * A(:, I);
iATA = inv(ATA);

out_i = [];

%% Iterate
while iter < maxiter
    iter = iter + 1;
    % J is inactive set
    J = setdiff(1 : N, union(I, out_i));
    % update direction
    sI = sign(C(I));
    d = zeros(N,1);
    d(I) = iATA * sI;
    v = A(:, I) * d(I);
    Atv = A' * v;
    % condition 1: residual correlation reaches lambda
    w = (lambda - C(J)) ./ (1 - Atv(J));
    gamma1 = min(w(w > 0));
    if ~isempty(gamma1)
        i1 = J(w == gamma1);
    else
        gamma1 = inf;
    end
    % condition 2: residual correlation reaches -lambda
    if isNonnegative
        gamma2 = inf;
    else
        w = (lambda + C(J)) ./ (1 + Atv(J));
        gamma2 = min(w(w > 0));
        if ~isempty(gamma2)
            i2 = J(w == gamma2);
        else
            gamma2 = inf;
        end
    end
    % condition 3: active coordinate reaches 0
    w = -x(I) ./ d(I);
    gamma3 = min(w(w > 0));
    if ~isempty(gamma3)
        i3 = I(w == gamma3);
    else
        gamma3 = inf;
    end
    % breakpoint is determined by the minimum gamma
    gamma = min([gamma1 gamma2 gamma3]);
    if isinf(gamma)
        break;
    end
    % update x
    x = x + gamma * d;
    % update lambda
    lambda = lambda - gamma;
    % update correlation
    C = C - gamma * Atv;
    C(I) = sign(C(I)) * lambda;
    % update active set
    if gamma == gamma1
        I = [I i1(1)];
        out_i = [];
    elseif gamma == gamma2
        I = [I i2(1)];
        out_i = [];
    elseif gamma == gamma3
        pos = find(I == i3);
        I([pos, end]) = I([end, pos]);
        I(end) = [];
        out_i = i3;
    end
    
    % determine keep going or not
    res = norm(A * x - y, 2);
    keep_going =  res > tol;
    
    if ~keep_going
        break;
    end
    
    % calculate (A_I' * A_I) ^ -1 with rank-1 update
    if gamma == gamma1 || gamma == gamma2
        % if an element is added
        v = A(:, I(end));
        u1 = A(:, I(1 : end - 1))' * v;
        u2 = iATA * u1;
        d = 1 / (v' * v - u1' * u2);
        u3 = d * u2;
        F = iATA + d * (u2 * u2');
        iATA = [F, -u3; -u3', d];
    else
        % if an element is removed
        iATA(:, [pos, end]) = iATA(:, [end, pos]);
        iATA([pos, end], :) = iATA([end, pos], :);
        F = iATA(1 : end - 1, 1 : end - 1);
        d = iATA(end, end);
        u = iATA(1 : end - 1, end);
        iATA = F - u * u' / d;
    end
end