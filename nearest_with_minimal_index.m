function [S,T,distance,time_seconds,Q,infotable] = nearest_with_minimal_index(A, B, min_index, maxiter, timemax_seconds, x0)
% Computes a (locally) nearest complex singular pencil T*x + S to the 
% square pencil B*x + A s.t. T*x + S is arbitrarily close to a pencil 
% having right minimal index 'min_index'. Generically, the right minimal 
% index of T*x + S is 'min_index'. 
% 
% Input:
%   A, B (matrix of size (n,n))
%       the coefficients of the pencil B*x + A in consideration
%   min_index (integer)
%        the right minimal index of the pencil T*x + S given as an output
%   maxiter (integer)
%        maximum amount of outer iterations, used in Manopt's
%        "trustregions" solver
%   timemax_seconds (integer)
%        maximum amount of time (in seconds) the algorithm can run before
%        the iteration is stopped
%   x0 (tensor of size (n,n,2))
%       contains the initial value of the optimization variables s.t. x0(:,:,1)
%       and x0(:,:,2) are unitary
%
% Output:
%   S, T (matrix of size (n,n))
%       the coefficients of the nearest singular pencil T*x + S with
%       minimal index 'min_index'
%   distance (vector)
%       the distance to the constructed singular pencil after each iteration
%   time_seconds (vector)
%       elapsed time after each iteration 
%   Q (tensor of size (n,n,2))
%       contains the final value of the optimization variables s.t. Q(:,:,1)
%       and Q(:,:,2) are unitary
%   infotable (table)
%       contains various types of information for diagnostic purposes
%
% Requirement: Manopt needs to be imported

n = length(B);

if not(exist('maxiter', 'var'))
    maxiter = 1000;
end
if not(exist('timemax_seconds', 'var'))
    timemax_seconds = 1000;
end
if not(exist('x0', 'var'))
    x0 = [];
end

% Rescale the pencil to be of norm 100
P_norm = norm([A B], 'f')*1e-2;
A = A / P_norm;
B = B / P_norm;

problem.M = stiefelcomplexfactory(n, n, 2);
problem.cost = @cost;
% The code uses the Euclidean gradient. Projection to 
% the tangent space of U(n) is handled automatically (see 
% stiefelcomplexfactory documentation)
problem.egrad = @egrad;
% Euclidean Hessian. Projection is handled automatically.
problem.ehess = @ehess;

options.tolgradnorm = 1e-10;
options.maxiter = maxiter;
options.maxtime = timemax_seconds;
options.verbosity = 2; % 2 = Default; 0 = No output; 

[Q, xcost, info, ~] = trustregions(problem, x0, options);

infotable = struct2table(info);
distance = sqrt(infotable.cost);
time_seconds = infotable.time;

% Construct the nearest singular pencil
T = triu(Q(:,:,1)*B*Q(:,:,2));
S = triu(Q(:,:,1)*A*Q(:,:,2));

T(min_index+1,min_index+1) = 0;
S(min_index+1,min_index+1) = 0;

T = Q(:,:,1)'*T*Q(:,:,2)';
S = Q(:,:,1)'*S*Q(:,:,2)';

% Squared distance to the singular pencil. Should be equal to xcost
assert(xcost -  (norm(B-T,'f')^2 + norm(A-S,'f')^2) < 1e-10)

% Rescale back
T = T*P_norm;
S = S*P_norm;
distance = distance*P_norm;

% ---------

function f = cost(Q)

    T = Q(:,:,1)*B*Q(:,:,2);
    S = Q(:,:,1)*A*Q(:,:,2);
    f = norm(tril(T,-1),'fro')^2 + norm(tril(S,-1),'fro')^2 + ...
        abs(T(min_index+1,min_index+1))^2 + abs(S(min_index+1,min_index+1))^2;
    
end

function g = egrad(Q)
    
    Q1 = Q(:,:,1);
    Q2 = Q(:,:,2);
    
    M11 = B*Q2;
    M01 = A*Q2;
    
    M12 = Q1*B;
    M02 = Q1*A;
    
    T = Q1*B*Q2;
    S = Q1*A*Q2;
    
    L1 = tril(T,-1);
    L0 = tril(S,-1);

    L1(min_index+1,min_index+1) = T(min_index+1,min_index+1);
    L0(min_index+1,min_index+1) = S(min_index+1,min_index+1);

    g = zeros(size(Q));
    g(:,:,1) = 2* L1 * M11' + 2* L0 * M01';
    g(:,:,2) = 2* M12' * L1 + 2* M02' * L0;

end

function H = ehess(Q, d)
    
    Q1 = Q(:,:,1);
    Q2 = Q(:,:,2);
    
    d1 = d(:,:,1);
    d2 = d(:,:,2);

    M11 = B*Q2;
    M01 = A*Q2;
    
    M12 = Q1*B;
    M02 = Q1*A;
    
    T = Q1*B*Q2;
    S = Q1*A*Q2;
    
    H = zeros(size(Q));

    L = tril(ones(size(Q1)),-1);
    L(min_index+1,min_index+1) = 1;

    L1 = L.*(d1*M11 + M12*d2);
    L0 = L.*(d1*M01 + M02*d2);

    H(:,:,1) = L1 * M11' + (L.*T) * d2' * B' ...
             + L0 * M01' + (L.*S) * d2' * A';

    H(:,:,2) = M12' * L1 + B' * d1' * (L.*T) ...
             + M02' * L0 + A' * d1' * (L.*S);

    % Scale by the omitted factor 2
    H = 2*H;

end

end
