function [ primalSteps, dualSteps, fVals, RVals, exitflag, auxOutput, iter, k ] = LCQP(Q, g, L, R, A, lbA, ubA, lb, ub, params)
%% Input Consistency
[nx, nc, ncomp] = PerformInputConsistencyTest(Q, g, A, lb, ub, lbA, ubA, L, R);

%% Parameter handling
% qpOASES mode ('mpc' could be faster but less precise)
options = qpOASES_options('default');
if (isfield(params, 'qpOASES_options'))
    options = params.qpOASES_options;
end

% Tolerence for complementarity violation
complementarityTolerance = 1e-10;
if (isfield(params, 'complementarityTolerance'))
    complementarityTolerance = params.complementarityTolerance;
    assert(complementarityTolerance > 0);
end

% Tolerence for stationarity violation
stationarityTolerance = options.terminationTolerance*10;
if (isfield(params, 'stationarityTolerance'))
    stationarityTolerance = params.stationarityTolerance;
    assert(stationarityTolerance > 0);
end

% Initial starting point
xk = zeros(nx,1);
if (isfield(params, 'x0'))
    % Initialize primals
    xk = params.x0;
    assert(all(size(xk) == [nx, 1]));
end

% Store all steps
storeSteps = false;
if (isfield(params, 'storeSteps'))
    storeSteps = params.storeSteps;
end

% Initial penalty parameter
Rk = 1;
if (isfield(params, 'R0'))
    Rk = params.R0;
    assert(Rk > 0);
end

% Solve initial QP with or without penalization
solveZeroPenaltyFirst = true;
if (isfield(params, 'solveZeroPenaltyFirst'))
    solveZeroPenaltyFirst = params.solveZeroPenaltyFirst;
    assert(islogical(solveZeroPenaltyFirst));
end

% Penalty update rule
penaltyUpdater = @(R) 2*R;
if (isfield(params, 'penaltyUpdater'))
    penaltyUpdater = params.penaltyUpdater;
end

% Penalty parameter value required to break the iteration.
Rbreak = 1000000;
if (isfield(params, 'Rbreak'))
    Rbreak = params.Rbreak;
    assert(Rbreak > 0);
end

% Maximum number of total iterations.
maxIter = 100000;
if (isfield(params, 'maxIter'))
    maxIter = params.maxIter;
    assert(maxIter > 0);
end

% Print output
printStats = false;
if (isfield(params, 'printStats'))
    printStats = params.printStats;
    assert(islogical(printStats));
end

% Switch to sparse matrices
useSparseMatrices = false;
if (isfield(params, 'useSparseMatrices'))
    useSparseMatrices = params.useSparseMatrices;
    assert(islogical(useSparseMatrices));
end

% Perform gradient perturbation method
useGradientPerturbation = false;
if (isfield(params, 'useGradientPerturbation'))
    useGradientPerturbation = params.useGradientPerturbation;
    assert(islogical(useGradientPerturbation));
end

% Retrieve stationarity type, i.e. S, M, C, W 
retrieveStationarityType = false;
if (isfield(params, 'retrieveStationarityType'))
    retrieveStationarityType = params.retrieveStationarityType;
    assert(islogical(retrieveStationarityType));
end

%% Data preperation

% Append rows of L and R to A
A = [A; L; R];
lbA = [lbA; zeros(2*ncomp,1)];
ubA = [ubA; inf(2*ncomp,1)];

% Add unbounded box constraints if none are passed
if (isempty(lb))
    lb = -inf(nx,1);
end
if (isempty(ub))
    ub = inf(nx,1);
end

% Switch to sparse matrices
if (useSparseMatrices)
    Q = sparse(Q);
    A = sparse(A);
    L = sparse(L);
    R = sparse(R);
end

%% Helpers and initializations
% Outer loop index
k = 0;

% Inner loop index
i = 0;

% Total iterate counter
iter = 0;

% Originial objective
F = @(primal) 1/2*primal'*Q*primal + g'*primal;

% Symmetric penalty matrix
C = L'*R + R'*L;

% Penalty function 
phi = @(primal) 1/2*primal'*C*primal;
grad_phi = @(primal) C*primal;

% Stationarity and complementarity
Stationarity = @(primal, pen, dual_xk, dual_gk) Q*primal + g + pen*grad_phi(primal) - dual_xk - A'*dual_gk;

% Merit function
merit_fun = @(primal, pen) F(primal) + pen*phi(primal);

% Initial return values
primalSteps = xk;
RVals = Rk;
fVals = [];
dualSteps = [];

%% Initiate qpOASES homotopy
if (solveZeroPenaltyFirst)
    [QP,xk,~,exitflag,iterQPO,lk,auxOutput] = qpOASES_sequence('i', Q, g, A, lb, ub, lbA, ubA, options, xk);
else
    [QP,xk,~,exitflag,iterQPO,lk,auxOutput] = qpOASES_sequence('i', Q, g + Rk*grad_phi(xk), A, lb, ub, lbA, ubA, options);
end
    
% return on failed QP solve
if (exitflag ~= 0)
    qpOASES_sequence('c', QP);
    return;
end

% Begin with full step
alpha = 1;
pk = xk - params.x0;

% Get duals
lam_xk = lk(1:length(xk));
lam_gk = lk((length(xk)+1):end);

% Main Loop
while ( true )        
    i = i + 1;

    iter = iter + 1;

    %% Store steps and update output
    if (storeSteps)
        primalSteps = [primalSteps, xk];
        dualSteps = [dualSteps, lk];
        fVals = [fVals, F(xk)];
        RVals = [RVals; Rk];
    else
        primalSteps = xk;
        dualSteps = lk;
        fVals = F(xk);
        RVals = Rk;
    end    


    %% Termination checks       
    % Evaluate stationarity and complementarity
    stat = norm(Stationarity(xk, Rk, lam_xk, lam_gk));

    % Iteration stats
    if (printStats)
        compl = abs(phi(xk));
        PrintStats(k, i, stat, compl, Rk, alpha, merit_fun(xk, Rk), norm(pk), iterQPO);
    end
    
    % Need this alternative stationarity check due to qpOASES termination
    % criterion
    if (stat <= stationarityTolerance)
        if (Rk > Rbreak)
            fprintf("Exceeded maximum penalty value!\n"); 
            qpOASES_sequence('c', QP);
            return;
        end
        compl = abs(phi(xk));
        if (compl < complementarityTolerance)
            qpOASES_sequence('c', QP);       
            
            L_indices = nx+nc+1:nx+nc+ncomp;
            R_indices = nx+nc+ncomp+1:nx+nc+2*ncomp;
            lk(L_indices) = lk(L_indices) - Rk*R*xk;
            lk(R_indices) = lk(R_indices) - Rk*L*xk;         
            
            if (retrieveStationarityType)
                GetStationarityType(xk, lk, L, R, L_indices, R_indices, complementarityTolerance);
            end
            
            return;
        end
        % Update penalty parameter
        k = k + 1;
        i = 0;
        Rk = penaltyUpdater(Rk);
    end
    if (iter > maxIter)
        disp('Maximum number of iterations reached. Exiting without convergence');
        
        exitflag = 2;
        qpOASES_sequence('c', QP);
        return;
    end
    
    %% Step computation
    % Update linearization
    d = grad_phi(xk);
    
    % Add perturbation if desired (could help avoid spurious solutions)
    if (useGradientPerturbation)
         d = d + (rand(nx, 1) - 0.5)*4*eps;
    end

    % Solve convex QP
    [x_new,~,exitflag,iterQPO,lambda,auxOutput] = qpOASES_sequence('h', QP, g + Rk*d, lb, ub, lbA, ubA, options);      
    
    % Terminate if qpOASES fails
    if (exitflag ~= 0)
        qpOASES_sequence('c', QP);
        return;
    end

    %% Step length computation
    pk = x_new - xk;
    alpha = GetOptimalStepLength(C, g, Rk, Q, pk, xk);
        
    %% Step update
    xk = xk + alpha*pk;
    lam_xk = lambda(1:nx);
    lam_gk = lambda((nx+1):end);
    lk = [lam_xk; lam_gk];       
end
end