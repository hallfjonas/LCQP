function [ alpha ] = GetOptimalStepLength(C, g, Rk, Q, pk, xk)

% Outer loop objective's Hessian
Qk = Q + Rk*C;

% Quadratic term
qk = pk'*Qk*pk;

% Linear term
lk = pk'*(Qk*xk + g);

% Obtain optimal step length
if (qk > eps && lk < 0)
    alpha = min(-lk/qk, 1);
else
    alpha = 1;
end

end