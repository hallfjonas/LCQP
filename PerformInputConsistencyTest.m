function [nx, nc, ncomp] = PerformInputConsistencyTest(Q, g, A, lb, ub, lbA, ubA, L, R)

[nQ, mQ] = size(Q);
[ng, mg] = size(g);
[nA, mA] = size(A);
[nlb, mlb] = size(lb);
[nub, mub] = size(ub);
[nlbA, mlbA] = size(lbA);
[nubA, mubA] = size(ubA);
[nL, mL] = size(L);
[nR, mR] = size(R);

% Assertions on objective
assert(nQ == mQ && mQ > 0);
nx = nQ;
assert(ng == nx && mg == 1);

% Assertions on A
nc = nA;
assert(nlbA == nc);
assert(nubA == nc);
if (nc > 0)    
    assert(mA == nx);
    assert(mlbA == 1);
    assert(mubA == 1); 
end

% Bound assertions
assert(nlb == nx && mlb == 1);
assert(nub == nx && mub == 1);

% Complementarity assertions
assert(nL > 0 && mL == nx);
ncomp = nL;
assert(nR == ncomp && mR == nx);

end