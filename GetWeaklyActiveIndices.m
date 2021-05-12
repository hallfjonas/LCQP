function[indices] = GetWeaklyActiveIndices(x, L, R, tol)

sqrtTol = sqrt(tol);
indices = [];

for i = 1:size(L,1)
    if (L(i,:)*x <= sqrtTol && R(i,:)*x <= sqrtTol)
        indices = [indices, i];
    end
end

end
