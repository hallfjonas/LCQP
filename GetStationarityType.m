function [ val, msg ] = GetStationarityType(x, lam, L, R, L_indices, R_indices, tol)

W_indices = GetWeaklyActiveIndices(x, L, R, tol);

% Flags for types
isSstationary = true;
isMstationary = true;

% Check M-,C-,W-stationarity
for i = 1:length(W_indices)
    dualProd = lam(L_indices(W_indices(i)))*lam(R_indices(W_indices(i)));
    dualMin = min(L_indices(W_indices(i)), lam(R_indices(W_indices(i))));
    
    if (dualMin < 0)
        isSstationary = false;
    end
    
    % Failure of M-Stationarity => C-,W-stationarity
    if ( abs(dualProd) >= tol && dualMin <= 0 )
        
        % Failure of C-stationarity => W-stationarity
        if (dualProd <= tol)
            val = 3;
            msg = fprintf("W-Stationary point found\n");
            return;
        end
        
        isMstationary = false;        
    end
end

% Strongly stationary point
if (isSstationary)
    val = 0;
    msg = fprintf("S-Stationary point found\n");
    return;    
end

% Mordukhovich stationary point
if (isMstationary)
    val = 1;
    msg = fprintf("M-Stationary point found\n");
    return;
end

% If M-stationarity failed, then must be C-stationary
val = 2;
msg = fprintf("C-Stationary point found\n");

end
