function [] = PrintStats(k, i, stat, compl, R, alphak, Mk, norm_pk, iterQPO)

if (mod(i, 10) == 1)
    fprintf("------|-------|--------------|-------------|---------|-------------|-----------|-----------|--------\n");
    fprintf("    k |     i | stationarity | complement. | penalty | step length | merit val |  norm pk  | qpIter \n");
    fprintf("------|-------|--------------|-------------|---------|-------------|-----------|-----------|--------\n");
end
fprintf("%5d | %5d | %12.2g | %11.2g | %7.2g | %11.2g | %9.6g | %9.2g | %6d \n", k, i, stat, compl, R, alphak, Mk, norm_pk, iterQPO);
end