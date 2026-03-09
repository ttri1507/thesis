function A = FUNC_RandomClustering(numSAT, numUAV, numUE, N_U)
numUAVeachSAT = numUAV/numSAT;
numUEeachSAT = numUE/numSAT;

% Initialize A: using diagonal line assign 1
A = zeros(numUAV,numUE);
for s = 1:4
    Us = numUAVeachSAT*(s-1) + 1:numUAVeachSAT*s;
    Ks = numUEeachSAT*(s-1) + 1:numUEeachSAT*s;
    for k = 1:numUEeachSAT
        % Find no-overload UAVs
        SumUs = sum(A(Us,:),2); IdU_ok = find((SumUs<N_U)==1);
        idu = randi(length(IdU_ok)); u = IdU_ok(idu);
        indUAV = Us(u);
        indUE = Ks(k);
        A(indUAV,indUE) = 1;
    end
end
end

