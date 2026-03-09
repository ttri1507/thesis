function A = FUNC_NearestClustering(MODEL, numSAT, numUAV, numUE, N_U)
numUAVeachSAT = numUAV/numSAT;
numUEeachSAT = numUE/numSAT;
Distance = sqrt( (MODEL.UAV(1,:)'-MODEL.UE(1,:)).^2 + ...
                 (MODEL.UAV(2,:)'-MODEL.UE(2,:)).^2 + ...
                 (MODEL.UAV(3,:)'-MODEL.UE(3,:)).^2 );
% Initialize A: using diagonal line assign 1
A = zeros(numUAV,numUE);
for s = 1:4
    Us = numUAVeachSAT*(s-1) + 1:numUAVeachSAT*s;
    Ks = numUEeachSAT*(s-1) + 1:numUEeachSAT*s;
    for k = 1:numUEeachSAT
        % Find no-overload UAVs
        SumUs = sum(A(Us,:),2); IdU_ok = find((SumUs<N_U)==1);
        Us_ok = Us(IdU_ok); 
        [dist_min, id_min] = min(Distance(Us_ok,k));
        indUAV = Us_ok(id_min);
        indUE = Ks(k);
        A(indUAV,indUE) = 1;
    end
end
end

