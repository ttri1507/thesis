function [P_sat, P_uav_uav, P_uav_ue, ASU, AUU] = FUNC_PA(MODEL, A, B, ...
numSAT, numUAV, numUE, numUAVeachSAT, numUEeachUAV, PmaxSAT, PmaxUAV)
% Power allocation for SATs
ASU = zeros(numSAT, numUAV);
% Power allocation for UAVs-to-UAVs
AUU = zeros(numUAV, numUAV);
for k = 1:numUE
    f = MODEL.RP(k);
    s = ceil(k/(numUAVeachSAT*numUEeachUAV));
    uk = find(A(:,k)==1);
    Us = numUAVeachSAT*(s-1) + 1:numUAVeachSAT*s;
    % If package f isn't prestored at any UAV, f must be transmited from s
    ASU(s,uk) = ASU(s,uk) + ~any(B(Us, f) == 1);
    
    % If package f isn't prestored at primary UAV, but it is in some secondary
    % UAV
    if((B(uk, f)==0)&&any(B(Us, f) == 1))
        % Find the secondary UAV having h max and lowest number of serving
        % UAVs xxxx In this code focuse on h max only
        Up = numUAV/4*(s-1) + find(B(Us,f)==1);
        H_norm_Up_uk = abs(MODEL.H_UU(Up, uk));
        [h_max, id_umax] = max(H_norm_Up_uk);
        AUU(Up(id_umax),uk) = AUU(Up(id_umax),uk) + 1;
    end
end
% The matrix after multiplying channel matrix and precoding
H_sqSU = abs(sum( conj(MODEL.H_SU).*MODEL.PrecodingS, 3)).^2;
WeighS = ASU./H_sqSU;
DenominatorS = sum(WeighS, 2); 
% Delete zero value
R0S = find(DenominatorS==0); DenominatorS(R0S) = 1;
DenominatorS = repmat(DenominatorS,1,numUAV);
P_sat = (PmaxSAT*ones(numSAT, numUAV))./DenominatorS.*WeighS;
P_sat(R0S,:) = 0;

H_sqUU = abs(MODEL.H_UU).^2;
WeighU = AUU./H_sqUU;
DenominatorU = sum(WeighU, 2);
% Delete zero value
R0U = find(DenominatorU==0); DenominatorU(R0U) = 1;
DenominatorU = repmat(DenominatorU,1,numUAV);
P_uav_uav = (PmaxUAV*0.25*ones(numUAV, numUAV))./DenominatorU.*WeighU;
P_uav_uav(R0U,:) = 0;
    
% Power allocation for UAVs-to-UEs
AUK = A;
H_sqUK = abs(MODEL.H_UK).^2;
WeighUK = AUK./H_sqUK;
DenominatorUK = sum(WeighUK, 2);
% Delete zero value
R0UK = find(DenominatorUK==0); DenominatorUK(R0UK) = 1;
DenominatorUK = repmat(DenominatorUK,1,numUE);
P_uav_ue = (PmaxUAV*0.75*ones(numUAV, numUE))./DenominatorUK.*WeighUK;
P_uav_ue(R0UK,:) = 0;
end