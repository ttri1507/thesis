function [MODEL] = FUNC_CoalitionGame(MODEL,B,t_thres,theta,Area,IdClusUE_RC,NumUE_in_SAT_RC)
numSAT = size(MODEL.SAT,2); numUE = size(MODEL.UE,2);
% Norm 2 of the channel
CHAN = sqrt(abs(sum( conj(MODEL.Channel).*MODEL.Channel ,3)));
H_total = sum(sum(CHAN));
% Distance from SATs to users, row-SAT, column-UE
Distance = sqrt( (MODEL.SAT(1,:)'-MODEL.UE(1,:)).^2 + ...
                 (MODEL.SAT(2,:)'-MODEL.UE(2,:)).^2 + ...
                 (MODEL.SAT(3,:)'-MODEL.UE(3,:)).^2 );
T_total = sum(sum(Distance))/(3e8);
% Convert from t_thres to d_thres
d_thres = t_thres * 3e8;
% Find cannot-links
C_not = Distance > d_thres; % if C_not(m,u)=1 then m-u is cannot-link
IdClusUE_old = ones(1,numUE);
% Initialize randomly all clusters
IdClusUE = IdClusUE_RC;
NumUE_in_SAT = NumUE_in_SAT_RC;

%% Loop
numChange_old = numUE; numTimechange = 0; numLoop = 0;
while(~all(IdClusUE_old==IdClusUE))
    numLoop = numLoop + 1;
    IdClusUE_old = IdClusUE;
    for u = 1:numUE
        % mp is m prime
        mp = IdClusUE(u); nUEmp = NumUE_in_SAT(mp);
        vmp_old = func_utility(Distance,CHAN,mp,IdClusUE,theta,T_total,H_total);
        % The set of SATs which can connet to u
        IdClusUE_ok = IdClusUE; NumUE_in_SAT_ok = NumUE_in_SAT; gap = 0;
        M_ok_u = find(C_not(:,u)==0);
        % delete mp
        if(~isempty(find(M_ok_u==mp)))
            M_ok_u(find(M_ok_u==mp)) = [];
        end
        for im = 1:size(M_ok_u,2)
            IdClusUE_temp = IdClusUE; NumUE_in_SAT_temp = NumUE_in_SAT;
            m = M_ok_u(im); nUEm = NumUE_in_SAT_temp(m);
            vm_old = func_utility(Distance,CHAN,m,IdClusUE,theta,T_total,H_total);
            if(nUEm == B) % SWAP
                Um = find(IdClusUE==m);
                for iup = 1:nUEm
                    up = Um(iup);
                    if(C_not(mp,up))
                        continue
                    end
                    %Swap clusters of u and up
                    % The number does not change
                    IdClusUE_temp = IdClusUE;
                    IdClusUE_temp(u) = m; IdClusUE_temp(up) = mp;
                    vmp_new = func_utility(Distance,CHAN,mp,IdClusUE_temp,theta,T_total,H_total);
                    vm_new = func_utility(Distance,CHAN,m,IdClusUE_temp,theta,T_total,H_total);
                    gap_new = (vmp_old + vmp_old)-(vmp_new + vmp_new);
                    if(gap_new>gap)
                        IdClusUE_ok = IdClusUE_temp;
                        gap = gap_new;
                    end
                end
            else % JOIN
                IdClusUE_temp(u) = m; NumUE_in_SAT_temp(mp)=NumUE_in_SAT_temp(mp)-1;
                NumUE_in_SAT_temp(m)=NumUE_in_SAT_temp(m)+1;
                vmp_new = func_utility(Distance,CHAN,mp,IdClusUE_temp,theta,T_total,H_total);
                vm_new = func_utility(Distance,CHAN,m,IdClusUE_temp,theta,T_total,H_total);
                gap_new = (vmp_old + vmp_old)-(vmp_new + vmp_new);
                if(gap_new>gap)
                    IdClusUE_ok = IdClusUE_temp;
                    NumUE_in_SAT_ok = NumUE_in_SAT_temp;
                    gap = gap_new;
                end
            end
        end
        IdClusUE = IdClusUE_ok; NumUE_in_SAT = NumUE_in_SAT_ok;
    end
    % % if the change less than 10% then break
    Change = (IdClusUE_old==IdClusUE);
    numChange_new = size(find(Change==0),2);
%     fprintf("The number of change is %d \n",numChange_new);
    if((numChange_new-numChange_old)==0)
        numTimechange = numTimechange + 1;
    else
        numTimechange = 0;
    end
    if(numTimechange>2||numLoop==100)
        break
    end
    numChange_old = numChange_new;
end
MODEL.IdClusUE = IdClusUE;
MODEL.NumUE_in_SAT = NumUE_in_SAT;
end


%% Mini functions
function value = func_utility(Distance,CHAN,cluster,IdClusUE,theta,T_total,H_total)
% cluster is related to SAT m
IdUE = find(IdClusUE==cluster);
value = 0;
for iu = 1:size(IdUE,2)
    u = IdUE(iu);
    value = value + theta*Distance(cluster,u)/3e8/T_total...
            - (1-theta)*CHAN(cluster,u)/H_total;
end
end


