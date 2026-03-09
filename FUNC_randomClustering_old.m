function [MODEL] = FUNC_randomClustering_old(MODEL,B,t_thres,theta,Area)
numSAT = size(MODEL.SAT,2); numUE = size(MODEL.UE,2);
% Distance from SATs to users, row-SAT, column-UE
Distance = sqrt( (MODEL.SAT(1,:)'-MODEL.UE(1,:)).^2 + ...
                 (MODEL.SAT(2,:)'-MODEL.UE(2,:)).^2 + ...
                 (MODEL.SAT(3,:)'-MODEL.UE(3,:)).^2 );
% Convert from t_thres to d_thres
d_thres = t_thres * 3e8;
% Find cannot-links
C_not = Distance > d_thres; % if C_not(m,u)=1 then m-u is cannot-link
IdClusUE = zeros(1,numUE);
% Initialize randomly all clusters
NumUE_in_SAT = zeros(1,numSAT);
numLoop = 0;
while(1)
    numLoop = numLoop + 1;
    U_free = find(IdClusUE==0);
    if(isempty(U_free))
        break
    end
    for iu = 1:size(U_free,2)
        u = U_free(iu);
        % Consider UEs that have not been served yet
        if(IdClusUE(u)==0)
            % Find a free near SAT and the number of UEs conneting less
            % than B
            NearSAT = find(Distance(:,u)<d_thres);
            Free_SAT = find(NumUE_in_SAT < B);
            [OkSAT,pos]=intersect(NearSAT,Free_SAT);
            OkSAT = OkSAT';
            if(isempty(OkSAT)) % worst case
                % Connect to the nearest SAT in the worst case
                [dismin,ind_SAT] = min(Distance(:,u));
                % Delete a UE in the nearest SAT
                U_in_indSAT = find(IdClusUE==ind_SAT);
                IdClusUE(U_in_indSAT(1)) = 0;
                % The u-th UE is served by (ind_SAT)-th SAT
                IdClusUE(u) = ind_SAT;
            else % random case
                ind_SAT = OkSAT(randi([1,size(OkSAT,2)],1));
                % The number UEs associate with (ind_SAT)-th SAT
                NumUE_in_SAT(ind_SAT) = NumUE_in_SAT(ind_SAT) + 1;
                % The u-th UE is served by (ind_SAT)-th SAT
                IdClusUE(u) = ind_SAT;
            end
        end
    end
    if(numLoop == 1000)
        error("Cannot clustering due to strict constraints");
    end
end
% % show random clustering
% MODEL.IdClusUE = IdClusUE;
% MODEL.NumUE_in_SAT = NumUE_in_SAT;
% FUNC_Plot_model_2D(MODEL, Area, "cluster")

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
