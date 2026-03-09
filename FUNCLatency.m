function [tk, case_tk] = FUNCLatency(MODEL, k, f, A, B, P_sat, P_uav_uav, P_uav_ue,...
              Area, Q, Noise_var, BW, ASU, AUU)
% n_request: number packages requested simultaneously
numUAV = size(MODEL.UAV, 2);
xk = MODEL.UE(1,k); yk = MODEL.UE(2,k);
numUAVeachSAT = size(MODEL.UAV,2)/size(MODEL.SAT,2);
numUEeachUAV = size(MODEL.UE,2)/size(MODEL.UAV,2);
s = ceil(k/(numUAVeachSAT*numUEeachUAV));
xs = MODEL.SAT(1,s); ys = MODEL.SAT(2,s); zs = MODEL.SAT(3,s);
Us = numUAV/4*(s-1)+1:numUAV/4*s;
u = find(A(:,k)==1);

xu = MODEL.UAV(1,u); yu = MODEL.UAV(2,u); zu = MODEL.UAV(3,u);
d_uk = sqrt((xu-xk)^2 + (yu-yk)^2 +(zu-0)^2);
h_uk = MODEL.H_UK(u,k);
p_uk = P_uav_ue(u,k);
H_other = MODEL.H_UK(u,:); H_other(k) = [];
P_other = P_uav_ue(u,:); P_other(k) = [];

% for u = 1:numUAV
% choose one out of three scenario
if(B(u,f))
    % Case 1: GU k is served by primary UAV
    t_pri = funcLatencyUK(d_uk, h_uk, p_uk, H_other, P_other, Q, Noise_var, BW);
    tk = B(u,f)*t_pri;
    case_tk = 1;
    % Case 2: GU k is served by a secondary UAV
elseif((B(u, f)==0)&&any(B(Us, f) == 1))
    Up = numUAV/4*(s-1) + find(B(Us,f)==1); % set of UAVs up has beta(u,f)=1
    H_norm_Up_u = abs(MODEL.H_UU(Up, u));
    [h_max, id_umax] = max(H_norm_Up_u);
    up = Up(id_umax);
    n_request = AUU(up, u);
    %t_sec = 1e6;
    %for idup = 1:size(Up)
%         up = Up(idup);
        xup = MODEL.UAV(1,up); yup = MODEL.UAV(2,up); zup = MODEL.UAV(3,up);
        d_upu = sqrt((xup-xu)^2 + (yup-yu)^2 +(zup-zu)^2);
        h_upu = MODEL.H_UU(up,u);
        p_upu = P_uav_uav(up, u);
        t_sec_temp = funcLatencyUpK(d_uk, h_uk, p_uk, H_other, P_other, ...
                d_upu, h_upu, p_upu, Q, Noise_var, BW, n_request);
%         if((t_sec_temp~=0) && (t_sec_temp<t_sec))
%             t_sec = t_sec_temp;
%         end
%     end
%     tk = (1 - B(u,f))*t_sec;
        tk = (1 - B(u,f))*t_sec_temp;
        case_tk = 2;

    % Case 3: GU k is served by satellite s
else
    d_su = sqrt((xs-xu)^2 + (ys-yu)^2 +(zs-zu)^2);
    h_su = MODEL.H_SU(s,u,:); precoding_su = MODEL.PrecodingS(s,u,:);
    p_su = P_sat(s,u);
    n_request = ASU(s, u);
    t_sat = funcLatencySK(d_uk, h_uk, p_uk, H_other, P_other, ...
        d_su, h_su, precoding_su, p_su, Q, Noise_var, BW, n_request);
    tk = (1 - B(u,f))*(1 - max(B(Us,f)))*t_sat;
    case_tk = 3;
end
if((tk==Inf)||(isnan(tk)))
    tk = 1;
end
% end
end

% function s = funcFindSforK(Area, xk, yk)
% x_area = Area(1); y_area = Area(2);
% if(xk < x_area/2)
%     if(yk < y_area/2)
%         s = 1;
%     else
%         s = 4;
%     end
% else
%     if(yk < y_area/2)
%         s = 2;
%     else
%         s = 3;
%     end
% end
% end