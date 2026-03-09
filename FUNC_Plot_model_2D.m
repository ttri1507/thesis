function FUNC_Plot_model_2D(MODEL, Area, mode)
SAT = MODEL.SAT; UAV = MODEL.UAV; UE = MODEL.UE;
if(mode=="normal")
    x_area = Area(1); y_area = Area(2);
    % plot the considered area
    plot([0 x_area x_area 0 0],[0 0 y_area y_area 0],'k','markersize',8,'Linewidth',4);
    hold on;
    for m = 1:size(SAT,2)
        plot(SAT(1,m),SAT(2,m),'r-^','markersize',10,'Linewidth',4);
    end
    for us = 1:size(UAV,2)
        plot(UAV(1,us),UAV(2,us),'b^','markersize',9,'Linewidth',2);
    end
    plot(UE(1,:),UE(2,:),'k*','markersize',8,'Linewidth',1.5);
    grid on
elseif(mode=="cluster")
    x_area = Area(1); y_area = Area(2);
    % plot the considered area
    plot([0 x_area x_area 0 0],[0 0 y_area y_area 0],'k','markersize',8,'Linewidth',4);
    hold on;
    plot(UE(1,:),UE(2,:),'b^','markersize',8,'Linewidth',1.5);
    IdClusUE = MODEL.IdClusUE;
    NumUE_in_SAT = MODEL.NumUE_in_SAT;
    for u = 1:size(UE,2)
        m = IdClusUE(u);
        plot([UE(1,u) SAT(1,m)], [UE(2,u) SAT(2,m)],'k','Linewidth',2);
    end
    for m = 1:size(SAT,2)
        plot(SAT(1,m),SAT(2,m),'r-^','markersize',10,'Linewidth',4);
    end
    grid on
end
end