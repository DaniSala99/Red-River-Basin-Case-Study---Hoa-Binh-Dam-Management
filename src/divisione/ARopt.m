function [stima,best] = ARopt(xc ,xv, mc, mv, sdc, sdv, qc, qv, p_max)

best=0;
R2_AR = nan*zeros(p_max,1);

for i = 1:p_max %per ogni colonna valuto ciascun ordine di esogeneità (i) per trovare l'ottimale
    
    %Calibrazione
    
    [teta_c, M_AR , y_cal] = stima_ARp(xc,i);
    % y_cal è la stima DESTAGIONALIZZATA del HoaBinh cumulative inflow 
    q_stima = y_cal.*sdc + mc;
    % q_stima è la stima RISTAGIONALIZZATA del HoaBinh cumulative inflow 
    
    %Validazione
    
    [M_AR_V, y_val] = validaz_ARp(xv, xv, teta_c);
    % y_val è la stima DESTAGIONALIZZATA del HoaBinh cumulative inflow (in
    % validaz)
    q_stima_v = y_val.*sdv + mv;
    % q_stima_v è la stima RISTAGIONALIZZATA del HoaBinh cumulative inflow
    % (in validaz)
    
    R2_AR(i) = 1 - sum((qv(i+1:end)-q_stima_v(i+1:end)).^2) / sum((qv(i+1:end)-mv(i+1:end)).^2);
    
    if R2_AR(i) > best
        best = R2_AR(i); 
        x = [xc;xv];
        sd = [sdc;sdv];
        m = [mc;mv];
        [M, y] = validaz_ARp(x,x,teta_c);
        stima = y.*sd + m;
    end
end
