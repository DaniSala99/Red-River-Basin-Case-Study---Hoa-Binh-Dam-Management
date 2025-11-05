function g_hyd = g_hydropower( r, h )

% set parameters
k = 86400/(3.6*10^6)*1000; % numerical coefficient
nu_n = 0.85; % nominal efficiency
qmax_T = 2360; % turbine capacity [m3/s]
qmin_T = 38; % minimum flow for turbine operation [m3/s]

z_t = 0.000000000003663570691434010*r.^3 - 0.000000136377708978325000000*r.^2 ...
    + 0.002087770897833130000000000*r + 11.660165118679700000000000000;

dh_t = h - z_t;

nu_t = nu_n.*(-0.0007476023419322190*dh_t.^2 + 0.1370692867959150000*dh_t ...
    + 3.0285414473892400000);

qt_T = nan.*r;
for i = 1:length(r)
    if r(i) > qmax_T
        qt_T(i) = qmax_T;
    elseif r(i) <= qmin_T
        qt_T(i) = 0; 
    else
        qt_T(i) = r(i);
    end     
end

% compute objective
g_hyd = k * nu_t .* qt_T .* dh_t;

end