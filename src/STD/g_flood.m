function g_flo = g_flood(h_HN)

h_F = 950; % flood threshold [cm]

g_flo = (h_HN - h_F).^2;
g_flo(h_HN <= h_F) = 0;
    
end
