function r = std_operating_policy(h, policy)

% -- Get policy parameters --
c = [policy(1) policy(3) policy(5) policy(7)]; % dimensionality of a mean
b = [policy(2) policy(4) policy(6) policy(8)]; % dimensionality of a sd   
w_in = [policy(9) policy(10) policy(11) policy(12)];
w = w_in./sum(w_in);
k = policy(13);
% DOMANDARE SE SI PUO' INSERIRE DIRETTAMENTE UN VETTORE DI PESI

A = length(c); % number of radial basis functions

fi = nan*zeros(A,1);
u = nan*zeros(A,1);
in = h;
% -- Construct the policy using rbf
for i = 1:A
    fi(i) = exp(-sum((in - c(i))^2./(b(i)^2)));
    u(i) = w(i)*fi(i);
end
u_fin = sum(u);
r = k * u_fin;
% -- Verify release are non-negative  --
r( r < 0 )  = 0; 

end