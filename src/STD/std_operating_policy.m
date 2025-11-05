function r = std_operating_policy(h, policy)

% -- Get policy parameters --
h1 = policy(1);
h2 = policy(2);
m1 = policy(3);
m2 = policy(4);
w = policy(5);

% -- Construct the policy using piecewise linear functions --
% water saving
L1 = w + m1 * ( h - h1 ); 
% flood control
L2 = w + m2 * ( h - h2 );
% release
r  = max( [ min( L1 , w ) ; L2 ] );

% -- Verify release are non-negative  --
r( r < 0 )  = 0; 

end