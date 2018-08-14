function [ beta1vec, beta2vec ] = IntersectMajorants( b1_a, b2_a, b1_b, b2_b, n_interp )

%% find lower bound for beta1
b1_min = max( b1_a(1), b1_b(1) );
b1_max = max( b1_a(end), b1_b(end) );

%% interpolate data on common (fine) grid
if( ~exist('n_interp', 'var') )
    n = 5*max( length(b1_a), length(b1_b) );
else
    n = n_interp;
end
% interpolation of b2 from vectors a and b

beta1vec = linspace(b1_min, b1_max, n);
b2_a_int = interp1( b1_a, b2_a, beta1vec );
b2_b_int = interp1( b1_b, b2_b, beta1vec );
beta2vec = max( b2_a_int, b2_b_int );
