preamble; 
R = randn(4, 7) 
R([1, 2, 4], [3, 6, 7]) 
R(2:2:end, 1:2:end)
R(R(:, 1) > 0.5, :) 
R(:, R(end, :) < -0.1) 
S = R; 
S(S < 0.6) = -S(S < 0.6)
S(S < 0) = -inf
S(S > 1.5) = NaN
find(isinf(S) | isnan(S))'