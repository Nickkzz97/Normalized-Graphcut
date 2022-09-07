function ncut = NcutValue(t, U2, W, D)
% NcutValue - 2.1 Computing the Optimal Partition Ncut. eq (3)
%
%  ncut = NcutValue(T, U2, D, W);
%
% Inputs ([]s are optional)
%  (scalar) t        splitting point (threshold)
%  (vector) U2       N x 1 vector representing the 2nd smallest
%                     eigenvector computed at step 2.
%  (matrix) W        N x N weight matrix
%  (matrix) D        N x N diagonal matrix
%
% Outputs ([]s are optional)
%  (scalar) ncut     The value calculated at the right term of eq (10).
%                    This is used to find minimum Ncut.


%% NcutValue - Computing the Optimal Partition Ncut
x = (U2 > t);
x = (2 * x) - 1;  
d = diag(D);
k = sum(d(x > 0)) / sum(d);
b = k / (1 - k);
y = (1 + x) - b * (1 - x);
ncut = (y' * (D - W) * y) / ( y' * D * y );
end