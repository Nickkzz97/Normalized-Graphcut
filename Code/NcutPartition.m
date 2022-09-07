function [Seg Id Ncut no_iteration] = NcutPartition(I, W, sNcut, sArea, id,no_iteration,sb)
% NcutPartition - Partitioning
%
% Synopsis
%  [sub ids ncuts] = NcutPartition(I, W, sNcut, sArea, [id])
%
% Description
%  Partitioning. This function is called recursively.
%
% Inputs ([]s are optional)
%  (vector) I        N x 1 vector representing a segment to be partitioned.
%                    Each element has a node index of V (global segment).
%  (matrix) W        N x N matrix representing the computed similarity
%                    (weight) matr
%                    W(i,j) is similarity between node i and j.
%  (scalar) sNcut    The smallest Ncut value (threshold) to keep partitioning.
%  (scalar) sArea    The smallest size of area (threshold) to be accepted
%                    as a segment.
%  (string) [id]     A label of the segment (for debugg)
%
% Outputs ([]s are optional)
%  (cell)   Seg      A cell array of segments partitioned.
%                    Each cell is the each segment.
%  (cell)   Id       A cell array of strings representing labels of each segment.
%                    IDs are generated as children based on a parent id.
%  (cell)   Ncut     A cell array of scalars representing Ncut values
%                    of each segment.
%
% Requirements
%  NcutValue

N = length(W);
d = sum(W, 2);
D = spdiags(d, 0, N, N);    % diagonal matrix
% Step 2 and 3. Solve generalized eigensystem (D -W)*S = S*D*U (12).
warning off;                % let me stop warning
[U, S] = eigs(D-W, D, 2, 'sm');
% 2nd smallest (1st smallest has all same value elements, and useless)
U2 = U(:, 2);
% Bipartition the graph at point that Ncut is minimized.
t = mean(U2);
t = fminsearch('NcutValue', t, [], U2, W, D);
A = find(U2 > t);       %%Partiotioning
B = find(U2 <= t);
if isreal(U2)
H=histogram(U2);        %stability condition
Hmin=min(H.Values);
Hmax=max(H.Values);
Hratio=Hmin/Hmax;
else 
    Hratio=0;
end

ncut = NcutValue(t, U2, W, D); %ncut value for checking minumum cut

% Step 5. Decide if the current partition should be divided
% if either of partition is too small, stop recursion.
% if Ncut is larger than threshold, stop recursion.



if (length(A) < sArea || length(B) < sArea) || (ncut > sNcut) || (no_iteration > sb-3) || (Hratio > 0.06)
    Seg(1)  = {I};
    Id(1)   = {id};          % for debugging
    Ncut(1) = {ncut};        % for duebuggin
    
    return;
end
no_iteration=no_iteration+1;
% Seg segments of A
[SegA IdA  NcutA   no_iteration] = NcutPartition(I(A), W(A, A), sNcut, sArea, [id '-A ncut'],no_iteration,sb);
% Seg segments of B
[SegB IdB NcutB no_iteration] = NcutPartition(I(B), W(B, B), sNcut, sArea, [id '-B ncut'],no_iteration,sb);
% concatenate cell arrays
Seg   = [SegA SegB];
Id   = [IdA IdB];
Ncut = [NcutA NcutB];
end
