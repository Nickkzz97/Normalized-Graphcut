function [Inc2 Inc Knc] = N_weightmatrix(I,SI,SX,r,sNcut,sArea,sb)

%% ncutImageSegment
[nRow, nCol,c] = size(I);                  % size of the image
N = nRow * nCol;                           %making th iage into 1 column so that it can be represented as a graph
V = reshape(I, N, c);                      %  Vertices of Graph
if c==3
    m=1:3;
else
    m=1;
end
%% ncutComputeW  
% Step 1. Compute weight matrix W, and D
W = sparse(N,N);  
% if c == 3  %% define feature vector
% F = F2(I);
% else
% F = F1(I);
%end% size of weight matrix
F = reshape(I, N, 1, c);                   % col vector % Spatial Location
X = cat(3, repmat((1:nRow)', 1, nCol), repmat((1:nCol), nRow, 1));
X = reshape(X, N, 1, 2);                   % col vector
%% calculation of weight matrix
for ic=1:nCol                              % Future Work: Reduce computation to half. It can be done because W is symmetric mat
    for ir=1:nRow                          % matlab tricks for fast computation (Avoid 'for' loops as much as possible, instead use repmat.)
        
        % This range satisfies |X(i) - X(j)| <= r (block distance)
        jc = (ic - floor(r)) : (ic + floor(r)); % vector
        jr = ((ir - floor(r)) :(ir + floor(r)))';
        jc = jc(jc >= 1 & jc <= nCol);
        jr = jr(jr >= 1 & jr <= nRow);
        jN = length(jc) * length(jr);

        % index at vertex. V(i)
        i = ir + (ic - 1) * nRow;
        j = repmat(jr, 1, length(jc)) + repmat((jc -1) * nRow, length(jr), 1);
        j = reshape(j, length(jc) * length(jr), 1); % a col vector

        % spatial location distance (disimilarity)
        XJ = X(j, 1, :);
        XI = repmat(X(i, 1, :), length(j), 1);
        DX = XI - XJ;
        DX = sum(DX .* DX, 3); % squared euclid distance

        % |X(i) - X(j)| <= r (already satisfied if block distance measurement)
        constraint = find(sqrt(DX) <= r);
        j = j(constraint);
        DX = DX(constraint);

        % feature vector disimilarity
        FJ = F(j, 1, :);
        FI = repmat(F(i, 1, :), length(j), 1);
        DF = FI - FJ;
        DF = sum(DF .* DF, 3); % squared euclid distance ( DF = sum(abs(DF), 3); % block distance)
        W(i, j) = exp(-DF / (SI*SI)) .* exp(-DX / (SX*SX)); %weight matrix
    end
end

%% ncutPartition
Seg = (1:N)';                             % Step 5. recursively repartition
%seg=Seg;
id = 'ROOT';                              % the first segment has whole nodes. [1 2 3 ... N]'
% Compute D
N = length(W);
d = sum(W, 2);
D = spdiags(d, 0, N, N); % diagonal matrix
% Step 2 and 3. Solve generalized eigensystem (D -W)*S = S*D*U (12).
warning off; %  stop warning
[U,S] = eigs(D-W, D, 2, 'sm');
% 2nd smallest (1st smallest has all same value elements, and useless)
U2 = U(:, 2);
% Bipartition the graph at point that Ncut is minimized.
t = mean(U2);
ncut = NcutValue(t, U2, W, D);
t = fminsearch('NcutValue', t, [], U2, W, D);
A = find(U2 > t);
B = find(U2 <= t);

% Step 4. Decide if the current partition should be divided

x = (U2 > t);
x = (2 * x) - 1;
d = diag(D);
k = sum(d(x > 0)) / sum(d);
b = k / (1 - k);
y = (1 + x) - b * (1 - x);
ncut = (y' * (D - W) * y) / ( y' * D * y );
if (length(A) < sArea || length(B) < sArea) || ncut > sNcut
     Seg  = {Seg};
     %Seg(1)=seg;
     Id(1)   = {id};   % for debugging
     Ncut(1) = {ncut}; % for duebugging
     return
 end
no_iteration=0;
% Seg segments of A
[SegA  IdA  NcutA no_iteration] = NcutPartition(Seg(A), W(A, A), sNcut, sArea, [id '-A'],no_iteration,sb);
% Seg segments of B
%no_iteration=2;
[SegB  IdB  NcutB  no_iteration] = NcutPartition(Seg(B), W(B, B), sNcut, sArea, [id '-B'],no_iteration,sb);
% concatenate cell arrays
Seg  = [SegA SegB];
Id   = [IdA IdB];
Ncut = [NcutA NcutB];

%% show
Inc  = zeros(size(I),'uint8');
for k=1:length(Seg)
 [r, c] = ind2sub(size(I),Seg{k});
 for i=1:length(r)
 Inc(r(i),c(i),1:3) = uint8(round(mean(V(Seg{k}, :))));
 end
end
Knc = length(Seg);
%figure(2), imshow(Inc);  title(['NormalizedCut',' : ',num2str(Knc)]);

for k=1:length(Seg)
 Inc2  = zeros(size(I),'uint8');
 [r, c] = ind2sub(size(I),Seg{k}); %gives indices of k
 for i=1:length(r)
 Inc2(r(i),c(i),m) = I(r(i),c(i),m);
 end
 %figure(),
 subplot(1,Knc,k); imshow(Inc2), title(['segmentation',num2str(k)]) 
end
Knc2 = length(Seg);
disp(Id)
disp(Ncut)
end

% needs:
% NcutPartition
% NcutValue
% F1 - F4: Compute a feature vector F
% F = F1(I) % for point sets
% F = F2(I) % intensity
% F = F3(I) % hsv, for color
% F = F4(I) % DOOG
%
function F = F1(I)
F = (I == 0);
end
function F = F2(I)
% intensity, for gray scale
F = I;
end
function F = F3(I)
F = I; % raw RGB converting to hsv and find values future work
end
%Doog filter future work
%I have done color segemntation using RGB value but it can also be done by
%designing feature vector using HSV values as mentioned in the Normalised
%segmentation paper by Shimallik.