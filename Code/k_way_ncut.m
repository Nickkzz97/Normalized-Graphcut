%function [Inc2 Inc Knc] = K_way_cut(I,SI,SX,r,sNcut,sArea,sb)
image_orginal=imread('C:\Users\Nicky Nirlipta\Desktop\Nicky Nirlipta Sahoo_EE19S042-ISP PROJECT\SEGMENTED OUTPUT\color 26-50\189080.JPG');
IM = image_orginal;
I=imresize(IM,0.5);
SI=7;
SX   = 8;                       % Spatial similarity

r    =1.5;                        % Spatial threshold (less than r pixels apart)
sNcut = 0.21;                   % The smallest Ncut value (threshold) to keep partitioning
sArea = 120;%
sb=4;
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
[U,S] = eigs(D-W, D, 5, 'sm');
U2 = U(:, 2);
t2 = mean(U2);
ncut = NcutValue(t2, U2, W, D);
t2 = fminsearch('NcutValue', t2, [], U2, W, D);
A1 = find(U2 > t2);
B1 = find(U2 <= t2);
no_iteration=0;
% Seg segments of A
[SegA1  IdA1  NcutA1 no_iteration] = NcutPartition(Seg(A1), W(A1, A1), sNcut, sArea, [id '-A1'],no_iteration,sb);
% Seg segments of B
%no_iteration=2;
[SegB1  IdB1  NcutB1  no_iteration] = NcutPartition(Seg(B1), W(B1, B1), sNcut, sArea, [id '-B1'],no_iteration,sb);
% concatenate cell arrays

U3 = U(:, 3);
t3 = mean(U3);
ncut = NcutValue(t3, U3, W, D);
t3 = fminsearch('NcutValue', t3, [], U3, W, D);
A2 = find(U3 > t3);
B2 = find(U3 <= t3);
no_iteration=0;
% Seg segments of A
[SegA2  IdA2  NcutA2 no_iteration] = NcutPartition(Seg(A2), W(A2, A2), sNcut, sArea, [id '-A2'],no_iteration,sb);
% Seg segments of B
%no_iteration=2;
[SegB2  IdB2  NcutB2  no_iteration] = NcutPartition(Seg(B2), W(B2, B2), sNcut, sArea, [id '-B2'],no_iteration,sb);
% concatenate cell arrays

U4 = U(:, 4);
t4 = mean(U4);
ncut = NcutValue(t4, U4, W, D);
t4 = fminsearch('NcutValue', t4, [], U4, W, D);
A3 = find(U4 > t4);
B3 = find(U4 <= t4);
no_iteration=0;
% Seg segments of A
[SegA3  IdA3  NcutA3 no_iteration] = NcutPartition(Seg(A3), W(A3, A3), sNcut, sArea, [id '-A3'],no_iteration,sb);
% Seg segments of B
%no_iteration=2;
[SegB3  IdB3  NcutB3  no_iteration] = NcutPartition(Seg(B3), W(B3, B3), sNcut, sArea, [id '-B3'],no_iteration,sb);
% concatenate cell arrays

U5 = U(:, 5);
t5 = mean(U5);
ncut = NcutValue(t5, U5, W, D);
t5 = fminsearch('NcutValue', t5, [], U5, W, D);
A4 = find(U5 > t5);
B4 = find(U5 <= t5);
no_iteration=0;
% Seg segments of A
[SegA4  IdA4  NcutA4 no_iteration] = NcutPartition(Seg(A4), W(A4, A4), sNcut, sArea, [id '-A4'],no_iteration,sb);
% Seg segments of B
%no_iteration=2;
[SegB4  IdB4  NcutB4  no_iteration] = NcutPartition(Seg(B4), W(B4, B4), sNcut, sArea, [id '-B4'],no_iteration,sb);
% concatenate cell arrays
Seg  = [SegA1 SegB1 SegA2 SegB2 SegA3 SegB3 SegA4 SegB4 ];
Id   = [IdA1 IdB1 IdA2 IdB2 IdA3 IdB3 IdA4 IdB4 ];
Ncut = [NcutA1 NcutB1 NcutA2 NcutB2 NcutA3 NcutB3 NcutA4 NcutB4];

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
%end

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



