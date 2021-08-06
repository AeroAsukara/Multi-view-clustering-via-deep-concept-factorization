function [W, V] = preSetKMeans(X, innerDim)

outerDim = size(X, 2);

labels = zeros(outerDim, 1);
labels = litekmeans(X', innerDim, 'Replicates', 20);
G = zeros(outerDim, innerDim);
for j = 1:innerDim
    G(:, j) = (labels(:) == j * ones(outerDim, 1));
end

Vt = G + 0.1 * ones(outerDim, innerDim);
Dw = diag(sum(G, 1))^-1;
W = (G + 0.1 * ones(outerDim, innerDim)) * Dw;
%[W, Vt] = NormalizerWV(Ker, W, Vt, NormV, Norm, innerDim);
V = Vt;


end

