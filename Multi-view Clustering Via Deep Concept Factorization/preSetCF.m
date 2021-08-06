function [W, V] = preSetCF(X, innerDim, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
outerDim = size(X, 2);

%K = X' * X;
K = constructKernel(X', X', options);
%KG = gpuArray(K);
V = rand(outerDim, innerDim);
W = rand(outerDim, innerDim);
%VG = gpuArray(V);
%WG = gpuArray(W);
%XG = gpuArray(X);

for i = 1 : 80
    
    numerator = K * V;
    denominator = K * W * (V' * V);
    W = W .* (numerator ./ max(denominator, 1e-9)); 
    numerator = K * W;
    denominator = V * W' * K * W;
    V = V .* (numerator ./ max(denominator, 1e-9)); 
    
    
end
%W = gather(WG);
%V = gather(VG);

end

