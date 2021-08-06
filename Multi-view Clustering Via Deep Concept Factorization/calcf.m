function o = calcf(X, U_M1, V, U, VVt, VXt)
% compute 1/2( Tr[UVV^{T}U^{T}] - 2Tr[UVX^{T}] ) which is more efficient
% than computing the least square
% 
%auxi = U_M1*U*V;
a1 = U_M1 * U * VVt * U' * U_M1';
a2 = U_M1 * U * VXt;
o =  trace(a1) - trace(a2);