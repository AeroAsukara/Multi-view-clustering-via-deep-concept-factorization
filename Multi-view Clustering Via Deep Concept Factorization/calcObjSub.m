

function obj = calcObjSub(X, U_M1, V, U, alpha, VVt, VXt)
obj = calcf(X, U_M1, V, U, VVt, VXt);
obj = obj + alpha * sum(max(U));
