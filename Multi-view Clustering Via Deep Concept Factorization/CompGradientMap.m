function newU = CompGradientMap( X, U_all, Phi, Vv, options)
%COMPGRADIENTMAP Summary of this function goes here
%   Detailed explanation goes here
% U是一个cell，存储了所有的base matrix
eta_d = 2;
eta_u = 2;
alpha = options.alphaT;
differror = options.error;
maxIter = options.submaxIter;
meanFitRatio = options.meanFitRatio;
Lt = options.L0;
debug = options.debug;
V = Vv';
layerNum = length(U_all);
U_M1 = X * Phi;%U_M-1

U = U_all{layerNum};%????U_M
[mFea,K]=size(U);

%Pre-compute XVt and VVt
UM1_XVt = U_M1'*X * V';
UM1_T_UM1 = U_M1'*U_M1;
VVt = V * V';
VXt = V * X';

objhistory = calcObjSub(X, U_M1, V, U, alpha, VVt, VXt);
if debug == 1
    disp(['  comp_grad_map: start optimization. Sub-objective is ',num2str(objhistory)]);
end
if isempty(maxIter)
    meanFit = objhistory*10;
end

maxErr = 1;
TLU = zeros(mFea,K); % used for storing temporary U^{t+1}
zeroconst = zeros(mFea,1);
t = 0;
while maxErr > differror 
    flag = 0;
    Grad = UM1_T_UM1*U*VVt - UM1_XVt;
    L = Lt;
    fobj_t = calcf(X, U_M1, V, U, VVt, VXt);%目标函数的f(U)部分
    gradnorm = sum(sum(Grad .* Grad));%梯度二范数的平方
    fobj = 0;
    while flag == 0
        % optimize m_L(U^t;U)
        B = U - (Grad/L);
        % First set TLU = 0 where B <= 0
        %%%TLU(B <= 0) = 0;
        % Then compute the part for which B > 0
%         if size(B,1) ~= size(TLU,1)
%             disp('bug');
%         end
        %disp(['size of B is ',num2str(size(B)),'size of TLU is ',num2str(size(TLU))]);
        if alpha == 0
            TLU = B;
        else
            for i = 1:K
                %disp(['here! i=',num2str(i)]);
                b = B(:,i);
                %disp(['size of b is ',num2str(size(b))]);
                posb = b(b>0);
                if (isempty(posb))
                    TLU  (:,i) = zeroconst;
                else
                    theta = computeTheta(posb,alpha/L);
                    TLU(:,i) = b - max(b - theta,0);
                end
                %disp(['here! i=',num2str(i)]);
            end
        end
        TLU(B < 0) = 0;
        
        % compare obj(TLU) and m_L(TLU)
        fobj = calcf(X, U_M1, V, TLU, VVt, VXt);
        %disp(fobj);
        tmpdiff = TLU - B;
        if (fobj - fobj_t - sum(sum(tmpdiff .* tmpdiff))*L*0.5 + gradnorm*0.5/L) <= 0%衡量L是否合适，使辅助函数满足条件
            flag = 1;
        else
            L = L*eta_u;
        end
    end
    U = TLU;
    Lt = max(options.L0, L /eta_d);
    
    newobj = fobj + alpha * sum(max(U));
    objhistory = [objhistory newobj]; %#ok<AGROW>
    
    t = t+1;
    if debug == 1
        if newobj - objhistory(end-1) <= 0
            disp(['  comp_grad_map: ',num2str(t),'-th iteration completed. Sub-objective is ',num2str(newobj),'. Diff with last iteration: ',num2str(newobj - objhistory(end-1))]);
        else
            warning(['  comp_grad_map: ',num2str(t),'-th iteration completed. Sub-objective is ',num2str(newobj),'. Diff with last iteration: ',num2str(newobj - objhistory(end-1))]);
        end
    end
    
    if isempty(maxIter)
        meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
        maxErr = (meanFit-newobj)/meanFit;
    else
        maxErr = 1;
        if t >= maxIter
            maxErr = 0;
        end
    end
end

newU = U;

