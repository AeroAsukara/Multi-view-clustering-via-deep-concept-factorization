function [Vall_, obj_] = multiViewClusteringVarDeepConceptFactorization(data, options, layers)


differror = options.error;
maxIter = options.maxIter;
minIter = options.minIter;
Rounds = options.round;
alpha = options.alpha;
beta = options.beta;
gamma = options.gamma;
pi = options.pi;
V_all = zeros();
obj = zeros();
Norm = 2;
NormV = 1;
layerNum = numel(layers);
viewNum = numel(data);

for i = 1:viewNum
    
    tempDim = size(data{i}, 2);
    Ker{i} = constructKernel(data{i}', data{i}', options);
    
    TempWt = constructW(data{i}', options);
    if isfield(options, 'NormWeight') && strcmpi(options.NormWeight, 'NCW')
        
        D_mhalf = sum(TempWt, 2).^-.5;
        D_mhalf = spdiags(D_mhalf, 0, tempDim, tempDim);
        TempWt = D_mhalf * TempWt * D_mhalf;
        
    end
    Wt{i} = TempWt;
    
    DCol = full(sum(Wt{i}, 2));
    D{i} = spdiags(DCol, 0, speye(size(Wt{i}, 1)));
    TempL = D{i} - Wt{i};
    if isfield(options, 'NormLap') && options.NormLap
        D_mhalf = DCol.^-.5;
        tempD_mhalf = repmat(D_mhalf, 1, nSmp);
        TempL = (tmpD_mhalf .* TempL) .* tmpD_mhalf';
        clear D_mhalf tmpD_mhalf;
        TempL = max(TempL, TempL');
    end
    L{i} = TempL;
    
    clear Kt TempWt DCol TempL;
    
end

iR = 0;
while iR < Rounds
    
    iR = iR + 1;
    W = cell(viewNum, layerNum);
    V = cell(viewNum, layerNum);
    Vcon = zeros();
    % -------------------------------------------------------------------init W,V-------------------------------------------------------------
    for v_ind = 1:viewNum
        
        X = data{v_ind};
        for i_layer = 1:numel(layers)
            
            if i_layer == 1
                H = X;
            else
                H = V{v_ind, i_layer - 1}';
            end
            
            [W{v_ind, i_layer}, V{v_ind, i_layer}] = preSetCF(H, layers(i_layer), options);
            %[W{v_ind, i_layer}, V{v_ind, i_layer}] = preSetKMeans(H, layers(i_layer));
            [W{v_ind, i_layer}, V{v_ind, i_layer}] = ...
                NormalizeWV(Ker{v_ind}, W{v_ind, i_layer}, V{v_ind, i_layer}, NormV, Norm, layers(i_layer));
            
            
        end
        Vcon = Vcon + (pi(v_ind)) * V{v_ind, numel(layers)};
        clear Kt Dw G tlabel
    end
    
       
    iter = 0;
    selectInit = 1;
    count = 0;
    while (iter < maxIter)
        
        iter = iter + 1;
        VC = zeros();
        DVcV = zeros();
        for v_ind = 1:viewNum
            
            V_tilde{v_ind, numel(layers)} = eye([size(data{v_ind}, 2), size(data{v_ind}, 2)]);
            for i_layer = numel(layers) - 1 : -1 : 1
                
                V_tilde{v_ind, i_layer} = W{v_ind, i_layer + 1} * V{v_ind, i_layer + 1}' * V_tilde{v_ind, i_layer + 1};
                
            end
            for i = 1 : numel(layers)
                
                if i == 1
                    numerator = Ker{v_ind} * V_tilde{v_ind, i}' * V{v_ind, i};
                    denominator = Ker{v_ind} * W{v_ind, i} * V{v_ind, i}' * (V_tilde{v_ind, i} * V_tilde{v_ind, i}') * V{v_ind, i};
                    W{v_ind, i} = W{v_ind, i} .* (numerator ./ max(denominator, 1e-9));
                    
                    numerator = V_tilde{v_ind, i} * Ker{v_ind} * W{v_ind, i} + alpha(i) * Wt{v_ind} * V{v_ind, i};
                    denominator = ...
                        (V_tilde{v_ind, i} * V_tilde{v_ind, i}') * V{v_ind, i} * W{v_ind, i}' * Ker{v_ind} * W{v_ind, i} + alpha(i) * D{v_ind} * V{v_ind, i};
                    V{v_ind, i} = V{v_ind, i} .* (numerator ./ max(denominator, 1e-9));
                else
                    numerator = Phi' * Ker{v_ind} * V_tilde{v_ind, i}' * V{v_ind, i};
                    denominator = Phi' * Ker{v_ind} * Phi * W{v_ind, i} * V{v_ind, i}' * (V_tilde{v_ind, i} * V_tilde{v_ind, i}') * V{v_ind, i};
                    W{v_ind, i} = W{v_ind, i} .* (numerator ./ max(denominator, 1e-9));
                end
                
                if i == numel(layers)
                    numerator = Ker{v_ind} * Phi * W{v_ind, i} + alpha(i)  * Wt{v_ind} * V{v_ind, i} + beta * (pi(v_ind)) * Vcon;
                    denominator = V{v_ind, i} * W{v_ind, i}' * Phi' * Ker{v_ind} * Phi * W{v_ind, i} + alpha(i) * D{v_ind} * V{v_ind, i} + beta * (pi(v_ind)) * V{v_ind, i}; 
                    V{v_ind, i} = V{v_ind, i} .* (numerator ./ max(denominator, 1e-9));
                   
                    VC = VC + (pi(v_ind)) * V{v_ind, i};
                else
                    if i ~= 1
                        numerator = V_tilde{v_ind, i} * Ker{v_ind} * Phi * W{v_ind, i} + alpha(i) * Wt{v_ind} * V{v_ind, i};
                        denominator = (V_tilde{v_ind, i} * V_tilde{v_ind, i}') * V{v_ind, i} * W{v_ind, i}' * Phi' * Ker{v_ind} * Phi * W{v_ind, i} + alpha(i) * D{v_ind} * V{v_ind, i};
                        V{v_ind, i} = V{v_ind, i} .* (numerator ./ max(denominator, 1e-9));
                    end
                end
                
                [W{v_ind, i}, V{v_ind, i}] = NormalizeWV(Ker{v_ind}, W{v_ind, i}, V{v_ind, i}, NormV, Norm, layers(i));
               
                if i == 1
                    
                    Phi = W{v_ind, i} * V{v_ind, i}';
                else
                    Phi = Phi * W{v_ind, i} * V{v_ind, i}';
                    
                end
               
            end
                   
        end
        Vcon = VC;
        %disp(size(Vcon));
%         for i = 1 : viewNum
%             
%             DVcV(i) = norm(V{i, numel(layers)} - Vcon, 'fro')^2;
%             dnorm_w(i) = (gamma * (DVcV(i)))^(1 / (1 - gamma));
%             
%         end
%         
%         for v_ind = 1 : viewNum
%             pi(v_ind) = dnorm_w(v_ind) / sum(dnorm_w);
%         end
        
        if (selectInit)
            [obj_NMF, obj_Lap, obj_VVc] = CalculateObj(data, W, V, L,Vcon, options, pi, viewNum, layers);
            new_obj = obj_NMF + obj_Lap + obj_VVc;
            objhistory = new_obj;
            selectInit = 0;
        else
            [obj_NMF, obj_Lap, obj_VVc] = CalculateObj(data, W, V, L,Vcon, options, pi, viewNum, layers);
            new_obj = obj_NMF + obj_Lap + obj_VVc;
        end
        if iter == 1
            Vall = Vcon;
            Veach = V(:, numel(layers));
        end
        if new_obj < objhistory(end)
            
            Vall = Vcon;
            Veach = V(:, numel(layers));
            objhistory = [objhistory new_obj];
        end
        
        if iter == 1
            new_error = 1;
        else
            new_error = abs(new_obj-objhistory(end-1));
        end
            
        
        
        if (new_error < differror)
            count = count + 1;
            if (count == 5)
                iter = maxIter;
            end
        end
        clear new_obj;
    
    end

    if iR == 1
        Vall_ = Vall;
        Veach_ = Veach;
        newobj_ = objhistory(end);
        obj_ = objhistory;
        %newobj_
    end
%     if objhistory(end) < newobj_
%         Vall_ = Vall;
%         Veach_ = Veach;
%         newobj_ = objhistory(end);
%         obj_ = objhistory;
%         %newobj_
%         Valll = Vall;
%         printResult(Valll, truelabel{1}, numC, options.clusteringFlag);
%     end
    if iR < Rounds
        fprintf('restart...\n');
        %clear objhistory;
    else
        fprintf('get the result:\n');
    end


end

