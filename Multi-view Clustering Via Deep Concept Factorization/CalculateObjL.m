
function [L] = CalculateObjL(data, W, V, L,Vcon, options, view_num, layers)

tempNMF = zeros(1,view_num);
tempLap = zeros(1,view_num);
tempVVc = zeros(1,view_num);

for i = 1:view_num

    XWV = data{i};
    
    for j = 1 : numel(layers)
        
        XWV = XWV * W{i, j} * V{i, j}';
        
    end
    
    for j = 1 : numel(layers)
        
        VL = V{i, j}' * L{i};
        VLV = VL * V{i, j};
        tempLap(i) = tempLap(i) + options.alpha(j) * trace(VLV);
    
    end
    
    tempNMF(i) = tempNMF(i) + norm(data{i}-XWV,'fro')^2;
    
    VtVc = V{i, numel(layers)} - Vcon;
    
    tempVVc(i) = tempVVc(i) + options.beta * norm(VtVc,'fro')^2;
    
    %tempWn(i) = tempWn(i) + options.eta * sum(max(W{i, numel(layers)}));
    
    clear KW KWV VL VLV;
end
obj_NMF = sum(tempNMF);
obj_Lap = sum(tempLap);
obj_VVc = sum(tempVVc);

L = obj_NMF + obj_Lap + obj_VVc;
%obj_Wn = sum(tempWn);
end



