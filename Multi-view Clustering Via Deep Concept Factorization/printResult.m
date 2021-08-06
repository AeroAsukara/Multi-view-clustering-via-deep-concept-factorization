function [ac,nmi_value,f_score,ari, indic] = printResult(X, label, K, clusteringFlag)
Xt = X;
if clusteringFlag == 1
    indic = litekmeans(Xt, K, 'Replicates',20);
else
    [~, indic] = max(Xt, [] ,2);
end
[ac, nmi_value, f_score, ari, cnt] = CalcMeasures(label, indic); % including bestMap()
fprintf('ACC:%0.4f\t%d/%d\tNMI:%0.4f\tFScore:%0.4f\tARI:%0.4f\t\n', ac, cnt, length(label), nmi_value,f_score, ari);