clear;

addpath(genpath('Dataset/'));
addpath(genpath('Results/'));
dataDirectory = 'Dataset/';
resultDirectory = 'Results/';

if(~exist('Results', 'file'))
    mkdir('Results');
    addpath(genpath('Results/'));
end

dataName = {'3Sources', 'BBC', 'BBCSport', '20NewsGroups', 'Yale', 'HW2sources'};
numberOfDatasets = length(dataName);

for objectDataset = 1 : 1
    
    options = [];
    options.WeightMode = 'HeatKernel';
    options.tt = 10;
    options.NormWeight = 'NCW';
    options.k = 5;
    options.KernelType = 'Gaussian';
    options.maxIter = 200;
    options.minIter = 20;
    options.round = 1;
    options.repeat = 10;
    options.error = 1e-2;
    options.clusteringFlag = 1;
    options.beta = 100;
    options.gamma = 0.5;
    options.pi = zeros();
    options.PiFlag = 1;
    options.alpha = [1 10];
    layers = [100 50];
    
    odata = objectDataset;
    formatData = [dataDirectory, cell2mat(dataName(odata))];
    load(formatData);
    fprintf('Performing on dataset: %s\n', cell2mat(dataName(odata)));
    
    %normalize dataset
    numberOfView = numel(data);
    
    for i = 1 : numberOfView
        
        data{i} = NormalizeData(data{i}, 2);
        options.pi(i) = 1 / numberOfView; 
%         options.pi(i) = 1;
        %qqqqqqqqqqqqqqqqqqqsssssssddddfffffggggghhhhhjjjjjkkkkkkllll;;;;;'''///.
        
    end
    
    numberOfCluster = length(unique(truelabel{1}));
    numberOfFeature = size(data{1}, 2);
    fprintf('Dataset feature: %d\n', numberOfFeature);
    
    ACC = zeros();
    NMI = zeros();
    FScore = zeros();
    ARI = zeros();
    
    Vcon = cell(1, options.repeat);
    %resultLabel = cell(1, options.repeat);
    objectiveFunctionValue = cell(1, options.repeat);
    
    for iter = 1 : options.repeat
        
        fprintf('Perform %d-th times of %d...', iter, options.repeat);
        
        [Vcon{iter}, objectiveFunctionValue{iter}] = ...
        multiViewClusteringVarDeepConceptFactorization(data, options, layers);
        [ACC(iter), NMI(iter), FScore(iter), ARI(iter), ~] = ...
            printResult(Vcon{iter}, truelabel{1}, numberOfCluster, options.clusteringFlag);
        
    end
    
    Result = zeros(10, options.repeat);
    
    Result(1,:) = ACC;
    Result(2,:) = NMI;
    Result(3,:) = FScore;
    Result(4,:) = ARI;
    Result(5,1) = mean(ACC);
    Result(5,2) = mean(NMI);
    Result(5,3) = mean(FScore);
    Result(5,4) = mean(ARI);
    Result(6,1) = std(ACC);
    Result(6,2) = std(NMI);
    Result(6,3) = std(FScore);
    Result(6,4) = std(ARI);
    
    save([resultDirectory,char(dataName(odata)),'resultLit.mat'],'Result','Vcon','objectiveFunctionValue');
    
    %plot(objectiveFunctionValue{1}(1,1:end),'LineWidth',4,'Color','r'); title 'HW2Sources'; xlabel('Iteration'); ylabel('Objective Function Value');
    %fprintf('mean(ACC):%0.4f\n',mean(ACC));
    fprintf('mean(NMI):%0.4f\n',mean(NMI));
    fprintf('mean(FScore):%0.4f\n',mean(FScore));
    fprintf('mean(ARI):%0.4f\n',mean(ARI));
    
end









































