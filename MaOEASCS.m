function MaOEASCS(Global)
% <algorithm> <H-N>
%A Many-objective Algorithm based on Staged Coordination Selection
%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [W,Global.N] = UniformPoint(Global.N,Global.M);

    Population = Global.Initialization();
    flag=1;


    %%DTLZ
    threshold1=0.005;
    threshold2=10086;
    flagth2=0;

    %%WFG
    %threshold1=0.5;
    %threshold2=1.5;


    while Global.NotTermination(Population)
        
        MatingPool = MatingSelection(Population.objs,W,flag);      
        %MatingPool = randi(Global.N,1,Global.N);
        Offspring  = Global.Variation(Population(MatingPool));
        [Sd1,~] = tNDSort(Population.objs,W);
        D1=mean(Sd1,2);
        [Population,flag,threshold1,threshold2,flagth2] = EnvironmentalSelection([Population,Offspring],W,Global.N,D1,flag,threshold1,threshold2,flagth2);
    end
end

function [Sd1,tFrontNo] = tNDSort(PopObj,W)
    N  = size(PopObj,1);
    NW = size(W,1);
    normP  = sqrt(sum(PopObj.^2,2));
    Cosine = 1 - pdist2(PopObj,W,'cosine');
    d1     = repmat(normP,1,size(W,1)).*Cosine ;
    d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2);   
    [~,class] = min(d2,[],2);
    tFrontNo = zeros(1,N);%初始一个1*N的0矩阵
    Sd1 = zeros(1,N);%初始一个1*N的0矩阵
    for i = 1 : NW %从每个向量上找出同一级的个体
        C = find(class==i);
        [~,rank] = sort(d1(C,i));
        tFrontNo(C(rank)) = 1 : length(C);
        Sd1(C)=d1(C,i);
    end  
end