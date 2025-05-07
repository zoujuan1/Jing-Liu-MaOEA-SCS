function MatingPool = MatingSelection(PopObj,W,flag)

% The mating selection of MaOEA-CSS

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------


    [N,~] = size(PopObj);
    if flag==1
        [ASF,Front]  = tNDSort(PopObj,W); 
    else
        [ASF,Front]  = tNDSort2(PopObj,W); 
    end

    MatingPool = zeros(1,N);
    for i = 1 : N
        p = randperm(N,2);
        if Front(p(1)) < Front(p(2))
            p = p(1);
        elseif Front(p(1)) > Front(p(2))
            p = p(2);
        else
            if ASF(p(1)) < ASF(p(2))
                p = p(1); 
            else
                p = p(2);
            end
        end
        MatingPool(i) = p;
    end
end

function [Sd1,tFrontNo] = tNDSort(PopObj,W)
    N  = size(PopObj,1);
    NW = size(W,1);

    normP  = sqrt(sum(PopObj.^2,2));
    Cosine = 1 - pdist2(PopObj,W,'cosine');
    d1     = repmat(normP,1,size(W,1)).*Cosine;
    d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2); 
    

    [~,class] = min(d2,[],2);
    
    %% Sort
    tFrontNo = zeros(1,N);
    Sd1 = zeros(1,N);
    for i = 1 : NW 
        C = find(class==i);
        [~,rank] = sort(d1(C,i));
        tFrontNo(C(rank)) = 1 : length(C);
        Sd1(C)=d2(C,i);
    end  
end


function [Sd1,tFrontNo] = tNDSort2(PopObj,W) 
    N  = size(PopObj,1);
    NW = size(W,1);

    normP  = sqrt(sum(PopObj.^2,2));
    Cosine = 1 - pdist2(PopObj,W,'cosine');
    d1     = repmat(normP,1,size(W,1))./Cosine;  
    d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2); 
    
    [~,class] = min(d2,[],2);
    
    tFrontNo = zeros(1,N);
    Sd1 = zeros(1,N);
    for i = 1 : NW
        C = find(class==i);
        [~,rank] = sort(d2(C,i));
        tFrontNo(C(rank)) = 1 : length(C);
        Sd1(C)=d2(C,i);
    end  
end