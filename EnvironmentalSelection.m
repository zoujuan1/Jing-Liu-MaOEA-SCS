function [Population,flag,threshold1,threshold2,flagth2] = EnvironmentalSelection(Population,W,N,D1,flag,threshold1,threshold2,flagth2)
% The environmental selection of MaOEA-SCS

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------


%     [QObj,z,znad] = Normalization(Population.objs,z,znad);
%     OBJ=Population.objs
%     Threshold=threshold1
    if flag==1 ;      
        St=(1:2*N+1);    

        [Sd1,tFrontNo] = tNDSort(Population.objs,W);


        MaxFNo    = find(cumsum(hist(tFrontNo,1:max(tFrontNo)))>=N,1);
        LastFront = find(tFrontNo==MaxFNo);
        

        AH=Sd1(LastFront);
        AQ=[AH',LastFront'];
        AW=sortrows(AQ);
        AE=flip(AW);
        AR=AE(:,2);
        LastFront=AR';
        tFrontNo(LastFront(1:sum(tFrontNo<=MaxFNo)-N)) = inf;
        Next      = St(tFrontNo<=MaxFNo);
        Population = Population(Next);
        [Sd2,~] = tNDSort(Population.objs,W);
        D2=mean(Sd2,2);
        if abs(D2-D1)<threshold1
            flag=0;
            flagth2=flagth2+1;
            threshold1=threshold1*1.4;
            if threshold1>1
                %%DTLZ
                threshold1=0.01;
                %%WFG
%                 threshold1=0.5;
            end
        end
        
%         TY=flagth2
            
%         OBJ2=Population.objs
    else
        
%         TT=flagth2
        
        %% Environmental selection2
        PopObj=ObjectiveNormalization(Population);        
        [tFrontNo,MaxFNo] = NDSort(PopObj,N);
        Next = tFrontNo <= MaxFNo;
        Population2 = Population(Next);
        PopObj2=ObjectiveNormalization(Population2);
        [~,FrontNo]  = tNDSort2(PopObj2,W); 
        MaxFNo =find(cumsum(hist(tFrontNo,1:max(FrontNo)))>=N,1);       
        Next = tFrontNo <= MaxFNo;             
        if sum(Next) < N
            Next = true;
        elseif sum(Next) > N
            Del  = Truncation(PopObj2,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [tFrontNo,MaxFNo] = NDSort(Population.objs,N);
%         Next = tFrontNo <= MaxFNo; 
%         if sum(Next) < N
%             Next = true;
%         elseif sum(Next) > N
%             Del  = Truncation(Population(Next).objs,sum(Next)-N);
%             Temp = find(Next);
%             Next(Temp(Del)) = false;
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Population for next generation
        Population = Population(Next);
        [Sd2,~] = tNDSort(Population.objs,W);
        D2=mean(Sd2,2);
%         ABS2=abs(D2-D1)
        
        if flagth2==1
            threshold2=max(0.5,3*abs(D2-D1));
%             TY=abs(D2-D1)
%             threshold2=1.5;
            flagth2=flagth2+1;
            
%             fprintf(1,'---------------------------------------\n');
        end

        if abs(D2-D1)>threshold2 
            flag=1;  
            threshold2=threshold2*0.95;
            if threshold2<0.05
                threshold2=0.1;
            end
        end
    end
end

function Del = Truncation(PopObj,K)
%     Distance = pdist2(PopObj,PopObj);
%     Distance(logical(eye(length(Distance)))) = inf;
%     Del = false(1,size(PopObj,1));
%     while sum(Del) < K
%         Remain   = find(~Del);
%         Temp     = sort(Distance(Remain,Remain),2);
%         [~,Rank] = sortrows(Temp);
%         Del(Remain(Rank(1))) = true;
%     end

% % Select part of the solutions by truncation   
    N = size(PopObj,1);    
    %% Calculate the shifted distance between each two solutions
    Distance = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    
    %% Truncation
    Del = false(1,N);
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
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
    tFrontNo = zeros(1,N);
    Sd1 = zeros(1,N);
    for i = 1 : NW 
        C = find(class==i);  
        [~,rank] = sort(d1(C,i));
        tFrontNo(C(rank)) = 1 : length(C);
        Sd1(C)=d1(C,i);
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

%     theta = zeros(1,NW) + 5;
%     theta(sum(W>1e-4,2)==1) = 1e6;
    
    tFrontNo = zeros(1,N);
    Sd1 = zeros(1,N);
    for i = 1 : NW 
        C = find(class==i);
        [~,rank] = sort(d2(C,i));
        tFrontNo(C(rank)) = 1 : length(C);
        Sd1(C)=d1(C,i);
    end  
end