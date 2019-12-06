function [nest,fitness,score]=csb1(nPop,data,lamda,miu)
%input:
%   nPop           pop number
%   data        gene expression data n by m
%   p           found and abandon Probability
%   threshold   
%   initRow,initCol
%output:
%   nest        nPop*(n+m)   0/1 bits
%   socre       nPop*5 (resi,bic,gene and sample Volume,var)

    if data == 1
        disp('import data: leukemia_big.csv')
        file = importdata('../data/leukemia_big.csv');
        data = file.data; %size 7128x72
        threshold = 0.05;
        initRow = 200;
        initCol = 10;    
    elseif data==2
        disp('import data: yeast.matrix.txt')
        data = importdata('../data/yeast.matrix.txt');
        threshold = 300;    %size 2884x17
        initRow = 200;
        initCol = 8;
    elseif data==3
        disp('import data: breast_Matlab.txt')
        data = importdata('../data/breast_Matlab.txt'); %size 13666x117
        threshold = 0.12;
        initRow = 200;
        initCol = 15;
    elseif data == 4
        disp('import data: ratStrain.txt')
        data = importdata('../data/ratStrain.txt'); %size  8799x122
        threshold = 0.12;
        initRow = 200;
        initCol = 15;
    elseif data == 5
        disp('import data: DLBCL.txt')
        data = importdata('../data/DLBCL.txt'); %size  12625x21
        threshold = 1000;
        initRow = 200;
        initCol = 15;
    end
    ShowIterInfo = 1; 
    p = 0.25;
    
    iter = 1000;
    
    n = size(data,1);
    m = size(data,2);
    %initial nest
    nest1 = zeros(nPop, n);
    nest2 = zeros(nPop, m);
    nest1 = ranFlip(nest1,initRow);
    nest2 = ranFlip(nest2,initCol);
    nest = [nest1, nest2];
    fitness = calc_fit(nest,data,lamda,miu);
    %iter begin
    for i = 1:iter

        %cuckoo search
        [nest,fitness] = cuckoo_search(nest,fitness,data,lamda,miu);
        %abandon
        [nest,fitness] = abandoning(nest,fitness,p,data,lamda,miu);
        resi = calc_resi(nest,data);

        if ShowIterInfo ==1
            disp(['iter==>',num2str(i)])
            disp(['min fitness: ', num2str(min(fitness))] )
            % disp(['max fitness: ', num2str(max(fitness))] )
            disp(['min resi: ', num2str(min(resi))] )

        end
        if min(fitness)<threshold
            break;
        end

        % if min(resi)<threshold
        %     break;
        % end
    end
    resi = calc_resi(nest,data);
    disp(['iter==>',num2str(i)])
    disp(['min resi: ', num2str(min(resi))] )
    disp(['min fitness: ', num2str(min(fitness))] )
    vol = cumVol(nest,n);
    disp(['vol:', num2str(vol(nPop,:))])
    disp(['mean fitness: ', num2str(mean(fitness))] )
    disp(['meanvol:', num2str(mean(vol))])
    vari = calc_var(nest,data);
    score = [resi,vol,vari];
end

function [nest] = ranFlip(nest, mSize)
    for i = 1: size(nest,1)
        ran = randperm(size(nest,2),mSize);
        for j = 1:mSize
            nest(i,ran(j)) = ~nest(i,ran(j));
        end
    end
end

function [nest,fitness] = cuckoo_search(nest,fitness,data,lamda,miu)
    mSize = 2;
    maxNumber =10;
    %keep the best bicluster,not change
    [bestV,bestIn] = min(fitness);
    best = nest(bestIn,:);
    nest(bestIn,:)=[];
    fitness(bestIn)=[];
    for i=1:size(nest,1)
        nestIN = ranFlip(nest(i,:),mSize);
        fitN = calc_fit(nestIN,data,lamda,miu);
        if fitN < fitness(i) %better
            nest(i,:) = nestIN;
        else
            n =1;
            for j=1:maxNumber
                nestIN = ranFlip(nest(i,:),mSize);
                fitN = calc_fit(nestIN,data,lamda,miu);
                if fitN < fitness(i)%better
                    nest(i,:) = nestIN;
                    break;
                end   
                n = n+1;                         
            end
            if n >10
                nestIN = ranFlip(nest(i,:),8);
                nest(i,:) = nestIN;
            end
        end
    end
    %keep the best bicluster,not change
    nest = [nest;best];
    fitness = [fitness;bestV];
end

function [nest,fitness] = abandoning(nest,fitness,p,data,lamda,miu)
%input:
%   nest        nPop*(n+m)   0/1 bits        
%   fitness        nPop*1
%   p           1*1
%   data        n by m
%output:
%   nest        nPop*(n+m)   0/1 bits        
%   fitness        nPop*1
    num = fix(size(nest,1)*p)+1;
    [~,ind] = sort(fitness);
    ind = ind(1:num);

    abandon = nest(ind,:);
    fitnessAb = calc_fit(abandon,data,lamda,miu);
    [abandonN,fitnessAbN] = cuckoo_search(abandon,fitnessAb,data,lamda,miu);

    nest(ind,:) = [];
    fitness(ind,:) = [];

    nest = [nest;abandonN];
    fitness = [fitness;fitnessAbN];

end
