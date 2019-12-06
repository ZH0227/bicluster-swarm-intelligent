function data = choseData(data)
    if data == 1
        disp('import data: leukemia_big.csv')
        file = importdata('../data/leukemia_big.csv');
        data = file.data; %size 7128x72
    
    elseif data==2
        disp('import data: yeast.matrix.txt')
        data = importdata('../data/yeast.matrix.txt');%size 2884x17
                
    elseif data == 3
        disp('import data: DLBCL.txt')
        data = importdata('../data/DLBCL.txt'); %size  12625x21
        
    elseif data==4
        disp('import data: breast_Matlab.txt')
        data = importdata('../data/breast_Matlab.txt'); %size 13666x117

    elseif data == 5
        disp('import data: ratStrain.txt')
        data = importdata('../data/ratStrain.txt'); %size  8799x122        
    end

    
end