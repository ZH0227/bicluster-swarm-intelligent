function data = choseData(data)
    if data == 1
        disp('import data: gse2403_BCLL.csv')
        file = importdata('../data/gse2403_BCLL.csv');
        data = file.data; %size 18125x21
    
    elseif data==2
        disp('import data: gds2350_YC.csv')
        file = importdata('../data/gds2350_YC.csv');
        data = file.data; %size 5847x50
                
    elseif data == 3
        disp('import data: gse952_RAT.csv')
        file = importdata('../data/gse952_RAT.csv');
        data = file.data; %size 7751x122
        
    elseif data==4
        disp('import data: gse2034_PBC.csv')
        file = importdata('../data/gse2034_PBC.csv');
        data = file.data; %size 21225x286
    
    end
end