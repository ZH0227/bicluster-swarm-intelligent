function [fitness]=calc_fit(nest,data,lamda,miu)
    %input:
    %   nest        nPop*(n+m)   0/1 bits        
    %   lamda       volGene weight
    %   miu         VolCon weight
    %output:
    %   fitness   nPop*1   =resi + lamda./volGene + miu./volCond;
        n = size(data,1);
    
        resi = calc_resi(nest,data);
    
        vol = cumVol(nest,n);
        volGene = vol(:,2);
        volCond = vol(:,3);
        invGene = lamda./volGene;
        invCond = miu./volCond;

        fitness = resi+invGene+invCond;
    
        
end

