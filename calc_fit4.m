function [fitness]=calc_fit4(nest,data,lamda,miu,omega)
    %input:
    %   nest        nPop*(n+m)   0/1 bits        
    %   lamda       volGene weight
    %   miu         VolCon weight
    %   omega       Var weight
    %output:
    %   fitness   nPop*1   =resi + lamda./volGene + miu./volCond + omega./vari;
        n = size(data,1);
    
        resi = calc_resi(nest,data);
    
        vol = cumVol(nest,n);
        volGene = vol(:,2);
        volCond = vol(:,3);
        invGene = lamda./volGene;
        invCond = miu./volCond;
        vari = calc_var(nest,data);
        invVar = omega./vari;
        
        fitness = resi+invGene+invCond+invVar;
    
end
