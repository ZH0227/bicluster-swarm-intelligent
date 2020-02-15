function [fitness]=calc_fit4(nest,lamda,miu,omega)
    %input:
    %   nest        nPop*(data_n+data_m)   0/1 bits        
    %   lamda       volGene weight
    %   miu         VolCon weight
    %   omega       Var weight
    %output:
    %   fitness   nPop*1   =resi + lamda./volGene + miu./volCond + omega./vari;
        resi = calc_resi(nest);
    
        vol = cumVol(nest);
        volGene = vol(:,2);
        volCond = vol(:,3);
        invGene = lamda./volGene;
        invCond = miu./volCond;
        % vari = calc_var(nest);
        % invVar = omega./vari;
        
        fitness = resi+invGene+invCond;
    
end
