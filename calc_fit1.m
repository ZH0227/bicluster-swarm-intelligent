function [fitness]=calc_fit1(nest,data,lamda,~)
    %input:
    %   nest        nPop*(n+m)   0/1 bits        
    %   lamda       vol weight
    %output:
    %   fitness   nPop*1   = resi+lamda./volMat;
        n = size(data,1);
    
        resi = calc_resi(nest,data);
    
        vol = cumVol(nest,n);
        volMat = vol(:,1);
        invVolMat = lamda./volMat;

        fitness = resi+invVolMat;
        
end

