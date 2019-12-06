function [fitness]=calc_fit2(nest,data,~,~)
    %input:
    %   nest        nPop*(n+m)   0/1 bits        
    %   lamda       vol weight
    %output:
    %   fitness   nPop*1    = resi + 1./vari;
    
        resi = calc_resi(nest,data);
    
        vari = calc_var(nest,data);
        fitness = resi+1./vari;
        
end