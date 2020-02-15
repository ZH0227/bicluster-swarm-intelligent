function result = calc_bench(pop)
%pop    nPop*(data_n+data_m)   0/1 
%result nPop*5        residue,bicV,geneV,condV,var

    resi = calc_resi(pop);
    Vol = cumVol(pop);
    vari = calc_var(pop);

    result = [resi Vol vari];
end