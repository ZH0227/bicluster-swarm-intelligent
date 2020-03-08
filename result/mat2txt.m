rootPath = pwd;
datas = {'BCLL', 'PBC', 'RAT', 'YC'};
data_N_map = containers.Map;
data_N_map('BCLL') = 12185;
data_N_map('PBC') = 21225;
data_N_map('RAT') = 7751;
data_N_map('YC') = 5847;

algs = {'CS', 'FA', 'CSFA', 'PSO', 'QPSO'};

for i=1:length(datas)
    data = char(datas(i));
    N = data_N_map(data);
    for j=1:length(algs)
        alg = char(algs(j));
        algPath = char(fullfile(rootPath, data, alg));
        dirOut = dir(char(fullfile(algPath, '*.mat')));
        files = {dirOut.name};
        bic_cnt = fix(length(files) / 4);


        txtFile = [algPath, '.txt'];
        disp(txtFile)
        fid = fopen(txtFile, 'w');
        for k=1:bic_cnt
            load(fullfile(algPath,[num2str(k, '%04d'), '_', 'bic.mat']));
            load(fullfile(algPath,[num2str(k, '%04d'), '_', 'score.mat']));
            bic = eval(['bic', alg]);
            scores = eval(['score', alg]);
            [~, genes] = find(bic(1:N)==1);
            [~, conditions] = find(bic(N+1:length(bic))==1);
            %matlab index is start from 1, but python is 0
            genes = genes -1;
            conditions = conditions -1;

            fprintf(fid, '%f\n',scores(1));
            fprintf(fid, '%d\n',scores(3));
            fprintf(fid, '%d\n',scores(4));
            fprintf(fid, '%f\n',scores(5));
            fprintf(fid, '%s\n',num2str(genes));
            fprintf(fid, '%s\n',num2str(conditions));
        end
        fclose(fid);
    end
end