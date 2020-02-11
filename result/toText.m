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

        % Genes = ones(bic_cnt,1);
        % Conditions = ones(bic_cnt,1);
        % GVS = ones(bic_cnt,1);
        % CVS = ones(bic_cnt,1);
        % MSRS = ones(bic_cnt,1);
        % VarS = ones(bic_cnt,1);

        csvfile = [algPath, '.txt'];
        fid = fopen(csvfile, 'w');
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

            % fprintf(fid, ['%s', ',', '%s', ',', '%f', ',', '%d', ',', '%d', ',', '%f','\n'],...
            % num2str(genes), num2str(conditions), scores(1), scores(3), scores(4),scores(5));
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