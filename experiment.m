nPop=100;
runtime = 1;
saveIf =1;
global data c2bT Lb Ub data_n dim
c2bT = 0.5;         % decide 0 or 1 threshold
for no_data = 2:2
    if no_data == 1
        dataName = 'BCLL';
        lamda=1*10^4;
        miu=1*10^3;
        omega=0;
    elseif no_data ==2
        dataName = 'YC';
        lamda=2*10^5;
        miu=2*10^3;
        omega=0;
    elseif no_data==3
        dataName = 'RAT';
        lamda=6*10^5;
        miu=1*10^4;
        omega=0;
    elseif no_data==4 
        dataName = 'PBC';
        lamda=3*10^6;
        miu=4*10^4;
        omega=0;
    end
    data = choseData(no_data);
    data_n = size(data ,1);
    data_m = size(data ,2);
    dim = data_n + data_m;
    Lb = zeros(1,dim);
    Ub = 1.*ones(1,dim);
    for i = 2:2

        [bicCS,costCS,scoreCS,historyCS] = csb(nPop,lamda,miu,omega);
        [bicFA,costFA,scoreFA,historyFA] = fab(nPop,lamda,miu,omega);
        [bicPSO,costPSO,scorePSO,historyPSO] = psob(nPop,lamda,miu,omega);
        [bicQPSO,costQPSO,scoreQPSO,historyQPSO] = qpsob(nPop,lamda,miu,omega);
        [bicCSFA,costCSFA,scoreCSFA,historyCSFA] = csfa(nPop,lamda,miu,omega);

        % figure;
        % plot(historyCS,'LineWidth',2);hold on;
        % plot(historyFA,'LineWidth',2);hold on;
        % plot(historyPSO,'LineWidth',2);hold on;
        % plot(historyQPSO,'LineWidth',2);hold on;
        % plot(historyCSFA,'LineWidth',2);
        % legend('CS','FA','PSO','QPSO','CSFA')
        % set(get(gca, 'XLabel'), 'String', 'Number of iterations');
        % set(get(gca, 'YLabel'), 'String', 'Fitness value');

        if saveIf
            data_root = fullfile('./result', dataName);
            algs = ["FA","CS","CSFA","PSO","QPSO"];
            for k=1:size(algs,2)
                alg = char(algs(k));
                alg_path = fullfile(data_root, alg);
                if exist(alg_path, 'dir')==0
                    mkdir(alg_path);
                end
                records = ["bic","history","score","cost"];
                for j=1:size(records,2)
                    rec = char(records(j));
                    %./result/dataName/algName/*_*.mat
                    file_name = [alg_path,'/',num2str(i,'%04d'),'_',rec,'.mat'];
                    save(file_name, [rec,alg]);
                end
            end
            % PngFile = fullfile(data_root,)
            % print(gcf,'-dpng',PngFile)

        end %end of save

    end %end of times
    % close(figure(gcf));
end % end of dataset