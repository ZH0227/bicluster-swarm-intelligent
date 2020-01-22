nPop=100;
runtime = 1;
saveIf =0;
for data = 1:1
    if data == 1
        dataName = 'BCLL';
        lamda=500;
        miu=10;
        omega=0.03;
    elseif data ==2
        dataName = 'YC';
        lamda=1000000;
        miu=1000;  
        omega=300000;
    elseif data==3
        dataName = 'RAT';
        lamda=20000000;
        miu=10000;
        omega=30000000;
    elseif data==4 
        dataName = 'PBC';
        lamda=2000;
        miu=30;
        omega=0.1;
    end
    for i = 1:1

        [bicC,costC,scoreC,historyC] = csb(nPop,data,lamda,miu,omega);
        [bicF,costF,scoreF,historyF] = fab(nPop,data,lamda,miu,omega);
        [bicPSO,costPSO,scorePSO,historyPSO] = psob(nPop,data,lamda,miu,omega);
        [bicQPSO,costQPSO,scoreQPSO,historyQPSO] = qpsob(nPop,data,lamda,miu,omega);
        [bicCF,costCF,scoreCF,historyCF] = csfa(nPop,data,lamda,miu,omega);
        figure;

        plot(historyC,'LineWidth',2);hold on;
        plot(historyF,'LineWidth',2);hold on;
        plot(historyPSO,'LineWidth',2);hold on;
        plot(historyQPSO,'LineWidth',2);hold on;
        plot(historyCF,'LineWidth',2);
        legend('CS','FA','PSO','QPSO','CSFA')
        set(get(gca, 'XLabel'), 'String', 'Number of iterations');
        set(get(gca, 'YLabel'), 'String', 'Fitness value');
        if saveIf
            PngFile = ['./exper2/',dataName,num2str(i),'.png']
            print(gcf,'-dpng',PngFile)

            scoreFile = ['./exper2/',dataName,num2str(i),'_score.mat']
            score = [scoreC;scoreF;scoreCF];
            save(scoreFile, 'score')

            bicFile = ['./exper2/',dataName,num2str(i),'_bic.mat']
            bic = [bicC;bicF;bicCF];
            save(bicFile, 'bic')

            costFile = ['./exper/',dataName,num2str(i),'_cost.mat']
            cost = [costC;costF;costCF];
            save(costFile, 'cost')
        end
    end
    % close(figure(gcf));
end
% print -deps Le