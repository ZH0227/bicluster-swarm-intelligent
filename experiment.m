nPop=10;
runtime = 1;
saveIf =1;
for data = 3:4
    if data == 1
        dataName = 'Le';
        lamda=500;
        miu=10;
        omega=0.03;
    elseif data ==2
        dataName = 'Ye';
        lamda=1000000;
        miu=1000;  
        omega=300000;
    elseif data==3
        dataName = 'DL';
        lamda=20000000;
        miu=10000;
        omega=30000000;
    elseif data==4 
        dataName = 'Br';
        lamda=2000;
        miu=30;
        omega=0.1;
    end
    for i = 1:2

        [bicC,costC,scoreC,historyC] = csb(nPop,data,lamda,miu,omega);
        [bicF,costF,scoreF,historyF] = fab(nPop,data,lamda,miu,omega);
        [bicCF,costCF,scoreCF,historyCF] = csfa(nPop,data,lamda,miu,omega);
        figure;

        plot(historyC,'LineWidth',2);hold on;
        plot(historyF,'LineWidth',2);hold on;
        plot(historyCF,'LineWidth',2);
        legend('CS','FA','CSFA')
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