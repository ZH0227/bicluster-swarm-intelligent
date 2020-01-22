% 算法主程序开始 Start FA
function [bic,cost,score,history]=fab(nPop,data,lamda,miu,omega)
%   data    n*m
%   lamda   the weight of gene Volume
%   miu     the weight of condition Volume
    %% Problem Definiton
    data = choseData(data);
    n = size(data,1);
    m = size(data,2);
    d = n+m;        % Number of Unknown (Decision) Variables
    % The Flag for Showing Iteration Information
    ShowIterInfo = 1;
    costFun = @calc_fit4;

    Lb = zeros(1,d);
    Ub = 1.*ones(1,d);
    c2bT = 0.5;         % decide 0 or 1 threshold
    % n=number of fireflies
    % MaxGeneration=number of pseudo time steps
    % ------------------------------------------------
    % alpha=0.25;      % Randomness 0--1 (highly random)
    % betamn=0.20;     % minimum value of beta
    % gamma=1;         % Absorption coefficient
    % ------------------------------------------------
    n=nPop;
    MaxGeneration=300;
    early_stopping_cnt = 0;
    early_stopping_maxcnt = 30;
    early_stopping_threshold = 1.0*10^-4;
    alpha=0.25;
    betamin=0.20;
    gamma=1;

    history = zeros(MaxGeneration,1);
    % 初始化萤火虫位置 generating the initial locations of n fireflies
    [ns,zn] = init_ffa(n,d,Lb,Ub);
    % 对每个萤火虫计算目标函数值 Evaluate new solutions (for all n fireflies)
    for i=1:n
        bic = conti2bit(ns(i,:),c2bT);
        zn(i)=costFun(bic,data,lamda,miu,omega);
    end
    LightbestO = min(zn);
    for k=1:MaxGeneration % 迭代开始
        % 更新alpha（可选）This line of reducing alpha is optional
        alpha=alpha_new(alpha,MaxGeneration);

        % 根据亮度排序 Ranking fireflies by their light intensity/objectives
        [Lightn,Index]=sort(zn);
        ns_tmp=ns;
        for i=1:n
            ns(i,:)=ns_tmp(Index(i),:);
        end

        %% 找出当前最优 Find the current best
        nso=ns;         %Old
        Lighto=Lightn;  %Old
        nbest=ns(1,:);
        Lightbest=Lightn(1);

        % early stopping
        change = LightbestO - Lightbest;
        if change < early_stopping_threshold
            early_stopping_cnt = early_stopping_cnt + 1;
        else
            early_stopping_cnt = 0;
        end

        % 向较优方向移动 Move all fireflies to the better locations
        [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,alpha,betamin,gamma,Lb,Ub);
        
        % 对每个萤火虫计算目标函数值 Evaluate new solutions (for all n fireflies)
        for i=1:n
            bic = conti2bit(ns(i,:),c2bT);
            zn(i)=costFun(bic,data,lamda,miu,omega);
        end
        
        if ShowIterInfo
            disp(['Iteration ' num2str(k) ': Best Cost = ' num2str(Lightbest)]);
        end
        history(k) = Lightbest;
        LightbestO = Lightbest;

        if early_stopping_cnt > early_stopping_maxcnt
            disp('fab early stoping');
            break
        end
    end %end of iterations

    if k < MaxGeneration
        for i =k:MaxGeneration
            history(i) = Lightbest;
        end
    end
    bic = conti2bit(nbest,c2bT);
    cost = Lightbest;
    score = calc_bench(bic,data);
    
    % costs = zeros(n,1);
    % bics = zeros(n,d);
    % for i = 1:n
    %     bics(i,:) = conti2bit(ns(i,:),c2bT);
    %     costs(i) =  costFun(bics(i,:),data,lamda,miu,omega);
    % end
    % score = calc_bench(bics,data);
end
% ----- All the subfunctions are listed here ------------
% 初始化萤火虫位置 The initial locations of n fireflies
function [ns,Lightn]=init_ffa(n,d,Lb,Ub)
    ns=zeros(n,d);
    for i=1:n
        ns(i,:)=Lb+(Ub-Lb).*rand(1,d); % 则在取值范围内随机取值
    end

    % 初始化目标函数 initial value before function evaluations
    Lightn=ones(n,1)*10^100;
end
% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,alpha,betamin,gamma,Lb,Ub)
% 参数取值范围绝对值 Scaling of the system
    scale=abs(Ub-Lb);

    % 更新萤火虫 Updating fireflies
    for i=1:n
        % The attractiveness parameter beta=exp(-gamma*r)
        for j=1:n
            r=sqrt(sum((ns(i,:)-ns(j,:)).^2));
            % Update moves
            if Lightn(i)>Lighto(j) % 如果i比j亮度更强 Brighter and more attractive
                beta0=1;
                beta=(beta0-betamin)*exp(-gamma*r.^2)+betamin;
                tmpf=alpha.*(rand(1,d)-0.5).*scale;
                ns(i,:)=ns(i,:).*(1-beta)+nso(j,:).*beta+tmpf;
            end
        end % end for j
    end % end for i

    % 防止越界 Check if the updated solutions/locations are within limits
    [ns]=findlimits(n,ns,Lb,Ub);
end
% This function is optional, as it is not in the original FA
% The idea to reduce randomness is to increase the convergence,
% however, if you reduce randomness too quickly, then premature
% convergence can occur. So use with care.
% alpha参数更新函数 
function alpha=alpha_new(alpha,NGen)
    % alpha_n=alpha_0(1-delta)^NGen=10^(-4);
    % alpha_0=0.9
    delta=1-(10^(-4)/0.9)^(1/NGen);
    alpha=(1-delta)*alpha;
end
% 防止越界 Make sure the fireflies are within the bounds/limits
function [ns]=findlimits(n,ns,Lb,Ub)
    for i=1:n
        % Apply the lower bound
        ns_tmp=ns(i,:);
        I=ns_tmp<Lb;
        ns_tmp(I)=Lb(I);
        % Apply the upper bounds
        J=ns_tmp>Ub;
        ns_tmp(J)=Ub(J);
        % Update this new move
        ns(i,:)=ns_tmp;
    end
end