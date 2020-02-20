% 算法主程序开始 Start FA
function [bic,cost,score,history]=fab(nPop,lamda,miu,omega)
%   data    n*m
%   lamda   the weight of gene Volume
%   miu     the weight of condition Volume
    %% Problem Definiton
    % The Flag for Showing Iteration Information
    ShowIterInfo = 1;
    global costFun

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
    early_stopping_maxcnt = 40;
    alpha=0.25;
    betamin=0.20;
    gamma=1;

    history = zeros(MaxGeneration,1);
    % 初始化萤火虫位置 generating the initial locations of n fireflies
    [ns,~] = init_ffa(n);
    % 对每个萤火虫计算目标函数值 Evaluate new solutions (for all n fireflies)
    bic = conti2bit(ns);
    zn = costFun(bic,lamda,miu,omega);
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
        if change < LightbestO/10000
            early_stopping_cnt = early_stopping_cnt + 1;
        else
            early_stopping_cnt = 0;
        end

        % 向较优方向移动 Move all fireflies to the better locations
        [ns]=ffa_move(n,ns,Lightn,nso,Lighto,alpha,betamin,gamma);
        
        % 对每个萤火虫计算目标函数值 Evaluate new solutions (for all n fireflies)
        bic = conti2bit(ns);
        zn=costFun(bic,lamda,miu,omega);
        
        if ShowIterInfo
            disp(['fab=> ','Iteration ' num2str(k) ': Best Cost = ' num2str(Lightbest)]);
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
    bic = conti2bit(nbest);
    cost = Lightbest;
    score = calc_bench(bic);
    
end
% ----- All the subfunctions are listed here ------------
% 初始化萤火虫位置 The initial locations of n fireflies
function [ns,Lightn]=init_ffa(n)
    global Lb Ub dim
    ns=zeros(n,dim);
    for i=1:n
        ns(i,:)=Lb+(Ub-Lb).*rand(1,dim); % 则在取值范围内随机取值
    end

    % 初始化目标函数 initial value before function evaluations
    Lightn=ones(n,1)*10^100;
end
% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,ns,Lightn,nso,Lighto,alpha,betamin,gamma)
% 参数取值范围绝对值 Scaling of the system
    global Lb Ub dim
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
                tmpf=alpha.*(rand(1,dim)-0.5).*scale;
                ns(i,:)=ns(i,:).*(1-beta)+nso(j,:).*beta+tmpf;
            end
        end % end for j
    end % end for i

    % 防止越界 Check if the updated solutions/locations are within limits
    [ns]=findlimits(ns);
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
function s=findlimits(ns)
    global Lb Ub
    Max = max(Ub);
    Min = min(Lb);
    ns(ns>Max) = Max;
    ns(ns<Min) = Min;
    s = ns;
end