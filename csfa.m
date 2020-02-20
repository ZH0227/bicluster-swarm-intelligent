function [bic,cost,score,history]=csfa(nPop,lamda,miu,omega)
    global Lb Ub dim
    
    pa=0.25;
    alpha=0.25;
    ShowIterInfo =1;
    %% Change this if you want to get better results
    N_IterTotal=300;
    early_stopping_cnt = 0;
    early_stopping_maxcnt = 40;
    
    % Random initial solutions
    nest = zeros(nPop,dim);
    for i=1:nPop
        nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
    end
    
    % Get the current best
    fitness=10^10*ones(nPop,1);
    [fmin,best,fa,fitness]=get_best_nest(nest,nest,fitness,lamda,miu,omega);
    history = zeros(N_IterTotal,1);

    fmin_old = min(fitness);
    for i =1:N_IterTotal
        % [fmin_old,~,~,~]=get_best_nest(nest,nest,fitness,data,lamda,miu,omega,c2bT,costFun);
        [nest,fitnessC,bestC,minC]=cs_iter(fa,fitness,pa,lamda,miu,omega);
        [fa,fitness,bestF,minF]=fa_iter(nest,fitnessC,lamda,miu,omega,N_IterTotal,alpha);

        if minC < minF
            fmin = minC;
            best = bestC;
            %repalce
            [~,index] = min(fitness);
            fa(index,:) = bestC;
        else
            fmin = minF;
            best = bestF;
        end

        % early stopping
        change = fmin_old - fmin;
        if change < fmin_old/10000
            early_stopping_cnt = early_stopping_cnt + 1;
        else
            early_stopping_cnt = 0;
        end
        fmin_old = fmin;

        history(i) = fmin;
        if ShowIterInfo
            disp(['csfab=> ', 'Iteration ' num2str(i) ': Best Cost = ' num2str(fmin)]);
        end
        if early_stopping_cnt > early_stopping_maxcnt
            disp('csfab early stoping');
            break
        end
    end % end of iteration
    if i < N_IterTotal
        for k =i:N_IterTotal
            history(k) = fmin;
        end
    end
    bic = conti2bit(best);
    cost = fmin;
    score = calc_bench(bic);
end

%% --------------- All subfunctions are list below ------------------
function [nest,fitness,best,fmin]=cs_iter(nest,fitness,pa,lamda,miu,omega)
%
%
    [~,index]=min(fitness);
    bestnest = nest(index,:);
    new_nest=get_cuckoos(nest,bestnest);   
    [fmin,~,nest,fitness]=get_best_nest(nest,new_nest,fitness,lamda,miu,omega);

    % Discovery and randomization
    new_nest=empty_nests(nest,pa) ;

    % Evaluate this set of solutions
    [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,lamda,miu,omega);

    % Find the best objective so far  
    if fnew<fmin
    fmin=fnew;
    end
end

function [fa,fitness,best,fmin]=fa_iter(ns,zn,lamda,miu,omega,N_IterTotal,alpha)
    global costFun
    % 更新alpha（可选）This line of reducing alpha is optional
    alpha=alpha_new(alpha,N_IterTotal);
    betamin=0.20;
    gamma=1;
    n=size(ns,1);
    d=size(ns,2);
    % 根据亮度排序 Ranking fireflies by their light intensity/objectives
    [Lightn,Index]=sort(zn);
    ns_tmp=ns;
    for i=1:n
        ns(i,:)=ns_tmp(Index(i),:);
    end

    
    nso=ns;         %Old
    Lighto=Lightn;  %Old
    % nbest=ns(1,:);
    % Lightbest=Lightn(1);

    % 向较优方向移动 Move all fireflies to the better locations
    [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,alpha,betamin,gamma);

    % 对每个萤火虫计算目标函数值 Evaluate new solutions (for all n fireflies)
    for i=1:n
        bic = conti2bit(ns(i,:));
        zn(i)=costFun(bic,lamda,miu,omega);
    end
    %% 找出当前最优 Find the current best
    [fmin,index]=min(zn);
    fa=ns;
    fitness=zn;
    best=ns(index,:);
end


%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best)
    % Levy flights
    n=size(nest,1);
    % Levy exponent and coefficient
    % For details, see equation (2.21), Page 16 (chapter 2) of the book
    % X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    
    for j=1:n
        s=nest(j,:);
        % This is a simple way of implementing Levy flights
        % For standard random walks, use step=1;
        %% Levy flights by Mantegna's algorithm
        u=randn(size(s))*sigma;
        v=randn(size(s));
        step=u./abs(v).^(1/beta);
      
        % In the next equation, the difference factor (s-best) means that 
        % when the solution is the best solution, it remains unchanged.     
        stepsize=0.01*step.*(s-best);
        % Here the factor 0.01 comes from the fact that L/100 should the typical
        % step size of walks/flights where L is the typical lenghtscale; 
        % otherwise, Levy flights may become too aggresive/efficient, 
        % which makes new solutions (even) jump out side of the design domain 
        % (and thus wasting evaluations).
        % Now the actual random walks or flights
        s=s+stepsize.*randn(size(s));
       % Apply simple bounds/limits
       nest(j,:)=simplebounds(s);
    end
end    
    %% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness,lamda,miu,omega)
    % Evaluating all new solutions
    global costFun
    for j=1:size(nest,1)
        bic = conti2bit(newnest(j,:));
        fnew=costFun(bic,lamda,miu,omega);
        if fnew<=fitness(j)
           fitness(j)=fnew;
           nest(j,:)=newnest(j,:);
        end
    end
    % Find the current best
    [fmin,K]=min(fitness) ;
    best=nest(K,:);
end 
    %% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,pa)
    % A fraction of worse nests are discovered with a probability pa
    n=size(nest,1);
    % Discovered or not -- a status vector
    K=rand(size(nest))>pa;
    
    % In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
    % this cuckoo's egg is less likely to be discovered, thus the fitness should 
    % be related to the difference in solutions.  Therefore, it is a good idea 
    % to do a random walk in a biased way with some random step sizes.  
    %% New solution by biased/selective random walks
    stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
    new_nest=nest+stepsize.*K;
    for j=1:size(new_nest,1)
        s=new_nest(j,:);
      new_nest(j,:)=simplebounds(s);  
    end
end 
% Application of simple constraints
function s=simplebounds(s)
    global Lb Ub
    % Apply the lower bound
    ns_tmp=s;
    I=ns_tmp<Lb;
    ns_tmp(I)=Lb(I);
    
    % Apply the upper bounds 
    J=ns_tmp>Ub;
    ns_tmp(J)=Ub(J);
    % Update this new move 
    s=ns_tmp;
end
% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,alpha,betamin,gamma)
% 参数取值范围绝对值 Scaling of the system
    global Lb Ub
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
    [ns]=findlimits(n,ns);
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
function [ns]=findlimits(n,ns)
    global Lb Ub
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