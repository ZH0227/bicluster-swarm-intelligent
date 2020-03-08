function [bic,cost,score,history]=csb(nPop,lamda,miu,omega)
    global Lb Ub dim
    
    pa=0.25;
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
    [fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness,lamda,miu,omega);
    history = zeros(N_IterTotal,1);

    %% Starting iterations
    for iter=1:N_IterTotal
        % Generate new solutions (but keep the current best)
         new_nest=get_cuckoos(nest,bestnest);   
         [~,~,nest,fitness]=get_best_nest(nest,new_nest,fitness,lamda,miu,omega);

         % Discovery and randomization
          new_nest=empty_nests(nest,pa) ;
        
        % Evaluate this set of solutions
          [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,lamda,miu,omega);

        % early stopping
        change = fmin - fnew;
        if change < fmin/10000
            early_stopping_cnt = early_stopping_cnt + 1;
        else
            early_stopping_cnt = 0;
        end
        % Find the best objective so far  
        if fnew<fmin
            fmin=fnew;
            bestnest=best;
        end
        history(iter) = fmin;
        if ShowIterInfo
            disp(['csb=> ','Iteration ' num2str(iter) ': Best Cost = ' num2str(fmin)]);
        end

        if early_stopping_cnt > early_stopping_maxcnt
            disp('csb early stoping');
            break
        end

    end %% End of iterations
    if iter < N_IterTotal
        for i =iter:N_IterTotal
            history(i) = fmin;
        end
    end
    cost = fmin;
    bic = conti2bit(bestnest);
    score = calc_bench(bic);

end
%% --------------- All subfunctions are list below ------------------
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
    bics = conti2bit(newnest);
    fnews=costFun(bics,lamda,miu,omega);
    for j=1:size(nest,1)
        fnew = fnews(j);
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

