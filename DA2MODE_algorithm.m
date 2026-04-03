% DA2MODE_algorithm.m
% Official MATLAB implementation of the DA2MODE algorithm
%
% Paper:
% Reda, M., Onsy, A., Haikal, A.Y., and Ghanbari, A.
% DA2MODE: Dynamic Archive with Adaptive Multi-Operator Differential
% Evolution for numerical optimization.
% Swarm and Evolutionary Computation, 99 (2025), 102130.
% DOI: 10.1016/j.swevo.2025.102130
%
% The DA2MODE algorithm was validated on:
% - CEC2017/2018 benchmark functions
% - CEC2020/2021 benchmark functions
% - CEC2022 benchmark functions
%
% Usage:
%   [goalReached, GlobalBest, countFE] = DA2MODE_algorithm();
%   [goalReached, GlobalBest, countFE] = DA2MODE_algorithm(CECyear, fNo, nd, lb, ub);
%
% Default values:
%   CECyear = 2020
%   fNo     = 3
%   nd      = 20
%   lb      = -100
%   ub      = 100
%
% Notes:
% - Requires the corresponding CEC benchmark wrapper files and compiled MEX files
% - Requires the correct benchmark input-data folders in the repository
%
% Author:
% Dr. Mohamed Reda
% University of Central Lancashire, UK
% Mansoura University, Egypt
%
% Contact:
% mohamed.reda.mu@gmail.com
% mohamed.reda@mans.edu.eg

function [goalReached, GlobalBest, countFE] = DA2MODE_algorithm(CECyear, fNo, nd, lb, ub)

   %% Begin intialiazation
    %% CEC function paramters

     % Default values if no inputs are passed
    if nargin < 1 || isempty(CECyear)
        CECyear = 2020;
    end

    if nargin < 2 || isempty(fNo)
        fNo = 3;
    end

    if nargin < 3 || isempty(nd)
        nd = 20;
    end

    if nargin < 4 || isempty(lb)
        lb = -100;
    end

    if nargin < 5 || isempty(ub)
        ub = 100;
    end

    %% Cost Function
    if CECyear == 2020
        costFunction = @(sol) cost_cec2020(sol, fNo);
    elseif CECyear == 2022
        costFunction = @(sol) cost_cec2022(sol, fNo);
    elseif CECyear == 2017
        costFunction = @(sol) cost_cec2017(sol, fNo);
    else
        error('Unsupported CEC year. Use 2017, 2020, or 2022.');
    end

    %% Algorithm paramters
    AlgName = 'DA2MODE Algorithm' ;
    n_opr = 3 ; % number of operators 

    arch_rate = 1.3; % 1.3 (A) rate of the archive control the size of it
    mem_size_scale = 6 ; % memory size (H) = 20 * nd, where 20 is the memeory scale
    
    init_CR = 0.8; % Initial CR
    init_F = 0.3; % initial F

    MinPopSize = 4 ; % min pop size (linear reduction for population)
   
    pArchMin = 0.5; % min archive rate (pArch is the pArch% worst of the popuation that will be added to the archive)
    init_pArch = 1 ; % intial archive rate (intial pArch)
    alphaEMA = 0.6 ; % EMA smoothing factor: is the ratio of the iteration success rate (for operators) to the accumarated success rate over all iteration EMA

    % intial pop size 
    nPop =    ceil(60 * nd ^(2/3)); % ceil(75*dim^(2/3));
    InitPop = nPop ;

    % prob. of each DE operator
    probDE1= 1./n_opr .* ones(1,n_opr);

    %% Initialize Archive Data 
    archive.NP = arch_rate * nPop; % the maximum size of the archive
    cand = create_empty_individual(); 
    cand.Position = zeros(1, nd); % empty solution of nd dim
    cand.Cost = inf(1, 1); % empty cost value
    archPop = repmat(cand, 1, 1);  % popualtion of one individual
    archive.pop = archPop; 

    %% Initialize Adaptive Archive for CR and F
    hist_pos=1;
    memory_size=mem_size_scale*nd;
    archive_f = ones(1,memory_size).* init_F;
    archive_Cr = ones(1,memory_size).* init_CR;

    %% initialize operators scores (weights)
    % EMA (Exponential Moving Average) parameters
    EMA_S = zeros(1, n_opr); % Initial EMA scores for each operator
    %EMA_S_exp = zeros(1, n_opr_explor); % Initial EMA scores for each exploration operator

    %% Stopping criteria
    tol = 10^-8;
    if CECyear == 2020 || CECyear == 2022
        if nd == 10
            maxfe = 200000 ; 
        elseif nd == 20 
            maxfe = 500000;
        else
            fprintf('Dimensions must be 10 or 20 \n');
            return
        end
    elseif CECyear == 2017
        maxfe = nd * 10000;
    end

    %% Display iteration prompt
    print_flag = true;

    %%  Global variable to count number of function evaluations
    global countFE;
    countFE = 0 ;

    %% Initialize iteration counter 
    N_iter = 0;

    %% Goal reached flag
    goalReached = false; 

    %% Initialize Global best
    GlobalBest.Position = [];
    GlobalBest.Cost = Inf ; 

    %% Set the seed for random number generator
    rng('default');  % Resets to the default settings
    rng('shuffle'); % set it to shuffle
    
    %% Initialize population and update the global best
    population = repmat(create_empty_individual(), nPop, 1);

    % Generate initial positions using Latin Hypercube Sampling
    LB = repmat(lb, 1, nd);
    UB = repmat(ub, 1, nd);
    LHS_samples = lhsdesign(nPop, nd); % Generates an nPop x CEC.dim matrix of samples
    % Scale samples to the problem's bounds
    for d = 1:nd
        LHS_samples(:, d) = LB(d) + (UB(d) - LB(d)) .* LHS_samples(:, d);
    end
    
    for i = 1:nPop
        % Initialize Position with scaled LHS samples
        if rand <=1
            population(i).Position = LHS_samples(i, :);
        else
            population(i).Position = lb + (ub - lb) .* rand(1, nd);
        end
        
        % Evaluation of the cost
        population(i).Cost = costFunction(population(i).Position);

        % Update Global Best
        if population(i).Cost < GlobalBest.Cost
            GlobalBest = population(i);
        end
    end

    % create a random permuatation of the popualtion
    pop_old = population(randperm(nPop),:);

   %% begin algorithm loop 
    while (GlobalBest.Cost > tol)  && (countFE <= maxfe)
        %% update the generation
        N_iter=N_iter+1; 

        %% Update popuation size (Linear)
        % EIR = 1 - ((GlobalBest.Cost - tol) / (InitError - tol));
        % FERate = (countFE / maxfe) ; 
        % TR = 0.4 * EIR + 0.6 * FERate; 
        % newPopSize= round(MinPopSize + ((InitPop - MinPopSize) * (1 - ( TR).^1)));

        newPopSize = round((((MinPopSize - InitPop) / maxfe) * countFE) + InitPop);

        %% Update the popuation according to the new popuation size
        nPop = numel(population); % current popsize 
        if nPop > newPopSize
            % Calculate the number of individuals to remove
            reduction_ind_num = nPop - newPopSize;
            if nPop - reduction_ind_num < MinPopSize
                reduction_ind_num = nPop - MinPopSize;
            end
        
            % Remove the worst individuals
            for r = 1 : reduction_ind_num
                % Sort population based on Cost
                [~, sortedIdx] = sort([population.Cost], 'descend');
                % Remove the worst individual
                population(sortedIdx(1)) = []; % it removed from the original popualtion
            end

            % update the current popsize
            nPop = numel(population);
        
            %% Update archive size based on the new population size
            archive.NP = round(arch_rate * nPop);

            % If archive size exceeds its limit, randomly remove some individuals
            current_archive_NP = numel(archive.pop);
            if current_archive_NP > archive.NP
                rndpos = randperm(current_archive_NP);
                rndpos = rndpos(1 : archive.NP);
                archive.pop = archive.pop(rndpos);
            end
        end

        %% Initialize the archive of the CR and F
        mem_rand_index = ceil(memory_size * rand(nPop, 1));
        mu_sf = archive_f(mem_rand_index);
        mu_cr = archive_Cr(mem_rand_index);
        
        %%  generate CR   
        cr = normrnd(mu_cr, 0.1);
        term_pos = find(mu_cr == -1);
        cr(term_pos) = 0;
        cr = min(cr, 1);
        cr = max(cr, 0);
        % sort the cr
        [cr,~]=sort(cr);

        %% for generating scaling factor
        F = mu_sf + 0.1 * tan(pi * (rand(1,nPop) - 0.5));
        pos = find(F <= 0);
        
        while ~ isempty(pos)
            F(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(1,length(pos)) - 0.5));
            pos = find(F <= 0);
        end
        
        F = min(F, 1);
        F=F';

         %% Sort the popuation     
        Costs = [population.Cost]; % original costs of the original popuation
        [Costs, SortOrder] = sort(Costs);
        population = population(SortOrder);  

        %% **** Mutation Phase ****
        % combine the popuation with the archive population  
        popAll = [population; archive.pop];  
        
        %% generate mutation operator probablities for each individual in the population
        % Randomly decide the mutation strategy for each individual
        bb = rand(nPop, 1);
        
        % Retrieve probabilities for each strategy
        probiter = probDE1(1, :);
        l2 = sum(probDE1(1:2));
        
        % Determine which strategy to apply for each individual
        op_1 = bb <= probiter(1) * ones(nPop, 1);
        op_2 = (bb > probiter(1)) & (bb <= l2);
        op_3 = (bb > l2) & (bb <= 1);

        %% generate random integer numbers
        r0 = 1 : nPop;
        [r1, r2,r3] = gnR1R2(nPop, size(popAll, 1), r0);

        %% Choose top individuals (at least one) for DE operator 1 and 2
        pNP12 = max(round(0.25 * nPop), 1); % At least one or 25% of the population size
        randindex = ceil(rand(1, nPop) .* pNP12); % Select indices from the best subset
        randindex = max(1, randindex); % Ensuring indices are valid (not less than 1)
        phix12 = population(randindex, :);
        
        %%  Choose top individuals (at least two) for DE operator 3
        pNP3 = max(round(0.5 * nPop), 2); %% choose at least two best solutions
        randindex = ceil(rand(1, nPop) .* pNP3); %% select from [1, 2, 3, ..., pNP]
        randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
        phix3 = population(randindex, :);

        %% Initialize mutation vector
        cand = create_empty_individual();
        cand.Position = zeros(1, nd);
        newPop = repmat(cand, nPop, 1);

        for i = 1:nPop
            % apply mutation rule
            x_curr = population(i).Position;
            x_r1 = population(r1(i)).Position;
            x_r3 = population(r3(i)).Position;
            xx_r2 = popAll(r2(i)).Position;
            x_phi12 = phix12(i).Position;
            x_phi3 = phix3(i).Position;

             if op_1(i)
                % Strategy 1: Mutation based on the difference between two individuals
                newPop(i).Position = x_curr + F(i) * (x_phi12 - x_curr + x_r1 - xx_r2);
            elseif op_2(i)
                % Strategy 2 : Mutation using three random individuals
                newPop(i).Position = x_curr + F(i) * (x_phi12 - x_curr + x_r1 - x_r3);
            elseif op_3(i)
                % Strategy 3 (DE3)
                newPop(i).Position = F(i) * x_r1 + F(i) * (x_phi3 - x_r3);
             end

             % Applying boundary check for each individual
             newPop(i).Position = han_boun_individual(newPop(i).Position, ub, lb, x_curr);
        end

        %% *** Crossover ***
        % Initialize ui as an array of individuals
        cand = create_empty_individual();
        cand.Position = zeros(1, nd);
        newPop2 = repmat(cand, nPop, 1);

        for i = 1:nPop
            % Generate a random number to decide the crossover method
            if rand < 0.4
                % Binomial Crossover
                mask = rand(1, nd) > cr(i);
                jrand = floor(rand * nd) + 1; % Ensure at least one dimension is inherited from vi
                mask(jrand) = false;
                newPop2(i).Position = population(i).Position; % Start with parent position from origial population
                newPop2(i).Position(~mask) = newPop(i).Position(~mask); % Inherit from newPop where mask is false
            else
                % Exponential Crossover
                startLoc = randi(nd);
                L = 0;
                while (rand < cr(i) && L < nd)
                    L = L + 1;
                end
                idx = mod(startLoc-1:startLoc+L-2, nd) + 1; % Ensuring wrapping around dimensions
                newPop2(i).Position = population(i).Position; % Start with parent position
                newPop2(i).Position(idx) = newPop(i).Position(idx); % Inherit from vi for selected indices
            end

            %% evaluate the new cost
            newPop2(i).Cost = costFunction(newPop2(i).Position);   
        end

        %% *** Update the archives ****
        %% get the I label for the improved individuals
        newCosts = [newPop2.Cost] ;
        I = (newCosts < Costs ); % Logical index of improved solutions

        %% update the archive with the old bad solutions in the population
        pArch = round((((pArchMin - init_pArch) / maxfe) * countFE) + init_pArch);
        archive = updateArchive_basic(archive, population(I == 1), pArch) ; % pass the worst

        %% update probDE1 (operators probabilities). of each DE
        diff2 = max(0, (Costs - newCosts))./abs(Costs + eps); % Improvement metric, adding eps for stability

        % Calculate performance scores for this iteration
        count_S = zeros(1, n_opr);
        count_S(1)=max(0,mean(diff2(op_1==1)));
        count_S(2)=max(0,mean(diff2(op_2==1)));
        count_S(3)=max(0,mean(diff2(op_3==1)));

        % Check if there is any significant improvement across all operators
        if all(count_S <= eps) % If no significant improvement, reset to equal probabilities
            probDE1 = ones(1, n_opr) / n_opr;
        else
            % % Update EMA commulative scores over time
            for i = 1:n_opr
                EMA_S(i) = alphaEMA * count_S(i) + (1 - alphaEMA) * EMA_S(i);
            end
            % Scale EMA_S scores to magnify differences
            scalingFactor = 1000; % Example scaling factor, adjust based on your specific problem
            scaled_EMA_S = EMA_S * scalingFactor;
    
            % Apply softmax to scaled EMA_S scores
            probDE1 = exp(scaled_EMA_S) / sum(exp(scaled_EMA_S));
    
            % Optionally enforce minimum and maximum probabilities (e.g., between 0.1 and 0.9)
            % This step is typically not necessary right after softmax, but included for completeness
            probDE1 = max(0.1, min(0.9, probDE1));

        end
        %% calc. imprv. for Cr and F archives
        goodCR = cr(I == 1);
        goodF = F(I == 1);
        diff = abs(Costs - newCosts);
        if size(goodF,1)==1
            goodF=goodF';
        end
        if size(goodCR,1)==1
            goodCR=goodCR';
        end
        num_success_params = numel(goodCR);
        if num_success_params > 0
            weightsDE = diff(I == 1)./ sum(diff(I == 1));
            %% for updating the memory of scaling factor
            archive_f(hist_pos) = (weightsDE * (goodF .^ 2))./ (weightsDE * goodF);
            
            %% for updating the memory of crossover rate
            if max(goodCR) == 0 || archive_Cr(hist_pos)  == -1
                archive_Cr(hist_pos)  = -1;
            else
                archive_Cr(hist_pos) = (weightsDE * (goodCR .^ 2)) / (weightsDE * goodCR);
            end
            
            hist_pos= hist_pos+1;
            if hist_pos > memory_size;  hist_pos = 1; end
        else
            archive_Cr(hist_pos)=0.5;
            archive_f(hist_pos)=0.5;
        end

        %% update population with the good solution and update the oldPop with the bad old solutions
        pop_old(I == 1) = population(I == 1); %save the bad individual in popualtion in the old popuation
        population(I == 1) = newPop2(I == 1); % relace the bad individuals in popuation with the better individuals in newPop2

        %% sort the population and old population
        [~, sortedIndices] = sort([population.Cost]); % the new updated merged costs
        population = population(sortedIndices);
        pop_old = pop_old(sortedIndices);

        %% update global best
        localBest = population(1); % the sorted population, the top solution is the best cost solution
        if localBest.Cost < GlobalBest.Cost  
            GlobalBest = localBest;
        end

        %% check if maxfes is exceeded 
        if countFE > maxfe 
            break;
        end

        %% print the iteration number
        if print_flag            
            fprintf('%s | CEC%d_F%d_D%d | Iteration %d |FEs %d | Error %d\n', AlgName,  CECyear, fNo , nd, N_iter, countFE,GlobalBest.Cost);
        end

         %check tolerance /error
        if (GlobalBest.Cost <= tol)
            GlobalBest.Cost  = 0 ;
            disp('tol reached');
            goalReached = true ; 
            break ;  % not needed, becuase it will exit in the next while loop check
        end

    end
end

function individual = create_empty_individual()
    individual.Position = [];
    individual.Cost = Inf;
end

function archive = updateArchive_basic(archive, newPop, p)
    % Update the archive with input solutions
    % Input:
    %   archive - The existing archive with fields 'NP' and 'pop'
    %   newPop - Array of new individuals to add to the archive
    
    if archive.NP == 0, return; end

    nPop = numel(newPop); 
    n_selected_Arch = round(p*nPop);

    [~, sortedIdx] = sort([newPop.Cost],'descend');
    % Remove the worst individual
    worstP = newPop(sortedIdx(1:n_selected_Arch)) ;% add the p% worst of the popuation

    % Combine existing archive population with new population
    combinedPop = [archive.pop; worstP];

    % Randomly remove solutions if necessary to maintain archive size
    nA = numel(combinedPop); 
    Asize = archive.NP;
    if nA > Asize
        combinedPop = combinedPop(floor(nA-Asize+1):nA);     
    end
    % return the uprated archive
    archive.pop = combinedPop;
end


function [r1, r2,r3] = gnR1R2(NP1, NP2, r0)
  
    NP0 = length(r0);
    r1 = floor(rand(1, NP0) * NP1) + 1;
    
    for i = 1 : 99999999
        pos = (r1 == r0);
        if sum(pos) == 0
            break;
        else % regenerate r1 if it is equal to r0
            r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r1 in 1000 iterations');
        end
    end
    
    r2 = floor(rand(1, NP0) * NP2) + 1;
    %for i = 1 : inf
    for i = 1 : 99999999
        pos = ((r2 == r1) | (r2 == r0));
        if sum(pos)==0
            break;
        else % regenerate r2 if it is equal to r0 or r1
            r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end
    
    r3= floor(rand(1, NP0) * NP1) + 1;
    %for i = 1 : inf
    for i = 1 : 99999999
        pos = ((r3 == r0) | (r3 == r1) | (r3==r2));
        if sum(pos)==0
            break;
        else % regenerate r2 if it is equal to r0 or r1
             r3(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end
end

function x = han_boun_individual(x, ub, lb, x2)
    if isscalar(ub)
        ub = repmat(ub, 1, numel(x2));
    end
    if isscalar(lb)
        lb = repmat(lb, 1, numel(x2));
    end

    x_L = lb;
    pos = x < x_L;
    x(pos) = (x2(pos) + x_L(pos)) / 2;

    x_U = ub;
    pos = x > x_U;
    x(pos) = (x2(pos) + x_U(pos)) / 2;
end


