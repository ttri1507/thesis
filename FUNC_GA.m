function [B_full, BestFit] = FUNC_GA(P, Pc, Pm,  s, MODEL, A, B, Area, Q, Noise_var,...
    BW, numUAVeachSAT, numUEeachUAV, numM, numF, numSAT, numUAV, numUE, PmaxSAT, PmaxUAV)
% Stop condition
max_generation = 20;
tolerance = 5e-3; % second
max_tole_generation = 5; % stop after the distance between two gene..s < tole

% Initialize
% P = 40; % #elements in the population
% Pc = 0.7; % The probability of crossover
% Pm = 0.3; % The probability of mutation

% Initialize the population
Parent = func_GA_InitializePopulation(P, numUAVeachSAT, numM, numF);

% Evaluate
% s = 1; % considering satellite
B_full = B;
[Fitness, Fit_min, bestEle] = func_GA_Evaluate(Parent, P, MODEL, A, B_full,...
              Area, Q, Noise_var, BW, s, numUAVeachSAT, numUEeachUAV, numM, numSAT, numUAV, numUE, PmaxSAT, PmaxUAV);
BestFit = Fit_min;

terminal = 0; generation = 0; tole_generation = 0;
while ~terminal
    generation = generation + 1;
    % Selection
    Parent = func_GA_SelectParent(Parent, P, Fitness, bestEle);
    % Crossover to generate new population
    Parent = func_GA_Crossover(Parent, P, numUAVeachSAT, numM, numF, Pc, bestEle);
    % Mutation
    Parent = func_GA_Mutation(Parent, Pm, bestEle);
    
    % Evaluate again
    [Fitness, Fit_min, bestEle] = func_GA_Evaluate(Parent, P, MODEL, A, B_full,...
              Area, Q, Noise_var, BW, s, numUAVeachSAT, numUEeachUAV, numM, numSAT, numUAV, numUE, PmaxSAT, PmaxUAV);
    BestFit = [BestFit, Fit_min];
    % Check stop condition
    if(generation == max_generation)
        terminal = 1;
    elseif(generation>1)
        if(abs(Fit_min - BestFit(generation-1))<tolerance)
            tole_generation = tole_generation + 1;
            if(tole_generation==max_tole_generation)
                terminal = 1;
            end
        else
            tole_generation = 0;
        end
    end
end
B_opt = Parent(:,:,bestEle);
Us = numUAVeachSAT*(s-1) + 1:numUAVeachSAT*s;
B_full(Us,:) = B_opt;
end