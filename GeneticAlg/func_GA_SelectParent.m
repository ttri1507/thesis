function Parent_out = func_GA_SelectParent(Parent, P, Fitness, bestEle)
% Kill parent does not meet the constraint tk > 1e6
Id_die = find(Fitness > 1e6);
Id_live = 1:P; Id_live(Id_die) = [];
Fitness(Id_die) = [];
% Compute the probability and CDF matrix, Note: minimizing problem
Fitness_new = max(Fitness) - Fitness;
Prob = Fitness_new/sum(Fitness_new);
Cdf = zeros(1,P+1);
for par = 1:P %xxxx maybe can use the linear in ref
    Cdf(par+1) = Cdf(par) + Prob(par);
end
% Selection acording to CDF
Parent_out(:,:,1) = Parent(:,:,bestEle); % Always keep the best element
for par = 2:P
    Index = find(Cdf<rand); id_chosen = Id_live(Index(end));
    Parent_out(:,:,par) = Parent(:,:,id_chosen);
end
end