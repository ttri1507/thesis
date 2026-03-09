function YY = func_GA_InitializePopulation(P, numUAVeachSAT, numM, numF)
% P: population size
% numUAVeachSAT: #UAVs served by each SAT
% numM: The cache capacity (maximum # packages can be prestored at each UAV)
% numF: The number of possible packages
% YY: 3D population (UAV x cacheplacement x elements)

YY = zeros(numUAVeachSAT, numF, P);
for p = 1:P
    Element = zeros(numUAVeachSAT, numF);
    for u = 1:numUAVeachSAT
        numCache = randi(numM);
        IdCache = randperm(numF, numCache); % random ko trung
        Element(u, IdCache) = 1;
    end
    YY(:,:,p) = Element;
end
end
