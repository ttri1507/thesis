function Child = func_GA_Crossover(Parent, P, numUAVeachSAT, numM, numF, Pc, bestEle)
% n = # pairs of chromosomes crossover
% P: population size
% numUAVeachSAT: #UAVs served by each SAT
% numM: The cache capacity (maximum # packages can be prestored at each UAV)
% numF: The number of possible packages
% Parent: 3D population (UAV x numF x elements)
% Pc: the probability of crossover
% bestEle: the element has the best value
Child = zeros(size(Parent));
for par1 = 1:P
    if(par1==bestEle)
        Child(:,:,par1) = Parent(:,:,par1);
    else
        if(Pc>rand) % crossover happens with par with the probability of Pc
            % Choose a different parent to mate
            par2 = par1;
            while(par1==par2)
                par2 = 1 + floor(rand*P);
            end
            % Choose randomly one cut point, the same for all UAV xxxx
            cp = 1 + randi(numF-1);
            % Child 1
            Child(:,1:cp,par1) = Parent(:,1:cp,par1);
            Child(:,cp+1:end,par1) = Parent(:,cp+1:end,par2);
            % Child 2
            Child(:,1:cp,par2) = Parent(:,1:cp,par2);
            Child(:,cp+1:end,par2) = Parent(:,cp+1:end,par1);
        else % crossover doesnt happen
            Child(:,:,par1) = Parent(:,:,par1);
            Child(:,:,par2) = Parent(:,:,par2);
        end
    end
end
end