function B = FUNC_RandomCache(numUAV, numM, numF)

B = zeros(numUAV, numF);
for u = 1:numUAV
    numCache = randi(numM);
    IdCache = randperm(numF, numCache); % random ko trung
    B(u, IdCache) = 1;
end

end

