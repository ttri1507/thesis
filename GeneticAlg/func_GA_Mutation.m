function Parent_out = func_GA_Mutation(Parent, Pm, bestEle)
Parent_out = Parent;
[x,y,P] = size(Parent);
Prop = rand(1, P);
Par_mut = find(Prop<Pm); % Elements have mutation
% Get rid of bestEle if it is in this list
[logi, loc] = ismember(bestEle, Par_mut);
if(logi)
    Par_mut(loc) = [];
end
numbitchange = ceil(x*y*Pm);
for id = 1:size(Par_mut)
    par = Par_mut(id);
    IDbit_x = randperm(x, numbitchange);
    IDbit_y = randperm(y, numbitchange);
    for i = 1:numbitchange
        if(Parent(IDbit_x(i),IDbit_y(i),par)==0)
            Parent_out(IDbit_x(i),IDbit_y(i),par) = 1;
        else
            Parent_out(IDbit_x(i),IDbit_y(i),par) = 0;
        end
    end
end
end