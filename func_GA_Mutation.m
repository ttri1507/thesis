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
if(~isempty(Par_mut))
for id = 1:size(Par_mut)
    par = Par_mut(id);
    IDbit = randperm(x*y, numbitchange);
    for i = 1:numbitchange
        xi = mod(IDbit(i),x); 
        if(xi==0) 
            xi = x; 
        end
        yi = ceil(IDbit(i)/x);
        if(Parent(xi,yi,par)==0)
            Parent_out(xi,yi,par) = 1;
        else
            Parent_out(xi,yi,par) = 0;
        end
    end
end
end
end