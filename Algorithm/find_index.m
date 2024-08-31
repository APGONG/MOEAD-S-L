function index = find_index(pop, N)
    minv = inf;
    for i = 1 : N
        if pop(i).objs < minv
            minv = pop(i).objs;
            index = i;
        end
    end
end
