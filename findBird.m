function ind_row = findBird(bird_name)
    if strcmpi(bird_name,'pu1')
        ind_row = 3;
    elseif strcmpi(bird_name,'ye3')
        ind_row = 4;
    elseif strcmpi(bird_name,'or3')
        ind_row = 5;
    elseif strcmpi(bird_name,'bl3')
        ind_row = 6;
    elseif strcmpi(bird_name,'bl4')
        ind_row = 7;
    elseif strcmpi(bird_name,'or4')
        ind_row = 8;
    else
        error('bird name not known')
    end
end