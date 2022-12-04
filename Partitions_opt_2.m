function [set_part, time_part, time_opt] = Partitions_opt_2(value_opt, init_state, targets, sub_nu)
    tic;
    targets(targets==init_state)=[];
    nu = targets;
    time_opt=[];
    time_part = [];
    sub_n_c_ag = [];
    sub_nu_l = cellfun('length',sub_nu);
    cover_time_opt = -value_opt;
    time = cover_time_opt(init_state,:);
    for i=1:size(sub_nu,1)
        rem_l_j = length(nu) - sub_nu_l(i);
        to_sub_j = find(sub_nu_l==rem_l_j);
        for c_j=1:length(to_sub_j)        
            j = to_sub_j(c_j);
            if(length(unique([sub_nu{i,1}(1:sub_nu_l(i)),sub_nu{j,1}(1:sub_nu_l(j))]))==length(nu))
                if(isempty(time_opt))
                    time_opt = max([time(i),time(j)]);
                    time_part = [time(i),time(j)];
                    set_part_1 = sub_nu{i,1}(1:sub_nu_l(i));
                    set_part_2 = sub_nu{j,1}(1:sub_nu_l(j));
                elseif(max([time(i),time(j)]) < time_opt)
                    time_opt = max([time(i),time(j)]);
                    time_part = [time(i),time(j)];
                    set_part_1 = sub_nu{i,1}(1:sub_nu_l(i));
                    set_part_2 = sub_nu{j,1}(1:sub_nu_l(j));
                end
            end
        end
    end
    set_part = {set_part_1; set_part_2};
    toc;
end