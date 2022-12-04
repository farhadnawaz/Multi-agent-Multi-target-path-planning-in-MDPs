function [time_v, time_all] = VI_heur_multi(P,init_state,targets,set_part_heur,noA) % for VI_heur implementation
    thr = 10^(-10); 
    gamma = 0.01;
    time = []; % time taken by heuristic (voronoi) partitions (mean max)
    tic;
    nu = double(targets); % Targets of MDP c_g
    n_nu = length(nu);
    time_all=[];
    for c=1:100 % to compute averge cover time of VI_heur
        time_a=[];
        for ca=1:noA
            nu_a = set_part_heur{ca};
            n_nu_a = length(nu_a);
            V_t = repmat(0,size(P,1),1);
            states_v = init_state;
            curr_state = init_state;
            curr_subset = setdiff(nu_a,curr_state);
            % Until agent covers the targets
            while(~(length(intersect(states_v,nu_a))==n_nu_a) || ~all(intersect(states_v,nu_a)==unique(nu_a)))
                [V_t] = VI_heur_fn(P,curr_subset,curr_state,gamma,thr,V_t); % compute the heuristic value function
                % Implement the policy from heuristic value function until the agent visits a target
                [curr_subset,states_v,curr_state] = heur_impl_fn(P,V_t,curr_subset,curr_state,gamma,states_v);
            end 
            time_a = [time_a length(states_v)-1];
        end
        time_all = [time_all; time_a];
        time = [time max(time_a)];
    end
    time_v = mean(time); 
    % load("D:\UIUC\RA\Cover_Time\Matlab\Data_files\Random_MDPs_paper_ACC\Partitions_opt\MDP_10_all_part_3.mat"); % Optimal partition
    % time_opt = cell2mat(part(:,4));
    % sub_opt = ((time_v - time_opt)./time_opt)*100;
    toc;
end