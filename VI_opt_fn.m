function [value_opt, a_opt, sub_nu, sub_n] = VI_opt_fn(P, init_state, targets)
    n_S = size(P,1); % Number of states
    n_A = size(P,2); % Number of actions
    nu = double(targets);
    nu(nu==init_state)=[];
    n_nu = length(nu);
    sub_nu ={}; % All Subsets of nu
    sub_n =[];
    for i=1:n_nu
        sub_t = nchoosek(nu,i);
        sub_nu = [sub_nu; [mat2cell(sub_t,ones(1,size(sub_t,1)))]];
        sub_n = [sub_n; size(sub_t,1)];
    end
    n_sub_nu = size(sub_nu,1);
    value_opt = repmat(0, n_S, n_sub_nu);
    a_opt = zeros(size(value_opt,1),size(value_opt,2));
    delta = 1000;
    gamma=1;
    sf=0;
    iter=0;
    tic;
    while(sf==0)
        iter=iter+1;
        value_1 = value_opt;
        a_opt_1 = a_opt;
        for curr_s = 1:size(value_opt,1)
            for curr_sub_i = 1:size(value_opt,2)
                curr_subset = sub_nu{curr_sub_i,1};
                curr_subset(curr_subset==curr_s)=[];
                if(~isempty(curr_subset))
                    P_curr_s = reshape(P(curr_s,:,:),n_A,n_S);
                    to_actions = find(sum(P_curr_s,2)>0);
                    V_t_max = -inf;
                    for c_a=1:length(to_actions)
                        a = to_actions(c_a);
                        to_states = find(P(curr_s,a,:)>0);
                        add = 0;
                        for c_to_s = 1:length(to_states)
                            to_s = to_states(c_to_s);
                            to_subset = curr_subset;
                            to_subset(to_subset==to_s)=[];
                            l_to_subset = length(to_subset);
                            if(~isempty(to_subset))
                                if(l_to_subset==1)
                                    start_search=1;
                                else
                                    start_search = 1+sum(sub_n([1:l_to_subset-1]));
                                end
                                search = start_search;
                                end_search = sum(sub_n([1:l_to_subset]));
                                while(~all(sub_nu{search}==to_subset))
                                    search = search+1;
                                end
                                add = add + P(curr_s, a, to_s)*value_opt(to_s, search);
                            end
                        end
                        if(-1+gamma*add > V_t_max)
                            a_max = a;               
                            V_t_max = -1 + gamma*add;
                        end
                    end
                    value_opt(curr_s, curr_sub_i) = V_t_max;
                    a_opt(curr_s, curr_sub_i) = a_max;
                else
                    value_opt(curr_s, curr_sub_i) = 0;
                end
            end
        end
        if(all(a_opt==a_opt_1))
            sf=1;
        end
    end
    toc;
end