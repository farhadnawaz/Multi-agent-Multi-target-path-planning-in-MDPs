function [V_t] = VI_heur_fn(P,curr_subset,curr_state,gamma,thr,V_t)
    n_S = size(P,1); % Number of states
    n=sqrt(n_S);
    n_A = size(P,2); % Number of actions
    time=[];
    pi_one = [];
    l_curr_subset = length(curr_subset);   
    V_1 = V_t;
    %% Decision making
    delta = 1000;
    while(delta>thr)
        for i=1:n_S
            V_t_max = -100000;
            to_actions = find(any(P(i,1:n_A,:),3)>0);
            for c_a=1:length(to_actions)
                k = to_actions(c_a);
                add = 0;
                to_states = find(P(i,k,:)>0);
                for c_s=1:length(to_states)
                    j = to_states(c_s);
                    if(any(j==curr_subset))
                        reward_t = -(l_curr_subset - 1);
                    else
                        reward_t = -l_curr_subset;
                    end
                    add = add + P(i,k,j)*(reward_t + gamma*V_t(j));
                end
                if(all(P(i,k,:)==0))
                    add = -inf;
                end
                if(add > V_t_max)              
                    V_t_max = add;
                end
            end
            V_t(i) = V_t_max;
        end
        V_1 = [V_1 V_t];
        delta = max(abs(V_1(:,end-1) - V_1(:,end))) - min(abs(V_1(:,end-1) - V_1(:,end)));
    end
end