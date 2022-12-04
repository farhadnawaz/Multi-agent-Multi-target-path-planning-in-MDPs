function [curr_subset,states_v,curr_state] = heur_impl_fn(P,V_t,curr_subset,curr_state,gamma,states_v)
    n_S = size(P,1); % Number of states
    n=sqrt(n_S);
    n_A = size(P,2); % Number of actions
    l_curr_subset = length(curr_subset);    
    while(l_curr_subset == length(curr_subset))        
        to_actions = find(any(P(curr_state,1:n_A,:),3)>0);
        add_a_all = zeros(n_A,1);
        for a=1:n_A
            add_a=0;
            to_states = find(P(curr_state,a,:)>0);
            for r=1:n_S
                %r = to_states(c_r);
                if(any(r==curr_subset))
                    reward_c = -(length(curr_subset)-1);
                else
                    reward_c = -(length(curr_subset));
                end
                add_a = add_a + P(curr_state,a,r)*(reward_c + gamma*V_t(r));
            end
            if(all(P(curr_state,a,:)==0))
                add_a = -100000;
            end
            add_a_all(a) = add_a;
        end
        [V,I] = max(add_a_all);        
        prob = P(curr_state,I,:);
        curr_state = randsrc(1,1,[[1:n_S]; prob(:)']);
        ind_r = find(curr_subset == curr_state);
        curr_subset(ind_r) = [];
        states_v = [states_v curr_state];            
    end
    l_curr_subset = length(curr_subset);
end