function [value_opt_HT, a_opt_HT, W] = VI_HT_fn(P, init_state, targets)
    n_S = size(P,1); % Number of states
    n_A = size(P,2); % Number of actions
    % weights for hitting time graph
    W = zeros(n_S);
    nu = unique([init_state; double(targets)]);
    n_nu = length(nu);
    value_opt = repmat(-1000, n_S, n_S);
    a_opt = zeros(n_S, n_S);
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
            for curr_goto_s_i = 1:n_nu
                curr_goto_s = nu(curr_goto_s_i);
                curr_goto_s(curr_goto_s==curr_s)=[];
                if(~isempty(curr_goto_s))
                    P_curr_s = reshape(P(curr_s,:,:),n_A,n_S);
                    to_actions = find(sum(P_curr_s,2)>0);
                    V_t_max = -inf;
                    for c_a=1:length(to_actions)
                        a = to_actions(c_a);
                        to_states = find(P(curr_s,a,:)>0);
                        add = 0;
                        for c_to_s = 1:length(to_states)
                            to_s = to_states(c_to_s);
                            goto_s = curr_goto_s;
                            goto_s(goto_s==to_s)=[];
                            if(~isempty(goto_s))
                                add = add + P(curr_s,a,to_s)*value_opt(to_s, goto_s);
                            end
                        end
                        if(-1+gamma*add > V_t_max)
                            a_max = a;               
                            V_t_max = -1 + gamma*add;
                        end
                     end
                     value_opt(curr_s, curr_goto_s) = V_t_max;
                     W(curr_s,curr_goto_s) = -value_opt(curr_s, curr_goto_s);
                     a_opt(curr_s, curr_goto_s) = a_max;
                else
                    value_opt(curr_s, curr_goto_s) = 0;
                    W(curr_s,curr_goto_s) = -value_opt(curr_s, curr_goto_s);
                end
            end
        end
        if(all(a_opt==a_opt_1))
            sf=1;
        end
    end
    value_opt_HT = value_opt;
    a_opt_HT = a_opt;
    toc;
end