# Multi-agent-Multi-target-path-planning-in-MDPs

This is the code base for the paper Multi-agent, Multi-target path planning in Markov Decision Processes: https://arxiv.org/pdf/2205.15841.pdf. 

Run the code VI_heur_multi.m for multi-agent planning, and VI_heur_single.m for single-agent planning. The MDP, initial state, and targets should be given as inputs. The functionality of each code is given below. 

1. P_data_stoch.m: generate a random MDP/graph/grid-world. 

2. Partition_Transfers_Swaps.m: Creates a heuristic partition of the target states, given a complete graph and number of agents. 

3. Partitions_opt_k.m: Creates the optimal partition of the target states for the given complete graph and k number of agents. 

4. VI_HT_fn.m: Generates the complete graph that corresponds to a given MDP, initial state and targets, where the weights of the graph are the optimal hitting time between targets and initial state. 

5. VI_heur_fn.m: Computes the suboptimal value function used in the single agent heuristic in Algorithm 1.

6. VI_heur_multi.m: Implements the heuristic path planning algorithm for multiple agents (Algorithm 1), given the MDP, initial state, number of agents, and the assignment of targets to agents.

7. VI_heur_single.m: Implements the heuristic path planning algorithm for single agent (Algorithm 1), given the MDP, initial state, and targets. 

8. VI_opt_fn.m: Implements the optimal policy iteration method and computes the optimal value function for single agent. 

9. heur_impl_fn.m: The agent implements the planned path based on the suboptimal value function, until the agent reaches an unvisited target.
