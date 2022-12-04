clear all;
close all;
clc;
%n = 10;
n_S = 80;
n_A = 8;
n_goto_s = 5;
P = zeros(n_S,n_A,n_S); % Transition Probabilities
for i=1:n_S    
    for j=1:n_A
        goto_s = randperm(n_S,n_goto_s);
        r = rand(n_goto_s,1);
        r1 = r / sum(r);
        for k=1:n_goto_s
            P(i,j,goto_s(k)) = r1(k); % Transition dynamics of the MDP
        end
    end
end
%% Grid World
% i=n+2;
% prob = 1;
% k = 5;
% r_prob = [k*prob, prob];
% while(i<(n-1)*n)
%     for a=1:n_A
%         r = repmat(r_prob(2),n_A,1);
%         r(a) = r_prob(1);
%         r=r/sum(r);
%         neigh_s = [-1,-n,1,n];
%         for c_s=1:length(neigh_s)
%             P(i,a,i+neigh_s(c_s)) = r(c_s);
%         end
%     end
%     if(mod(i,n)==n-1)
%         i=i+3;
%     else
%         i=i+1;
%     end
% end
% i=1;
% while(i<=n_S)
%     if(mod(i,n)==1 && i+n<=n_S && i-n>0)        
%         neigh_s = [-n,1,n];
%         neigh_a = [2,3,4];        
%     elseif(mod(i,n)==0 && i+n<=n_S && i-n>0)
%         neigh_s = [-1,-n,n];
%         neigh_a = [1,2,4];
%     elseif(mod(i,n)>1 && i-n<0)
%         neigh_s = [-1,1,n];
%         neigh_a = [1,3,4];
%     elseif(mod(i,n)>1 && i+n>n_S)
%         neigh_s = [-1,-n,1];
%         neigh_a = [1,2,3];
%     elseif(i==1)        
%         neigh_s = [1,n];
%         neigh_a = [3,4];
%     elseif(i==n)
%         neigh_s = [-1,n];
%         neigh_a = [1,4];
%     elseif(i==n*(n-1)+1)        
%         neigh_s = [-n,1];
%         neigh_a = [2,3];
%     elseif(i==n^2)
%         neigh_s = [-1,-n];
%         neigh_a = [1,2];
%     end
%     for a=1:length(neigh_a)
%         r = repmat(r_prob(2),length(neigh_a),1);
%         r(a) = r_prob(1);
%         r=r/sum(r);
%         for c_s=1:length(neigh_s)
%             P(i,a,i+neigh_s(c_s)) = r(c_s);
%         end
%     end
%     if(mod(i,n)==1 && i~=1 && i~=(n-1)*n+1)
%         i=i+n-1;
%     else
%         i=i+1;
%     end
% end