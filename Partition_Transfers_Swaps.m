function set_part_heur = Partition_Transfers_Swaps(W, init_state, init_part, targets, noA, sub_n, sub_nu, value_opt)
    part_all=[];
    tc=0;
    sc=0;
    c=0;
    oc=[];
    n_S = size(W,1);
    curr_state = init_state;
    nu = double(targets);
    nu(nu==init_state)=[];
    targets = nu;
    tic;
    if(strcmp(init_part, 'k-central'))
        %% k-central vertex Initial partition
        G_c = targets(randi([1,length(targets)]));
        d = W(G_c, targets);
        dset = zeros(noA,1);
        for c_ag=1:noA-1
            [max_val, d_max] = max(d);
            d_max_v = targets(d_max);
            G_c = [G_c, d_max_v];
            for i=1:length(targets)
                d(i) = min([d(i),W(d_max_v,targets(i))]);
            end
        end
        r_targets = setdiff(targets,G_c);
        set_part={[init_state, G_c(1)]};
        for i=2:noA
            set_part= [set_part; {[init_state, G_c(i)]}];
        end
        for ci=1:length(r_targets)
            i = r_targets(ci);
            mind = W(G_c(1),i);
            minc = 1;
            for j=1:noA
                if(W(G_c(j),i)<mind)
                    mind = W(G_c(j),i);
                    minc = j;
                end
            end
            sett = cell2mat(set_part(minc));
            set_part(minc) = {[sett, i]};
        end
    elseif(strcmp(init_part,'sequential'))
        %% sequential initial partition
        set_part = {};
        for i=1:noA
            if(i==noA)
                sett = [init_state, targets(i:length(targets))];
                set_part = [set_part; {sett}];
            else
                sett = [init_state, targets(i)];
                set_part = [set_part; {sett}];
            end
        end
    else
        disp('Invalid choice of initial partition :(');
        set_part = [];
        return
    end
    toc;
    %% Trasnfers and Swaps
    pairs = nchoosek([1:noA],2);
    cht = ones(noA,1);
    chs = ones(noA,1);
    cha = ones(noA,1);
    chsi = ones(size(pairs,1),1);
    chti = chsi;
    chai = chti;
    ch=1;
    while(ch)
        c=c+1;
        for i=1:size(pairs,1)
            pair1 = pairs(i,1);
            pair2 = pairs(i,2);            
            %% Transfers
            if (cht(pair1)==1 || cht(pair2)==1 || chti(i)==1)
                W1t = W(:,cell2mat(set_part(pair1)));
                W1 = W1t(cell2mat(set_part(pair1)),:);
                W2t = W(:,cell2mat(set_part(pair2)));
                W2 = W2t(cell2mat(set_part(pair2)),:);
                size1 = sum(W1,'all')/(length(W1));
                size2 = sum(W2,'all')/(length(W2));
                if (size1<size2)
                    pair1 = pairs(i,2);
                    pair2 = pairs(i,1);
                end
                W1t = W(:,cell2mat(set_part(pair1)));
                W1 = W1t(cell2mat(set_part(pair1)),:);
                W2t = W(:,cell2mat(set_part(pair2)));
                W2 = W2t(cell2mat(set_part(pair2)),:);
                size1 = sum(W1,'all')/(length(W1));
                size2 = sum(W2,'all')/(length(W2));
                mins = max(size1,size2);
                minv=-1;
                set1 = [cell2mat(set_part(pair1))]';
                set1(find(set1==init_state))=[];
                set2 = [cell2mat(set_part(pair2))]';
                set2(find(set2==init_state))=[];
                for j=1:length(set1)
                    vj = set1(j);
                    del1w = sum(W(cell2mat(set_part(pair1)),vj),'all') + sum(W(vj,cell2mat(set_part(pair1))),'all');
                    del2w = sum(W(cell2mat(set_part(pair2)),vj),'all') + sum(W(vj,cell2mat(set_part(pair2))),'all');
                    sizet1 = (sum(W1,'all') - del1w)/((length(W1) - 1));
                    sizet2 = (sum(W2,'all') + del2w)/((length(W2)+1));
                    mint = max(sizet1,sizet2);
                    if (mins>mint)
                        mins = mint;
                        minv = vj;
                    end
                end
                if minv>-1
                    tc=tc+1;
                    set_part(pair1) = {setdiff(cell2mat(set_part(pair1)),minv)};
                    set_part(pair2) = {[cell2mat(set_part(pair2)), minv]};
                    cht(pair1)=1;
                    cht(pair2)=1;
                    chti(i)=1;
                else
                    cht(pair1)=0;
                    cht(pair2)=0;
                    chti(i)=0;
                end
            %% Swaps
            elseif (chs(pair1)==1 || chs(pair2)==1 || chsi(i)==1)
                W1t = W(:,cell2mat(set_part(pair1)));
                W1 = W1t(cell2mat(set_part(pair1)),:);
                W2t = W(:,cell2mat(set_part(pair2)));
                W2 = W2t(cell2mat(set_part(pair2)),:);
                size1 = sum(W1,'all')/(length(W1));
                size2 = sum(W2,'all')/(length(W2));
                mins = max(size1,size2);
                minvj=-1;
                set1 = [cell2mat(set_part(pair1))]';
                set1(find(set1==init_state))=[];
                set2 = [cell2mat(set_part(pair2))]';
                set2(find(set2==init_state))=[];
                for j=1:length(set1)
                    vj = set1(j);
                    delw1j = sum(W(cell2mat(set_part(pair1)),vj),'all') + sum(W(vj, cell2mat(set_part(pair1))),'all');
                    delw2j = sum(W(cell2mat(set_part(pair2)),vj),'all') + sum(W(vj, cell2mat(set_part(pair2))),'all');
                    for k=1:length(set2)
                        vk = set2(k);
                        delw1k = sum(W(cell2mat(set_part(pair1)),vk),'all') + sum(W(vk, cell2mat(set_part(pair1))),'all');
                        delw2k = sum(W(cell2mat(set_part(pair2)),vk),'all') + sum(W(vk, cell2mat(set_part(pair2))),'all');
                        sizet1 = (sum(W1,'all')-delw1j+delw1k-W(vj,vk)-W(vk,vj))/(length(W1));
                        sizet2 = (sum(W2,'all')+delw2j-delw2k-W(vj,vk)-W(vk,vj))/(length(W2));
                        mint = max(sizet1,sizet2);
                        if (mins>mint)
                            mins = mint;
                            minvj = vj;
                            minvk = vk;
                        end
                    end
              end
              if minvj>-1
                  sc=sc+1;
                  set_part(pair1) = {setdiff(cell2mat(set_part(pair1)),minvj)};
                  set_part(pair2) = {setdiff(cell2mat(set_part(pair2)),minvk)};
                  set_part(pair1) = {[cell2mat(set_part(pair1)), minvk]};
                  set_part(pair2) = {[cell2mat(set_part(pair2)), minvj]};
                  chs(pair1)=1;
                  chs(pair2)=1;
                  chsi(i)=1;
              else
                  chs(pair1)=0;
                  chs(pair2)=0;
                  chsi(i)=0;
              end
            end
      end
      ch = sum([sum(chs), sum(cht), sum(chti), sum(chsi)]);
    end
    for i=1:noA
        set_part{i} = setdiff(set_part{i},init_state);
        length_i = length(set_part{i});
        if(length_i==1)
            start_search = 1;
        else
            start_search = 1+sum(sub_n([1:length_i-1]));
        end
        end_search = sum(sub_n([1:length_i]));
        search = start_search;
        while(~all(sub_nu{search}==set_part{i}))
            search = search+1;
        end
        set_part{i,2} = search;
        set_part{i,3} = -value_opt(init_state, search);
    end
    set_part{noA+1,3} = max([set_part{:,3}]);
    set_part_heur = set_part;
    toc;
end