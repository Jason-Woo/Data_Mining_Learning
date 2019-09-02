function [C] = normalized_cut(adj,target_K,fea)

%adj = cal_adj(adj,fea);

K = 10;
c = kdivide(adj,K);
c_temp = c;
sample_number = size(c,1);

while K>target_K
    Ncut = zeros(K);
    for q = 1:K
        for p = 1:K
            c = c_temp;
            for i = 1:sample_number
                if c(i) == q
                    c(i) = p;
                end
            end     
            cutAB = zeros(sample_number,1);
            assocAV = zeros(sample_number,1);
            for i=1:sample_number
                for j=1:sample_number
                    if adj(i,j) == 1
                        assocAV(i) = assocAV(i) + 1;
                    end
                    if adj(i,j) == 1 && i~=j && c(i) ~= c(j)
                        cutAB(i) = cutAB(i) + 1;
                        
                    end
                end
            end
            for i = 1:sample_number
                if assocAV(i) ~= 0
                    Ncut(q,p) = Ncut(q,p) + cutAB(i)/assocAV(i);
                end
            end 
            clear c;
        end
    end
    
    [a,min_row] = min(Ncut);
    [~,min_col] = min(a);
    
    for i = 1:sample_number
        if c_temp(i) == min_col
             c_temp(i) = min_row(min_col);
        end
    end
    
    K = K - 1;
    clear Ncut; 
end

C = c_temp;