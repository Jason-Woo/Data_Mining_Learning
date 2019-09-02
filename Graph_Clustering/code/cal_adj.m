function [adj_weight] = cal_adj(adj,fea)

sigma = 1;
n = size(adj,1);
W = zeros(n,n);
for i = 1:n
    for j = n:-1:i
        if adj(i,j) == 1
            W(i,j) = exp(-norm(fea(i,:)-fea(j,:))^2/2/sigma);
        end
        W(j,i) = W(i,j);
    end
end

adj_weight = W;