function [c] = kdivide(adj,K)

n = size(adj,1);
W = adj;
%W = W - diag(diag(W));
D =diag(sum(W'));
 
L = D^(-.5)*W*D^(-.5);
 
[eigVectors,eigValues] = eig(L);
[~, ind] = sort(diag(eigValues), 'ascend');
nEigVec = eigVectors(:,ind(1:K));

U=zeros(size(nEigVec,1),K);
for i=1:size(nEigVec,1)
    tt = sqrt(sum(nEigVec(i,:).^2));    
    if tt~=0
        U(i,:) = nEigVec(i,:) ./ tt; 
    end
end

mat = U(:,1:K);
mat_tem = mat;

%K-MEANS
sample_num = size(mat, 1);       % 样本数量
sample_dimension = size(mat, 2); % 每个样本特征维度

% 选距离最远的点做初始类中心
clusters = zeros(K, sample_dimension);
clusters(1,:) = mat_tem(1,:);
for i = 2:K
    num = size(mat_tem,1);
    dist = zeros(num,1);%当前点到其他类中心的平均距离
    for j = 1:num
        total = 0;
        for k = 1:K
            total = total + norm(clusters(k,:)-mat_tem(j,:));
        end
        dist(j) = total/j;
    end
    [~,index] = max(dist);
    clusters(i,:) = mat_tem(index,:);
    mat_tem(index,:) = [];%删去避免重复
    clear dist num;
end
% clusters

c = zeros(sample_num, 1); % 每个样本所属簇的编号

PRECISION = 0.00001;

iter = 70; 
for i=1:iter
    % 遍历所有样本数据，确定所属簇。公式1
    for j=1:sample_num
        gg = repmat(mat(j,:), K, 1);
        gg = gg - clusters;   % norm:计算向量模长
        for n=1:K
            tt(:,n) = norm(gg(n,:));
        end
        [~, minIdx] = min(tt);
        c(j) = minIdx;
    end
%      c'
    
    % 遍历所有样本数据，更新簇心。公式2
    convergence = 1;
    for j=1:K
        up = 0;
        down = 0;
        for k=1:sample_num
            up = up + (c(k)==j) * mat(k,:);
            down = down + (c(k)==j);
        end
        if down==0
            new_cluster = clusters(j,:)+0.1;
        else
            new_cluster = up/down;
        end
        delta = clusters(j,:) - new_cluster;
        if (norm(delta) > PRECISION)
            convergence = 0;
        end
        clusters(j,:) = new_cluster;
    end
%     clusters
    if (convergence)
        break;
    end
end
disp('done');