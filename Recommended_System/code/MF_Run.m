function result = MF_Run(train,test)
answer = test(:,3);
score = MF(train,test);
[test_size,~] = size(test);
aver_answer = sum(sum(answer)) / test_size;
aver_score = sum(sum(score)) / test_size;
TP = 0;
FN = 0;
FP = 0;
TN = 0;
for i=1:test_size
    if(answer(i,1)>=aver_answer && score(i,1)>=aver_score)
        TP = TP + 1;
    elseif(answer(i,1)>=aver_answer && score(i,1)<aver_score)
        FN = FN + 1;
    elseif (answer(i,1)<aver_answer && score(i,1)>=aver_score)
        FP = FP + 1;
    elseif (answer(i,1)<aver_answer && score(i,1)<aver_score)    
        TN = TN + 1;
    end
end
precision = TP / (TP + FP);
recall = TP / (TP + FN);
f1 = (2 * TP) / (2 * TP + FP + FN);
result = [precision recall f1];
end

function result = MF(train,test)
[m,n] = size(train);
alpha = 1;
k = 20;
P = zeros(m,k);
Q = zeros(k,n);
[e,loss0] = calc_loss(train,P,Q);
loss1 = loss0 + 1;

while(abs(loss1-loss0)>0.001)    
    for i=1:m
        for j=1:n
            for l=1:k
            P(i,k) = P(i,k) + 2 * alpha * e(i,j) * Q(k,j);
            Q(k,j) = Q(k,j) + 2 * alpha * e(i,j) * P(i,k);
            end
        end
    end
    loss0 = loss1;
    [e,loss1] = calc_loss(train,P,Q);
end
R_est = (P*Q);
[test_size,~] = size(test);
result = zeros(test_size,0);
for i=1:test_size
    result(i,1) = R_est(test(i,1),test(i,2));
end
end

function [e,loss] = calc_loss(R,P,Q)
[m,n] = size(R);
[~,k] = size(P);
e = R;
e2 = zeros(m,n);
for i=1:m
    for j=1:n
        for l=1:k
            e(i,j) = e(i,j) - (P(i,l)*Q(l,j));
        end
        e2(i,j) = e(i,j)^2;
    end
end
loss = sum(sum(e2));
end