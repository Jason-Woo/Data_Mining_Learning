function result = IBCF_Run(train,test)
answer = test(:,3);
score = IBCF(train,test);
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

function score = IBCF(train,test)
[user_num,item_num] = size(train);
R_aver = mean(train,2);
S = zeros(item_num,item_num);
for i=1:item_num
    for j=i+1:item_num
        temp1=0;
        temp2=0;
        temp3=0;
        for k=1:user_num
            temp1 = temp1 + ((train(k,i)-R_aver(k))*(train(k,j)-R_aver(k)));
            temp2 = temp2 + ((train(k,i)-R_aver(k))^2);
            temp3 = temp3 + ((train(k,j)-R_aver(k))^2);
        end
        S(i,j) = temp1/(sqrt(temp2)*sqrt(temp3));
        S(j,i) = S(i,j);
    end  
end

[test_size,~] = size(test);
score = zeros(test_size,1);
for i=i:test_size
    temp1=0;
    temp2=0;
    for j=1:item_num
        temp1 = temp1 + (S(test(i,2),j) * train(test(i,1),j));
        temp2 = temp2 + S(test(i,2),j);
    end
    score(i,1) = temp1 / temp2;
end
end