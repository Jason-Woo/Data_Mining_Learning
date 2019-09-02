function result = Slope_One_Run(train,test)
answer = test(:,3);
score = Slope_One(train,test);
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

function score = Slope_One(train,test)
[user_num,item_num] = size(train);
dev = zeros(item_num,item_num);
for i=1:item_num
    for j=i+1:item_num
        user_cnt=0;
        for k=1:user_num
           if(train(k,i)~=0 && train(k,j)~=0)
               dev(i,j) = dev(i,j)+ (train(k,i) - train(k,j));
               user_cnt = user_cnt + 1;
           end
        end
        if(user_cnt~=0)
            dev(i,j) = dev(i,j) / user_cnt;
            dev(j,i) = dev(i,j);
        else
            dev(i,j) = 0;
            dev(j,i) = 0;
        end
    end
end
[test_num,~] = size(test);
score = zeros(test_num,1);
for i=1:test_num
    item_cnt = 0;
    for j=1:item_num
        if(train(test(i,1),j)~=0)
            item_cnt = item_cnt + 1;
            score(i,1) = score(i,1) + (dev(test(i,2),j) + train(test(i,1),j));
        end
    end
    if(item_cnt~=0)
        score(i,1) = score(i,1) / item_cnt;
    else
        score(i,1) = 0;
    end
end
end