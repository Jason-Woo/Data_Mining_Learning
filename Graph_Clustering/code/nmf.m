function [U,V] = nmf(A,k)

[n,~] = size(A);
U = abs(rand(n,k));  
V = abs(rand(n,k));  
lse0 = 0;
lse1 = LSE(A-U*V');
iter = 0;
while(abs(lse0 - lse1) > 0.0000000001)
    n_U = U .* (A*V) ./ ((U*V')*V);  
    n_V = V .* (A'*U) ./ ((V*U')*U); 
    U = n_U;
    V = n_V;
    lse0 = lse1;
    lse1 = LSE(A-U*V');
    disp(lse1);
    iter = iter +1;
end

[~,index]=max(V,[],2);
dlmwrite('E:\Jason\ังส๕\Data Mining\Lab\Lab2\code\result.txt',index);
disp(iter);

end

function lse = LSE(M)

lse = 0;
[n,~] = size(M);
for i=1:n
    for j=1:n
        lse = lse + M(i,j);
    end
end

end