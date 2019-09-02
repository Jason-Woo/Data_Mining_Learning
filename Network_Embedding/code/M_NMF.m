function U = M_NMF(A,k)
alpha = 5;
beta = 5;
lambda = 1000000000;
m = 4;
[n,~] = size(A);
S1 = A;
S2 = zeros(n,n);
for i=1:n
    for j=1:n
        S2(i,j) = (S1(i,:) * (S1(j,:))') / (norm(S1(i,:)) * norm(S1(j,:)));
    end
end
S = S1+5 .* S2;
e_2 = sum(sum(A));
degree = sum(A);
B = zeros(n,n);
B1 = zeros(n,n);
for i=1:n
    for j=1:n
        B(i,j) = A(i,j) - (degree(i) * degree(j)) / e_2;
        B1(i,j) = (degree(i) * degree(j)) / e_2;
    end
end

M = abs(rand(n,m));  
U = abs(rand(n,m));  
C = abs(rand(k,m));  
H = abs(rand(n,k));  
aim0 = 0;
aim1 = AIM(S,M,U,H,C,B,alpha,beta);
iter = 0;
while(abs(aim0 - aim1) > 0.0000000001)
    M0 = M .* ((S*U)./(M*(U')*U));
    U0 = U .* (((S')*M+alpha*H*C)./(U*((M')*M+alpha*(C')*C)));
    C0 = C .*(((H')*U)./(C*(U')*U));
    delta = (2*beta*(B1*H)) .* (2*beta*(B1*H)) + 16*lambda*(H*(H')*H) .* (2*beta*A*H+2*alpha*U*(C')+(4*lambda-2*alpha)*H);
    H0 = H .* sqrt((-2*beta*B1*H+sqrt(delta))./(8*lambda*H*(H')*H));
    M = M0;
    U = U0;
    C = C0;
    H = H0;
    aim0 = aim1;
    aim1 = AIM(S,M,U,H,C,B,alpha,beta);
    iter = iter +1;
end

end

function aim = AIM(S,M,U,H,C,B,alpha,beta)

aim = (sum_fro(S-M*(U')))^2 + alpha * (sum_fro(H-U*(C')))^2 - beta * trace((H')*B*H);

end

function fro = sum_fro(A)
[m,n] = size(A);
fro = 0;
for i=1:m
    for j=1:n
        fro = fro + A(i,j)^2;
    end
end
fro = sqrt(fro);
end