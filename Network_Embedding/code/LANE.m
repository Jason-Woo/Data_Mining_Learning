function H = LANE(d,G,A,Y)
alpha1 = 100;
alpha2 = 10;
epsion = 0.0001;
[n,~] = size(G);
SG = zeros(n,n);
for i=1:n
    for j=1:n
        if(norm(G(i,:))*norm(G(j,:))==0)
            SG(i,j) = 0;
        else
            SG(i,j) = (dot(G(i,:), G(j,:)))/(norm(G(i,:))*norm(G(j,:)));
        end
    end
end
SA = zeros(n,n);
for i=1:n
    for j=1:n
        if(norm(A(i,:))*norm(A(j,:))==0)
            SA(i,j) = 0;
        else
            SA(i,j) = (dot(A(i,:), A(j,:)))/(norm(A(i,:))*norm(A(j,:)));
        end
    end
end
SY = zeros(n,n);
YY = Y*(Y');
for i=1:n
    for j=1:n
        if(norm(YY(i,:))*norm(YY(j,:))==0)
            SY(i,j) = 0;
        else
            SY(i,j) = (dot(YY(i,:), YY(j,:)))/(norm(YY(i,:))*norm(YY(j,:)));
        end
     end
end
DG = zeros(n,n);
for i=1:n
    DG(i,i) = sum(SG(i,:))+0.1;
end
DA = zeros(n,n);
for i=1:n
    DA(i,i) = sum(SA(i,:))+0.1;
end
DY = zeros(n,n);
for i=1:n
    DY(i,i) = sum(SY(i,:))+0.1;
end
LG = (DG ^(-1/2)) * SG * (DG ^(-1/2));
LA = (DA ^(-1/2)) * SA * (DA ^(-1/2));
LY = (DY ^(-1/2)) * SY * (DY ^(-1/2));

UG = zeros(n,d);
UA = zeros(n,d);
UY = zeros(n,d);
H = zeros(n,d);
iter=0;
J1 = 1000;
flag = 1;
while(flag == 1)  
    M0 = LG + alpha1*UA*(UA') + alpha2*UY*(UY') + H*(H');
    [UG0,~] = eigs(M0,d);
    M1 = alpha1*LA + alpha1*UG*(UG') + H*(H');
    [UA0,~] = eigs(M1,d);
    M2 = alpha2*LY + alpha2*UG*(UG') + H*(H');
    [UY0,~] = eigs(M2,d);
    M3 = UG*(UG') + UA*(UA') + UY*(UY');
    [H0,~] = eigs(M3,d);
    UG = UG0;
    UA = UA0;
    UY = UY0;
    H = H0;
    J0 = J1;
    J1 = calc_J(alpha1,alpha2,UG,UA,UY,LG,LA,LY,H);
    iter=iter+1;
    if(abs(J1-J0) < epsion)
        flag = 0;
    end
end
end

function J = calc_J(alpha1,alpha2,UG,UA,UY,LG,LA,LY,H)
JG = trace(UG'*LG*UG);
JA = trace(UA'*LA*UA);
JY = trace(UY'*LY*UY);
rou1 = trace(UA'*UG*(UG')*UA);
rou2 = trace(UG'*H*(H')*UG);
rou3 = trace(UA'*H*(H')*UA);
rou4 = trace(UY'*H*(H')*UY);
Jcorr = rou2 + rou3 + rou4;

J=JG+alpha1*JA+alpha1*rou1+alpha2*JY+Jcorr;
end

