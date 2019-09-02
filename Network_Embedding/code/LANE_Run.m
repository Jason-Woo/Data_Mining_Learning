function ac  = LANE_Run(d,G,A,Y)
k = 50;
delta = 60;
[n,~] = size(G);
n_train = n*0.6;
G_train = G(1:n_train,1:n_train);
A_train = A(1:n_train,:);
A_test = A(n_train+1:n,:);
Y_train = Y(1:n_train);
Y_test_real = Y(n_train+1:n);
H_train = LANE(d,G_train,A_train,Y_train);
G1 = G(1:n_train,:);
G2 = G(n_train+1:n,:);
H_test = G2*pinv(pinv(H_train)*G1)+delta*A_test*pinv(pinv(H_train)*A_train);
Y_test = KNN(H_train,Y_train,H_test,k);
ac = classificationACC(Y_test_real,Y_test);
end
