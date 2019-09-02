function ac  = M_NMF_Run(A,L)
k = 80;
[n,~] = size(A);
n_train = n*0.6;
A_train = A(1:n_train,1:n_train);
A_test = A(n_train+1:n,n_train+1:n);
L_train = L(1:n_train);
L_test_real = L(n_train+1:n);
U_train = M_NMF(A_train,k);
U_test = M_NMF(A_test,k);
L_test = KNN(U_train,L_train,U_test,k);
ac = classificationACC(L_test_real,L_test);
end
