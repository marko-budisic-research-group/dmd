function test_svd(X)
% TEST_SVD Visualize differences in SVD for various algorithms

[U,S,V] = svd_lapack(X,0,'gesvd');
plot(diag(S))
semilogy(diag(S),'o','DisplayName','LAPACK gesvd');
hold on;
[U,S2,V] = svd(X,0);
hold on
semilogy(diag(S2),'x','DisplayName','MATLAB SVD0 full')
S3 = svd(X,0);
semilogy(S3,'s','DisplayName','MATLAB SVD0 values')
[U,S4,V] = svd_lapack(X,0,'gesdd');
semilogy(diag(S4),'*','DisplayName','LAPACK gesdd')
hold off;
legend
