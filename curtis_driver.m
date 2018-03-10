% Test a sequence of random members of O(n, R) or O(n, C).
n = 1000

% Test a random member of O(n).
printf("Testing O(%d).\n", n);
A = randn(n,n);
[Q, R] = qr(A);
[X, d, XDefl, XBlock, dDefl, dBlock] = curtis_unitary_eig(Q);
printf("\n");

% Test a random member of U(n).
printf("Testing U(%d).\n", n);
A = randn(n,n) + randn(n,n)*i;
[Q, R] = qr(A);
[X, d, XDefl, XBlock, dDefl, dBlock] = curtis_unitary_eig(Q);
printf("\n");

% Test a highly clustered member of U(400).
printf("Testing clustered U(400).\n");
A = randn(400,400);
[Q, R] = qr(A);
[X, d] = curtis_unitary_eig(Q*diag([ones(100,1);(ones(100,1)*1i+ones(100,1))/sqrt(2);-ones(200,1)])*Q');
printf("\n");

% Test I.
printf("Testing eye(%d).\n", n);
[X, d] = curtis_unitary_eig(eye(n,n));
printf("\n");
