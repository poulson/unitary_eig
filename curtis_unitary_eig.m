function [X, d, XDefl, XBlock, dDefl, dBlock] = curtis_unitary_eig(Q)

[m, n] = size(Q);

% If w is an eigenvector of the Hermitian part of Q and the projection of
% Qw into the orthogonal complement of w has norm less than the following
% value, then w is deemed an eigenvector of Q.
% NOTE: It's not yet clear how to best choose this value or the degree to
% which carefully forming/solving the 2x2 subproblems is critical.
kDeflationTol = 100 * n * eps^0.8;

% When selectively orthogonalizing members of the overcomplete set of
% eigenvectors, only eigenpairs whose eigenvalues have a distance less than
% or equal to the following value will interact.
kMaxEigValDist = eps^0.5;

% Use every other sorted eigenvalue?
kTrivialSubsample = 1;

% If a projected eigenvector has norm less than this value, it is skipped.
kResidThreshold = 0.9;

% The number of stages of subspace iteration to perform. It seems that,
% typically, the only added value is improved eigenvector orthogonality,
% and so any iteration is not worthwhile.
kNumSubspaceIter = 0;

% Turn this to one if you would like timing information printed.
kTime = 1;

% Increase this value for more debugging information.
kDebug = 2;

% Compute the eigenvectors of twice the Hermitian part of Q.
tic;
S = Q + Q';
[U, ~] = eig(S);
eighTime = toc;
if kTime,
  printf("The initial Hermitian eigensolver took %f seconds.\n", eighTime);
end

% Form an overcomplete set of eigenvectors.
tic;
XDefl = zeros(n, n) + i*zeros(n, n);
dDefl = zeros(n, 1);
QU = Q * U;
numDeflated = 0;
numDelayedPairs = 0;
maxDeflatedResidual = 0;
minUndeflatedResidual = 1;
maxUndeflatedResidual = 0;
W = zeros(n, 2*n);
for j=1:n,
  B = [U(:, j), QU(:, j)];
  % NOTE: If U(:, j) is an eigenvector, then we do not want to form a 
  % 2x2 Rayleigh quotient with a basis for the column space, as the
  % second column vector of the thin QR's Q will be outside of the 
  % target invariant subspace.
  [QB, RB] = qr(B, 0);
  residual = abs(RB(2,2));
  if residual > kDeflationTol,
    W(:,2*numDelayedPairs+1) = QB(:,1);
    W(:,2*numDelayedPairs+2) = QB(:,2);
    numDelayedPairs = numDelayedPairs + 1;
    minUndeflatedResidual = min(minUndeflatedResidual,residual);
    maxUndeflatedResidual = max(maxUndeflatedResidual,residual);
  else,
    XDefl(:,numDeflated+1) = B(:,1);
    % NOTE: We can reuse an above multiplication by Q since QB(:, 1) = U(:, j).
    dDefl(numDeflated+1) = B(:,1)' * B(:,2);
    if kDebug >= 3,
      printf("  Deflated for eigenvalue %e + %ei due to a residual of %e\n", ...
        dDefl(numDeflated+1), dDefl(numDeflated+1)/1i, residual)
    end
    numDeflated = numDeflated + 1;
    maxDeflatedResidual = max(maxDeflatedResidual,residual);
  end
end
if kDebug >= 2,
  printf("  Maximum deflated residual: %e\n", maxDeflatedResidual);
  printf("  Minimum undeflated residual: %e\n", minUndeflatedResidual);
  printf("  Maximum undeflated residual: %e\n", maxUndeflatedResidual);
end
XDefl = XDefl(:,1:numDeflated);
dDefl = dDefl(1:numDeflated);
W = W(:,1:2*numDelayedPairs);
QW = Q*W;
XBlock = zeros(n, 2*numDelayedPairs) + i*zeros(n, 2*numDelayedPairs);
dBlock = zeros(2*numDelayedPairs, 1) + i*zeros(2*numDelayedPairs, 1);
for jPair=1:numDelayedPairs,
  pairInd = [2*jPair-1, 2*jPair];
  % NOTE: We could clearly rearrange our computation to batch this operation.
  Z = W(:,pairInd)' * QW(:,pairInd);
  % Force Z to be unitary.
  [QZ, RZ] = qr(Z);
  phases = [RZ(1,1)/abs(RZ(1,1)), RZ(2,2)/abs(RZ(2,2))];
  Z = QZ * diag(phases);
  [UZ, LambdaZ] = eig(complex(Z));
  XBlock(:, pairInd) = W(:, pairInd) * UZ;
  dBlock(pairInd) = diag(LambdaZ);
end
formXOverTime = toc;
if kTime,
  printf("Deflated %d eigenvectors and delayed %d in %f seconds.\n", ...
    numDeflated, 2*numDelayedPairs, formXOverTime);
end

if kDebug >= 2,
  crossNorm = norm(XDefl'*XBlock, 'fro');
  printf(" || XDefl' XBlock ||_2 = %e\n", crossNorm);
  XDeflOrthogResid = norm(XDefl'*XDefl - eye(numDeflated), 'fro');
  printf(" || XDefl' XDefl - I ||_2 = %e\n", XDeflOrthogResid);
  XBlockOrthogResid = norm(XBlock'*XBlock - eye(2*numDelayedPairs), 'fro');
  printf(" || XBlock' XBlock - I ||_2 = %e\n", XBlockOrthogResid);
end

% Cyclically sort the eigenvalues on the unit circle.
tic;
[~, sortInd] = sort(arg(dBlock));
dBlock = dBlock(sortInd);
XBlock = XBlock(:,sortInd);
sortTime = toc;
if kTime,
  printf("Cyclically sorting took %f seconds.\n", sortTime);
end

% Find a subspace of spanning eigenvectors
X = zeros(n, n) + i*zeros(n, n);
d = zeros(n, 1) + i*zeros(n, 1);
tic;
if kTrivialSubsample,
  X(:,1:numDelayedPairs) = XBlock(:,1:2:end);
  d(1:numDelayedPairs) = dBlock(1:2:end);
  numSaved = numDelayedPairs;
else,
  maxResidNorm = 0.0;
  maxResidIndex = -1;
  numBackwardOrthog = 0;
  numForwardOrthog = 0;
  numSaved = 0;
  for j=1:2*numDelayedPairs,
    x = XBlock(:,j);
    xResid = x;
  
    % Get the list of eigenvectors just before this one that must be
    % orthogonalized against.
    backwardLowerBound = numSaved + 1;
    for jSaved=numSaved:-1:1,
      if abs(d(jSaved)-dBlock(j)) < kMaxEigValDist,
        backwardLowerBound = jSaved;
      else,
        break;
      end
    end
    if backwardLowerBound <= numSaved,
      XBack = X(:,backwardLowerBound:numSaved);
      numBackwardOrthog = numBackwardOrthog + numSaved - backwardLowerBound + 1;
      xResid = xResid - XBack * (XBack' * xResid);
    end

    % Get the list of eigenvectors on the other end of the circle that we
    % need to orthogonalize against.
    forwardUpperBound = 0;
    for jSaved=1:backwardLowerBound-1,
      if abs(d(jSaved)-dBlock(j)) < kMaxEigValDist,
        forwardUpperBound = jSaved; 
      else,
        break;
      end
    end
    if forwardUpperBound >= 1,
      XBeg = X(:,1:forwardUpperBound);
      numForwardOrthog = numForwardOrthog + forwardUpperBound;
      xResid = xResid - XBeg * (XBeg' * xResid);
    end

    residNorm = norm(xResid);
    if residNorm > kResidThreshold,
      X(:,numSaved+1) = xResid / residNorm;
      % TODO(jpoulson): Recompute Rayleigh quotient?
      d(numSaved+1) = dBlock(j);
      numSaved = numSaved + 1;
    else,
      if residNorm > maxResidNorm,
        maxResidNorm = residNorm;
        maxResidIndex = j;
      end
    end
    if numSaved+numDeflated == n,
      break;
    end
  end
  if kDebug >= 2,
    printf("  Total backwards orthogs: %d\n", numBackwardOrthog);
    printf("  Total forwards orthogs: %d\n", numForwardOrthog);
    printf("  Maximum residual norm was %e\n", maxResidNorm);
    if maxResidIndex >= 0,
      printf("  Maximum residual index was %d\n", maxResidIndex);
      printf("  Maximum residual eigenvalue was %e + %ei\n", ...
        dBlock(maxResidIndex), dBlock(maxResidIndex)/1i);
    end
  end
end
% Pack the deflated eigenpairs into the back of X and d.
X(:,n-numDeflated+1:n) = XDefl;
d(n-numDeflated+1:n) = dDefl;
% Complain if we haven't stored enough eigenpairs.
if numSaved+numDeflated ~= n,
  printf("Only saved %d and deflated %d vectors.\n", numSaved, numDeflated);
  return
end
customPivTime = toc;
if kTime,
  printf("Selecting a subset of eigenvectors took %f seconds.\n", ...
    customPivTime);
end

% Measure the initial error
if kDebug >= 1,
  residNorm = norm(Q*X - X*diag(d), 'fro');
  orthogResidNorm = norm(X * X' - eye(n,n), 'fro');
  printf("|| Q X - X D ||_2 = %e\n", residNorm);
  printf("|| X X' - I ||_2 = %e\n", orthogResidNorm);
end

% Perform naive subspace iteration.
for iter=1:kNumSubspaceIter,
  tic;
  [X, R] = qr(Q*X);
  B = Q*X;
  for j=1:n,
    d(j) = X(:,j)'*B(:,j);
  end
  subspaceTime = toc;
  if kTime,
    printf("Subspace iteration took %f seconds.\n", subspaceTime);
  end
  if kDebug >= 1,
    residNorm = norm(Q*X - X*diag(d), 'fro');
    orthogResidNorm = norm(X * X' - eye(n,n), 'fro');
    printf("|| Q X - X D ||_2 = %e\n", residNorm);
    printf("|| X X' - I ||_2 = %e\n", orthogResidNorm);
  end
end

return
