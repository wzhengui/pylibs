% [SVD1,SVD2,PC1,PC2,EXPVAR,Lambda] = CALSVD2(A,B,N) Compute SVDs
%
% Ref: H. Bjornson and S.A. Venegas: "A manual for EOF and SVD - 
%      Analyses of climatic Data" 1997
%================================================================
%
%  Guillaume MAZE - LPO/LMD - March 2004
%  gmaze@univ-brest.fr


   function [e1,e2,pc1,pc2,expvar,Lambda,dsumCF] = calsvd2(A,B,N);


%================================================================
% Ref: H. Bjornson and S.A. Venegas: "A manual for EOF and SVD - 
% Analyses of climatic Data" 1997 => p18

% Assume that A is (time*map) matrix
[n p]=size(A);

% Remove the mean of each column (ie the time mean in each station records)
S=detrend(A,'constant');
P=detrend(B,'constant');

% Form the covariance matrix:
C=S'*P;

% Find eigenvectors and singular values
[U,Lambda,V] = svds(C,N);

% PC
a=S*U;
b=P*V;

% Make them clear for output
for iN=1:N
    e1(iN,:) = squeeze( U(:,iN) )';
   pc1(iN,:) = squeeze( a(:,iN) )';
    e2(iN,:) = squeeze( V(:,iN) )';
   pc2(iN,:) = squeeze( b(:,iN) )';
end

% Amount of variance explained a 0.1 pres et en %
L2=Lambda.^2;
dsum=diag(L2)/trace(L2);
for iN=1:N
   expvar(iN)=fix( ( dsum(iN)*100/sum(dsum) )*10 ) /10;
end

