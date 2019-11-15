% [CONF,CONFLEVEL] = CONFEXPVAR(A,B,nsamples,N)
%
% A(Ntime,MmapA), B(Ntime,MmapB) are time-synchronised 
%  and time-constant sampled data sets.
% MmapA and MmapB are not required to be identic.
% nsamples is what it is (given at least a 100 factor 
%  will be nice to provide a trustable result).
% N is the max number of eigen vectors to be computed.
%
% CONFEXPVAR Compute the confidence limits for significant 
% eigenvalues in a SVD decomposition of data sets A and B.
% This generates nsamples random-in-time data set of B 
% and finds percent variance of N SVD modes computed 
% with data set A. 
% This is exactly the monte-carlo method but based on 
% a random-in-time only data set. This considerably 
% decrease the cputime.
%
% 2005/11. Guillaume Maze
% gmaze@univ-brest.fr

function [conf,lev,dsumCF] = confexpvar(A,B,nsamples,N)
  
% We try to find the p confidence levels of expvar.
% Keeping as it is the first field, we determine SVDs with the 
% second field arbitrary sorted in time.
% 1st field : A
% 2nd field : B

n = size(A,1);
if size(B,1)~=n
  disp('Error, A and B required to be of the same time length !');
  disp('See help confexpvar');
  stop
end
  
timer = cputime;
disp('Compute confidence levels of explained variance...');
for iR=1:nsamples
    if iR==fix(nsamples/2)
       disp(strcat(num2str(iR),'/',num2str(nsamples)));
       disp(['Please, still wait a little time ... around '...
             num2str(cputime-timer) ' cpu seconds']);
    end
    
    % Generate a random new time index:
    t=randperm(n);
    % and sort the 2nd field with it:
    Br=B(t,:);

    % Now we compute SVD and expvar:
    C=A'*Br;
    [U,Lambda,V] = svds(C,N);
    L2 = Lambda.^2;   dsum = diag(L2)/trace(L2);
    L  = Lambda;      dsum =  diag(L)/trace(L);
    dsumCF(iR,:) = dsum';

end %for iR
disp(['OK, this took ' num2str(cputime-timer) ' cpu seconds'])

% CF(nsamples,N) contains expvar of the nsamples realisations:
CF=sort(dsumCF,1);

% Now we just take the decadal significant levels:
lev=fix([nsamples/10:nsamples/10:nsamples]);
conf = CF(lev,:) ; 
lev=(lev./max(lev))*100;
