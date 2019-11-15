% [EOFs,PC,EXPVAR] = calEeof(M,N,METHOD,NLAG,DT) Compute EEOF
%
% => Compute N Extended EOF of a M(TIME*MAP) matrix.
% METHOD is the EOF method computing (see caleof.m)
% NLAG is the number of map you want
% DT is the increment between each map (not the time step
% of your data)
%
% The function can filter datas with a Lanczos band-pass 
% filter. Edit the calEeof.m file and change option in it. 
%
% See also: CALEOF
%
% Ref: Weare and Nasstrom 1982
%================================================================

%  Guillaume MAZE - LPO/LMD - March 2004
%  gmaze@univ-brest.fr

    function varargout = calEeof(M,N,METHOD,NLAG,DT)

% ---------------------------------
% Divers
tic
LAG = NLAG*DT;
foll = 1; % Turn this value to 0 if you don t want commentary along computing

% ---------------------------------
% LAG must be >0, LAG=1 correspond to initial position
if LAG<=0
   LAG = 1;
end

% ---------------------------------
% Usely for us, M is a TIME*MAP matrix so:
% We get dimensions
  [n p]=size(M);
% and put M in the correct form for computing
  M = M';

% ---------------------------------
% Eventualy time filtering of data
if 0==1
   if (foll~=0),disp('==> Time filtering...');end
   dt = 1;         % Real time step of your data (diff from DT)
   Fc  = 1/dt/8;   % This the low cut off frequency
   Fc2 = 1/dt/2;   % and this the high cut off frequency
   SIGNAL = M(1,:);
   nj = fix(length(SIGNAL)/10); % Nb of points of the window

   for ipt = 1 : p
       SIGNAL = M(ipt,:)';
       SIGNALF = lanczos(SIGNAL,Fc2,nj);
       SIGNALF = SIGNALF - lanczos(SIGNALF,Fc,nj);
       Y(:,ipt) = SIGNALF;
       if mod(ipt,10)==0 % We display a plot of filtered data
          clf;
          plot((0:n-1),SIGNAL,'k',(0:n-1),SIGNALF,'r');          
          drawnow;
       end
   end
   M = Y';
end, clear Fc Fc2 SIGNAL nj SIGNALF Y ipt, close


% ---------------------------------
% This is the matrix where we concat all submatrices from M.
% Size of each submatrix is : NLAG*p rows vs n-LAG+1 colums 
F = zeros(NLAG*p,n-LAG+1);

% ---------------------------------
% Let's go for F:
if (foll~=0),disp('==> Forming concatenated matrix...');end
if DT==1
% CLASSIC CASE
for ilag = 1 : LAG
       % Extract submatrix ilag from M
       Msub = M(:,ilag:n-LAG+ilag);
       % Concat matrix
       F( (ilag-1)*p+1 : ilag*p , : ) = Msub;
end
else
% DT>1
it=0;
for ilag = 1 : DT : LAG
       % Extract submatrix ilag from M
         Msub = M( : , ilag : n-LAG+ilag);
       % Concat matrix
         it = it + 1;
         F( (it-1)*p+1 : it*p , : ) = Msub;
         % imagesc(Msub);pause % for debugging
end
end

% ---------------------------------
% Compute EOFs by normal way
% (Don't forget that caleof is taking a TIME*MAP argument matrix,
% that's why we have F' !)
if (foll~=0),disp('==> Computing EOFs ...');end
[e,pc,expvar] = caleof(F',N,METHOD);

% ---------------------------------
% e is the matrix with EOFs inside. We should extract each map
% for different lags from it.
% e is NEOF*(LAG*MAP) and we construct 3D matrix CHP as NEOF*LAG*MAP
if (foll~=0),disp('==> Reshaping vectors ...');end
for ilag = 1 : NLAG
    E = e( : , (ilag-1)*p+1 : ilag*p );
    CHP(:,ilag,:) = E;
end


% ---------------------------------
% OUTPUT VARIABLES
switch nargout
  case 1
   varargout(1) = {CHP} ;
  case 2
   varargout(1) = {CHP} ;
   varargout(2) = {pc};
  case 3
   varargout(1) = {CHP} ;
   varargout(2) = {pc};
   varargout(3) = {expvar};
end

% ---------------------------------
if (foll~=0),toc,disp('==> That''s all folks !');end
