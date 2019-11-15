% dispEEOF(CHP,EXPVAR,DT,NLAG,MOD) Display few EEOFs.
%
% => DISPLAY FEW EEOFs.
% CHP contains all the EEOFs as EOF*LAG*X*Y.
% EXPVAR is a matrix with the explained variance of each
%  EEOFs in %. This is just for title.
% DT is time step between maps.
% NLAG is the number of LAG to display.
% MOD contains explicitly the EEOFs to display.
%
% Rq: it's using the plotm function
%
%================================================================

%  Guillaume MAZE - LPO/LMD - March 2004
%  gmaze@univ-brest.fr

function [] = dispEEOF(CHP,EXPVAR,DT,NLAG,MOD)

% Number of EEOF to display     
NMOD = length(MOD);

% Open figure and first guest
figur;
clf;hold on
typg=8;
width  = .9/NMOD;
height = .9/NLAG;
dleft = (.95-width*NMOD)/4;
left  = dleft;
load mapanom

% Let'go :

% We choose an EEOF
for imod = 1 : NMOD
    mod = MOD(imod);
    
    % and plot maps
    for lag=1:NLAG
      bottom = .95-height*lag;
      subplot('Position',[left bottom width height]);
      C = squeeze(real(CHP(mod,lag,:,:)));
%      C = C./xtrm(C); % Eventually normalise field
      plotm(C,typg);
      if(lag==1),caxis([-abs(xtrm(C)) abs(xtrm(C))]);cx=caxis;end;
      caxis(cx);
      colormap(mapanom);
      if(lag==1)
         titre=strcat('EOF',num2str(mod),'(',num2str(EXPVAR(mod)),'%)');
         titre=strcat(titre,'; DT=',num2str(DT),' ');
         title(titre);
      end;
    end
    left = left + width + dleft;

end
