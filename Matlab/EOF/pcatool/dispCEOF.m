% dispCEOF(CEOF,EXPVAR,PHI,MOD) Display few CEOFs.
%
% => DISPLAY FEW CEOFs.
% CEOF contains all the CEOFs as nCEOF*X*Y.
% EXPVAR is a matrix with the explained variance of each
%  CEOFs in %. This is just for title.
% PHI is phase when ploting maps.
% MOD contains explicitly the CEOFs to display.
%
% Rq: it's using the plotm function
%
%================================================================

%  Guillaume MAZE - LPO/LMD - June 2004
%  gmaze@univ-brest.fr

function [] = dispCEOF(CEOF,EXPVAR,PHI,MOD)

% Number of EEOF to display     
NMOD = length(MOD);
NANG = length(PHI);

% Open figure and first guest
figur;clf;hold on
typg=8;
width  = .9/NMOD;
height = .9/NANG;
dleft = (.95-width*NMOD)/4;
left  = dleft;
load mapanom2

% Let'go :

% We select a CEOF:
for imod = 1 : NMOD
    mod = MOD(imod);
    chp=squeeze(CEOF(mod,:,:));
    Rchp=real(chp);
    Ichp=-imag(chp);

    % 
    for ang=1:NANG
      C = cos(PHI(ang))*Rchp + sin(PHI(ang))*Ichp;
      m(ang)=abs(xtrm(C));
    end
    m=max(m);
    
    % and plot maps
    for ang=1:NANG
      C = cos(PHI(ang))*Rchp + sin(PHI(ang))*Ichp;     
      bottom = .95-height*ang;
      subplot('Position',[left bottom width height]);
      plotm(C,typg);
      m_text(10,-70,strcat('\phi=',num2str(PHI(ang)*180/pi),'^o'));
%      ylabel(strcat('\phi=',num2str(PHI(ang)*180/pi),'^o'),'Rotation',0);
      caxis([-m m]);
      colormap(mapanom);
      if(ang==1)
         titre=strcat('CEOF',num2str(mod),'(',num2str(EXPVAR(mod)),'%)');
         title(titre);
      end;
    end
    left = left + width + dleft;

end
