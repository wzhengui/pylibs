% dispEEOF_v2(CHP,PC,EXPVAR,DT,NLAG,MOD) Display few EEOFs with their
%                                     PC time evolution
%
% => DISPLAY FEW EEOFs with their PC time evolution.
% CHP contains all the EEOFs as EOF*LAG*X*Y.
% PC contains all the PCs.
% EXPVAR is a matrix with the explained variance of each
%  EEOFs in %. This is just for title.
% DT is time step between maps.
% NLAG is the number of LAG to display.
% MOD contains explicitly the EEOFs to display.
%
% Rq: it's using the plotm function
% See also: dispEEOF
%
%================================================================

%  Guillaume MAZE - LPO/LMD - June 2004
%  gmaze@univ-brest.fr

function [] = dispEEOF(CHP,PC,EXPVAR,DT,NLAG,MOD)

% Number of EEOF to display     
NMOD = length(MOD);

% Open figure and first guest
figur;
clf;hold on
typg=8;
width  = .9/NMOD;
height = .9/(NLAG+1);
dleft = (.95-width*NMOD)/4;
left  = dleft;
load mapanom

% Let'go :

% We choose an EEOF
for imod = 1 : NMOD
    mod = MOD(imod);
    
    % and plot maps
    for iplot=1:NLAG+1
      bottom = .95-height*iplot;
      subplot('Position',[left bottom width height]);
      if iplot==1 % Plot the pc's time serie
         plot(PC(imod,:));
         grid on;box on;
         set(gca,'XAxisLocation','top');axis tight;
         titre=strcat('EEOF',num2str(mod),'(',num2str(EXPVAR(mod)),'%)');
         titre=strcat(titre,'; DT=',num2str(DT),' ');
         title(titre);
         if (imod==1),ax1=get(gca,'YLim');else,set(gca,'YLim',ax1);end
      else % Plot the lag maps
         lag=iplot-1;
         C = squeeze(real(CHP(mod,lag,:,:)));
%         C = C./xtrm(C); % Eventually normalise field
         plotm(C,typg);
         if(lag==1),caxis([-abs(xtrm(C)) abs(xtrm(C))]);cx=caxis;end;
         caxis(cx);
         colormap(mapanom);
         %ylabel(strcat('{\bf',num2str((lag-1)*DT),'}'))
         m_text(10,-70,strcat('{\bf',num2str((lag-1)*DT),'}'));
      end %if

    end
    left = left + width + dleft;

    % Adjust plot width via the plot Aspect ratio
    ch=get(gcf,'Children');
    posiPC=get(ch(iplot),'PlotBoxAspectRatio');
    posiMAP=get(ch(iplot-1),'PlotBoxAspectRatio');
    set(ch(iplot),'PlotBoxAspectRatio',posiMAP);

end

