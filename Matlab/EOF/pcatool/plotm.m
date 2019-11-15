%PLOTM  Trace un champ de donnee geographique
% [CS,H]=PLOTM(C,[ZONE]) 
%   Trace le champ C(LONGITUDE,LATITUDE) sur une carte 
%   du monde en projection Mercator par defaut.
%   Le parametre optionnel ZONE peut prendre les valeurs
%     ZONE = 1 : Projection Stereographique centree sur le pole sud
%     ZONE = 2 : Projection Mercator de l'hemisphere sud
%     ZONE = 3 : Projection Mercator du globe
%     ZONE = 4 : Pas de projection (imagesc)
%     ZONE = 5 : Projection Mercator de l'hemisphere sud complet
%     ZONE = 6 : Projection Mercator du globe complet
%     ZONE = 7 : Projection Stereographique centree sur le pole nord
%     ZONE = 8 : Idem 2 mais sans label sur les axes
%     ZONE = 9 : Projection Mercator de l'hemisphere sud aux limites du Fmask
%   Returns contour matrix CS as described in CONTOURC and a vector H
%   of handles to PATCH objects (for use by CLABEL).
%
%   La definition des grilles (longitude,latitude) est celle du
%   modele par defaut, mais elle peut etre modifiee en passant en
%   global les nouvelles variables nommees lon et lat.
%
% NECESSITE LE PACKAGE M_MAP
%
%
% (nov. 2003 -> gmaze@univ-brest.fr)

function varargout = plotm(C,varargin)

if (nargin>2)|(nargin<1)
     help plotm.m
     error('plotm.m : Wrong number or bad parameter(s)')
     return
elseif (nargin==2)
   arg = varargin{:};
   if (arg>14)|(arg<1)
     help plotm.m
     error('plotm.m : Wrong number or bad parameter(s)')
     return
   end
end %if

%--------------------------------------------------
% Variabes
%--------------------------------------------------
[lx,ly]=size(C);
% Si lon et lat ne sont pas definie comme dans le modele
% on prend en global n'importe quel definition
if ((lx~=64)|(ly~=32)) 
   global lon lat
else
   lon= linspace(0.,360.,64);
   lat= [85.76 80.27 74.74 69.21 63.68 58.14 52.61 47.07 41.53 ...
      36.00 30.46 24.92 19.38 13.84 8.31 2.77 -2.77 -8.31 ...
      -13.84 -19.38 -24.92 -30.46 -36.00 -41.53 -47.07 -52.61 ...
      -58.14 -63.68 -69.21 -74.74 -80.27 -85.76];
end

% Le nombre de contour est par defaut a 10, mais on peut le
% modifier en declarant nbcont en global
global nbcont
if isempty(nbcont)==1
   nbcont=10;
end

%-----------------------------------------------
% Projection
%-----------------------------------------------
if (nargin==2)
   arg = varargin{:};
   switch arg(1)
     case 1 % Projection stereo 
       m_proj('stereo','lon',0,'lat',-90,'rad',71);
     case 2 % Projection mercator (HS)
       m_proj('miller','lon',[0 360],'lat',[-80 -20]);
       m_proj('miller','lon',[0 360],'lat',[-60 -20]);
     case 3 % Projection mercator (Tout)  -> Defaut
       m_proj('miller','lon',[0 360],'lat',[-80 80]);
     case 5 % Projection mercator (HS-complet)
       m_proj('miller','lon',[0 360],'lat',[-85 -5]);
     case 6 % Projection mercator (Tout-complet)
       m_proj('miller','lon',[0 360],'lat',[-85 85]);
     case 7 % Projection stereo 
       m_proj('stereo','lon',0,'lat',90,'rad',71);
     case 8 % Projection mercator (HS) sans Label sur les axes
       m_proj('miller','lon',[0 360],'lat',[-80 -20]);
       m_proj('miller','lon',[0 360],'lat',[-75 -30]);
     case 9 % Projection mercator (HS) sans Label sur les axes aux limtes du FMASK
       m_proj('miller','lon',[0 360],'lat',[-69 -25]);
     case 10 % Projection stereo recentree sur le pole sud
       m_proj('stereo','lon',0,'lat',-90,'rad',51);
     case 11 % Indien
       m_proj('lambert','long',[11 125],'lat',[-60 -30]);
     case 12 % Pacifique
       m_proj('lambert','long',[120 290],'lat',[-60 -30]);
     case 13 % Atlantique
       m_proj('lambert','long',[-70 20],'lat',[-60 -30]);
     case 14 % Centré Atlantique        
       m_proj('lambert','long',[-90 90],'lat',[-70 -20],'rectbox','on');
   end
else
% Projection mercator (Tout)  -> Defaut
  m_proj('mercator','lon',[0 360],'lat',[-80 80]);
end

%------------------------------------------------
% Affichage
%------------------------------------------------
load MODELCOAST
if length(lon)==65
  modelcoast=[modelcoast;modelcoast(1,:)];
end

hold on
if (nargin==2)
   arg = varargin{:};
   switch arg(1)
     case 1 % Projection stereo (Par defaut)
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_coast('color',[0 0 0]);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-30 -40 -50 -60 -70],'yticklabels',[],...
             'xtick',[0 60 120 180 -120 -60],'XaxisLocation','top')
     case 2 % Projection mercator (HS)
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_coast('color',[0 0 0]);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-80 -70 -60 -50 -40 -30 -20]);
     case 3 % Projection mercator (Tout)
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_coast('color',[0 0 0]);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-80 -60 -40 -20 0 20 40 60 80]);
     case 4 % Pas de projection
      pcolor(lon,lat,C');
      axis tight;box on
%      axis([lon(1) lon(64) lat(32) lat(1)]);
%      gc=gca;set(gc,'YDir','reverse');
     case 5 % Projection mercator (Hs-complet)
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_coast('color',[0 0 0]);
      COASTYPE='color';COASTCOLOR=[0 0 0];
%      COASTYPE='patch';COASTCOLOR=[0 0.5 0.6];
%      m_coast(COASTYPE,COASTCOLOR);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-80 -60 -40 -20 0]);
     case 6 % Projection mercator (Tout-complet)
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_coast('color',[0 0 0]);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-80 -60 -40 -20 0 20 40 60 80]);
     case 7 % Projection stereo
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_coast('color',[0 0 0]);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[30 40 50 60 70],'yticklabels',[],...
             'xtick',[0 60 120 180 -120 -60],'XaxisLocation','top')
     case 8 % Projection mercator (HS) sans label sur les axes
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_pcolor(lon,lat,C');shading interp
%      m_coast('color',[0 0 0]);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick','','xtick','');
     case 9 % Projection mercator (HS) sans label sur les axes aux limtes du FMASK
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_coast('color',[0 0 0]);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick','','xtick','');
     case 10 % Projection stereo recentree sur le pole Sud
      [cs,h]=m_contourf(lon,lat,C',nbcont);
%      m_coast('color',[0 0 0]);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-40 -50 -60 -70],'yticklabels',[],...
             'xtick',[0 60 120 180 -120 -60],'XaxisLocation','top')
     case 11 % Indien
      [cs,h]=m_contourf(lon,lat,C',nbcont);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-40 -50 -60 -70],'yticklabels',[],...
             'xtick',[30 60 90 120 150 180],'XaxisLocation','top')
     case 12 % Pacifique
      [cs,h]=m_contourf(lon,lat,C',nbcont);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-40 -50 -60 -70],'yticklabels',[],...
             'xtick',[120 150 180 210 240 270 300],'XaxisLocation','top')
     case 13 % Atlantique
      [cs,h]=m_contourf(lon,lat,C',nbcont);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('ytick',[-40 -50 -60 -70],'yticklabels',[],...
             'xtick',[-60 -30 0 30 60],'XaxisLocation','top')
     case 14 % Centré Atlantique
      [cs,h]=m_contourf(lon,lat,C',nbcont);
      m_contour(lon,lat,modelcoast',1,'k');
      m_grid('yticklabels',[],'XaxisLocation','top',...
             'fontsize',6,...
             'xtick',[0 60 120 180 -120 -60],'xticklabels',[180 -120 -60 0 60 120]);
            
   end % switch
else
  %  Projection Mercator (Par defaut)
     [cs,h]=m_contourf(lon,lat,C',nbcont);
%     m_coast('color',[0 0 0]);
     m_contour(lon,lat,modelcoast',1,'k');
     m_grid('ytick',[-80 -60 -40 -20 0 20 40 60 80]);
end



% OUTPUT VARIABLES
switch nargout
  case 1
   varargout(1) = {cs} ;
  case 2
   varargout(1) = {cs} ;
   varargout(2) = {h};
end
