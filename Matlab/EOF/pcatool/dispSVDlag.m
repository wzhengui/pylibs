% dispSVDlag(typg,CHP1,CHP2,N,expvar,LAGS,dt) Display SVDs lagged
%
% => Display SVDs lagged
% typg is plot type you want (see plotm.m).
% CHP contains SVDs and is of the form: CHP(iSVD,X,Y)
% N is number of SVD to display
% expvar contained the explained variance of SVDs
%
% Rq: it''s using the plotm function
%
%================================================================

%  Guillaume MAZE - LPO/LMD - December 2004
%  gmaze@univ-brest.fr

   function [] = dispSVD(typg,CHP1,CHP2,N,expvar,LAGS,dt)

load mapanom2
global ffile1 ffile2
% --------------------------------------------------
% AFFICHAGE
% --------------------------------------------------
if (typg~=1)&(typg~=7)&(typg~=3)&(typg~=6)
   CBAROR='horiz';
else
   CBAROR='vert';
end

   CHP1=squeeze(CHP1(:,N,:,:));
   CHP2=squeeze(CHP2(:,N,:,:));

Nlag=size(CHP1,1);

iw=round(Nlag/dt);jw=2;


f1=figur;clf; hold on
ilag=0;
for lag=1:dt:Nlag
  ilag=ilag+1;

    subplot(iw,jw,2*ilag-1)
C = squeeze(real(CHP1(lag,:,:)));
plotm(C,typg);caxis([-abs(xtrm(C)) abs(xtrm(C))]);
if(ilag==1);cx=caxis;end;caxis(cx);
colormap(mycolormap(mapanom,21));
c=colorbar(CBAROR);
titre1=strcat('lag:',num2str(LAGS(lag)),'y: ',ffile1);
title(titre1);

     subplot(iw,jw,2*ilag)
C = squeeze(real(CHP2(lag,:,:)));
plotm(C,typg);caxis([-abs(xtrm(C)) abs(xtrm(C))]);
if(ilag==1);cx=caxis;end;caxis(cx);
colormap(mycolormap(mapanom,21));
c=colorbar(CBAROR);
titre1=strcat('lag:',num2str(LAGS(lag)),'y: ',ffile2);
title(titre1);

end %for iN

suptitle(strcat('Lagged SVD, mode:',num2str(N)));
