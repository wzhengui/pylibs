% dispSVD(typg,CHP1,CHP2,N,expvar) Display SVDs
%
% => Display SVDs
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

function [] = dispSVD(typg,CHP1,CHP2,N,expvar)

load mapanom2
global ffile1 ffile2

if length(N)>1
  mod=N;
else
  mod=[1:N];
end  
Nx=length(mod);
% --------------------------------------------------
% AFFICHAGE
% --------------------------------------------------
if (typg~=1)&(typg~=7)&(typg~=3)&(typg~=6)
   CBAROR='horiz';
   iw=Nx;jw=2;
else
   CBAROR='vert';
   iw=Nx;jw=2;
end


f1=figur;clf; hold on
for iN=1:Nx
     subplot(iw,jw,2*iN-1)
C = squeeze(real(CHP1(mod(iN),:,:)));
%plotm(C./xtrm(C),typg);caxis([-abs(xtrm(C)) abs(xtrm(C))]);
plotm(C,typg);
caxis([-abs(xtrm(C)) abs(xtrm(C))]);
%if(iN==1);cx=caxis;end;caxis(cx);
colormap(mycolormap(mapanom,21));
c=colorbar(CBAROR);
titre1=strcat('SVD',num2str(mod(iN)),':',blanks(1),ffile1);
titre2=strcat(' (',num2str(expvar(mod(iN))),'%)');
tt=title(strvcat(strcat(titre1,titre2),char(1)));
set(tt,'fontweight','bold')

     subplot(iw,jw,2*iN)
C = squeeze(real(CHP2(mod(iN),:,:)));
%plotm(C./xtrm(C),typg);caxis([-abs(xtrm(C)) abs(xtrm(C))]);
plotm(C,typg);
caxis([-abs(xtrm(C)) abs(xtrm(C))]);
%if(iN==1);cx=caxis;end;caxis(cx);
colormap(mycolormap(mapanom,21));
c=colorbar(CBAROR);
titre1=strcat('SVD',num2str(mod(iN)),':',blanks(1),ffile2);
titre2=strcat(' (',num2str(expvar(mod(iN))),'%)');
tt=title(strvcat(strcat(titre1,titre2),char(1)));
set(tt,'fontweight','bold')

end %for iN


% $$$ for iN=1:2*N
% $$$   figure(f1);
% $$$   subplot(iw,jw,iN);
% $$$ %  colorbar(CBAROR);
% $$$ end
% $$$ refresh


%titre1=strcat(num2str(lat(il1)),'/',num2str(lat(il2)));
titre2=strcat(' Total explained variance:',num2str(sum(expvar(mod))),'%');
%titre3=strcat('(Analysed files: ',ffile1,' vs ',ffile2,')');
%titre=strvcat(titre2,titre3);
suptitle(titre2);

set(f1,'Position',[378 39 313 647]); % Laptop screen
%set(f1,'Position',[369 55 316 899]); % Desktop screen
%set(f1,'Name',strcat(ffile,'<>',titre1,'<> REAL PART'));
