% dispSVDnb(typg,CHP1,CHP2,N,expvar) Display SVDs
%
% => Display SVDs in white and black
% typg is plot type you want (see plotm.m).
% CHP contains SVDs and is of the form: CHP(iSVD,X,Y)
% N is number of SVD to display
% expvar contained the explained variance of SVDs
%
% Rq: it''s using the plotmNB function
%
%================================================================

%  Guillaume MAZE - LPO/LMD - December 2004
%  gmaze@univ-brest.fr

function [] = dispSVDnb(typg,CHP1,CHP2,N,expvar)

global ffile1 ffile2
% --------------------------------------------------
% AFFICHAGE
% --------------------------------------------------
if (typg~=1)&(typg~=7)&(typg~=3)&(typg~=6)
   CBAROR='horiz';
   iw=N;jw=2;
else
   CBAROR='vert';
   iw=N;jw=2;
end


f1=figur;clf; hold on
for iN=1:N
     subplot(iw,jw,2*iN-1)
C = squeeze(real(CHP1(iN,:,:)));
%[cs,h]=plotmNB(C,typg); whos cs h
[cs,h]=plotmNB(C,typg);
%clabel(cs,h,'labelspacing',100);
titre1=strcat('SVD',num2str(iN),':',ffile1);
titre2=strcat(' (',num2str(expvar(iN)),'%)');
title(strcat(titre1,titre2));

     subplot(iw,jw,2*iN)
C = squeeze(real(CHP2(iN,:,:)));
[cs,h]=plotmNB(C,typg);
%clabel(cs,h,'labelspacing',100);
titre1=strcat('SVD',num2str(iN),':',ffile2);
titre2=strcat(' (',num2str(expvar(iN)),'%)');
title(strcat(titre1,titre2));

end %for iN




%titre1=strcat(num2str(lat(il1)),'/',num2str(lat(il2)));
titre2=strcat(' Total explained variance:',num2str(sum(expvar)),'%');
%titre3=strcat('(Analysed files: ',ffile1,' vs ',ffile2,')');
%titre=strvcat(titre2,titre3);
suptitle(titre2);

set(f1,'Position',[378 39 313 647]); % Laptop screen
%set(f1,'Position',[369 55 316 899]); % Desktop screen
%set(f1,'Name',strcat(ffile,'<>',titre1,'<> REAL PART'));
