% dispEOF(typg,CHP,N,expvar) Display EOFs
%
% => Display EOFs
% typg is plot type you want (see plotm.m).
% CHP contains EOFs and is of the form: CHP(EOF,X,Y)
% N is number of EOF to display
% expvar contained the explained variance of EOFs
%
% Rq: it''s using the plotm function
%
%================================================================

%  Guillaume MAZE - LPO/LMD - March 2004
%  gmaze@univ-brest.fr

function [] = dispEOF(typg,CHP,N,expvar)

load mapanom2
global ffile
% --------------------------------------------------
% AFFICHAGE
% --------------------------------------------------
if (typg~=1)&(typg~=7)&(typg~=3)&(typg~=6)

   CBAROR='horiz';
if N>12,NN=12;else,NN=N;end
switch NN
  case 1,tblmagique=[1];                         iw=1;jw=1;
  case 2,tblmagique=[1 2];                       iw=2;jw=1;
  case 3,tblmagique=[1 2 3];                     iw=3;jw=1;
  case 4,tblmagique=[1 2 3 4];                   iw=2;jw=2;
  case 5,tblmagique=[1 2 3 4 5];                 iw=3;jw=2;
  case 6,tblmagique=[1 2 3 4 5 6];               iw=3;jw=2;
  case 7,tblmagique=[1 2 3 4 5 6 7];             iw=4;jw=2;
  case 8,tblmagique=[1 2 3 4 5 6 7 8];           iw=4;jw=2;
  case 9,tblmagique=[1 2 3 4 5 6 7 8 9];         iw=3;jw=3;
 case 10,tblmagique=[1 2 3 4 5 6 7 8 9 10];      iw=4;jw=3;
 case 11,tblmagique=[1 2 3 4 5 6 7 8 9 10 11];   iw=4;jw=3;
 case 12,tblmagique=[1 2 3 4 5 6 7 8 9 10 11 12];iw=4;jw=3;
end
iw=N;jw=1;

else

   CBAROR='vert';
if N>12,NN=12;else,NN=N;end
switch NN
  case 1,tblmagique=[1];                         iw=1;jw=1;
  case 2,tblmagique=[1 2];                       iw=2;jw=1;
  case 3,tblmagique=[1 2 3];                     iw=3;jw=1;
  case 4,tblmagique=[1 2 3 4];                   iw=2;jw=2;
  case 5,tblmagique=[1 2 3 4 5];                 iw=3;jw=2;
  case 6,tblmagique=[1 2 3 4 5 6];               iw=3;jw=2;
  case 7,tblmagique=[1 2 3 4 5 6 7];             iw=4;jw=2;
  case 8,tblmagique=[1 2 3 4 5 6 7 8];           iw=4;jw=2;
  case 9,tblmagique=[1 2 3 4 5 6 7 8 9];         iw=3;jw=3;
 case 10,tblmagique=[1 2 3 4 5 6 7 8 9 10];      iw=4;jw=3;
 case 11,tblmagique=[1 2 3 4 5 6 7 8 9 10 11];   iw=4;jw=3;
 case 12,tblmagique=[1 2 3 4 5 6 7 8 9 10 11 12];iw=4;jw=3;
end



end %if



f1=figur;clf; hold on

for iN=1:NN
     subplot(iw,jw,tblmagique(iN)); hold on
C = squeeze(real(CHP(iN,:,:)));
%plotm(C./xtrm(C),typg);caxis([-abs(xtrm(C)) abs(xtrm(C))]);
plotm(C,typg);caxis([-abs(xtrm(C)) abs(xtrm(C))]);
if(iN==1);cx=caxis;end;caxis(cx);
colormap(mycolormap(mapanom,21));

%c=colorbar(CBAROR);
titre1=strcat('EOF',num2str(iN));
titre2=strcat(' (',num2str(expvar(iN)),'%)');
title(strcat(titre1,titre2));

end %for iN

%titre1=strcat(num2str(lat(il1)),'/',num2str(lat(il2)));
titre2=strcat(' Total explained variance:',num2str(sum(expvar)),'%');
%titre3=strcat('(Analysed file: ',ffile,')');
%titre=strvcat(titre2,titre3);
suptitle(titre2);

set(f1,'Position',[378 39 313 647]); % Laptop screen
%set(f1,'Position',[369 55 316 899]); % Desktop screen
%set(f1,'Name',strcat(ffile,'<>',titre1,'<> REAL PART'));
set(f1,'Name','<> REAL PART <>');


% -------------------------------------------
if isreal(CHP)==0

f2=figur;clf; hold on

for iN=1:NN
     subplot(iw,jw,tblmagique(iN)); hold on

C = squeeze(imag(CHP(iN,:,:)));
plotm(C,typg);caxis([-abs(xtrm(C)) abs(xtrm(C))]);
%if(iN==1);cx=caxis;end;
caxis(cx);
colormap(mycolormap(mapanom,21));

c=colorbar(CBAROR);
titre1=strcat('EOF',num2str(iN));
titre2=strcat(' (',num2str(expvar(iN)),'%)');
title(strcat(titre1,titre2));

end %for iNN

%titre1=strcat(num2str(lat(il1)),'/',num2str(lat(il2)));
titre2=strcat(' Total explained variance:',num2str(sum(expvar)),'%');
%titre3=strcat('(Analysed file: ',ffile,')');
%titre=strvcat(titre2,titre3);
suptitle(titre2);

set(f2,'Position',[378 39 313 647]); % Laptop screen
%set(f2,'Position',[369 55 316 899]); % Desktop screen
set(f2,'Name','<> IMAGINARY PART <>');

end


