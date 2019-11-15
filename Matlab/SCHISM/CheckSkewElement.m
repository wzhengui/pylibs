%this script is to show the skew elements in the grid. 
clear;clc;close all;

%---input--
fname='G2.gr3';
angle_critical=[15,145]; %min. and max.


%---read grid---
[ne,np,lx,ly,dp,i34,elnode]=read_schism_hgrid(fname);

%----create table---
nx=nan(3,4,4);
for r1=3:4
    for r2=1:r1
        for r3=1:(r1-1)
            tmp=r2+r3;
            if tmp>r1
                tmp=tmp-r1;
            end
            nx(r3,r2,r1)=tmp;
        end
    end
end


%---compute angle--------
angle=nan(ne,4);
for r1=1:ne
    for r2=1:i34(r1)
        i1=elnode(r1,r2);
        i2=elnode(r1,nx(1,r2,i34(r1)));
        i3=elnode(r1,nx(2,r2,i34(r1)));
        
        V1=[lx(i1)-lx(i2),ly(i1)-ly(i2)]; %p2->p1
        V2=[lx(i3)-lx(i2),ly(i3)-ly(i2)]; %p2->p1
        d1=sqrt(sum(V1.*V1));d2=sqrt(sum(V2.*V2));
        theta=acos(sum(V1.*V2)/(d1*d2))*180/pi; 
        angle(r1,r2)=theta;        
    end
end

angle_min=min(angle(:));
angle_max=max(angle(:));
disp(['Minmium angle is: ',num2str(angle_min)]);
disp(['Maximum angle is: ',num2str(angle_max)]);

%----find angle beyond critical--
Emin=[]; Emax=[];
for r1=1:4
    tmin=find(angle(:,r1)<=angle_critical(1));
    tmax=find(angle(:,r1)>=angle_critical(2));    
    Emin=[Emin; tmin];
    Emax=[Emax; tmax];
end
Emin=unique(Emin);
Emax=unique(Emax);
EV=unique([Emin;Emax]);

%---plot grid------
ind1=i34==3; ind2=i34==4;
patch('Faces',elnode(ind1,1:3),'Vertices',[lx,ly],'FaceVertexCdata',dp,...
  'FaceColor','none','EdgeColor','k'); hold on;
patch('Faces',elnode(ind2,1:4),'Vertices',[lx,ly],'FaceVertexCdata',dp,...
  'FaceColor','none','EdgeColor','k'); hold on;

%--EV
flag=1; ds=2e4;
for r1=1:length(EV)
    ei=EV(r1);
    xi=lx(elnode(ei,1:i34(ei)));
    yi=ly(elnode(ei,1:i34(ei)));
    xi=[xi;xi(1)];
    yi=[yi;yi(1)];
    patch(xi,yi,'r'); hold on;  
    if flag==1
        set(gca,'xlim',[mean(xi)-ds,mean(xi)+ds],'ylim',[mean(yi)-ds,mean(yi)+ds]);
        pause(1);
    end
end



