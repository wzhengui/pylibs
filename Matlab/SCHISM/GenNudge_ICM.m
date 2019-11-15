%function GenNudge_ICM()
%generate ICM nudge files
%the code is similar to 3D boundary condition(salt3D.th, temp3D.th,and elev2D.th);

%--------------input parameters---------------
StartT=datenum('2012-01-01');
EndT=datenum('2014-12-31');

%---start time in *.txt
DataStartT=datenum('2012-01-01');
Out_h=[0:-1:-35];  %depth in *.txt

%--add library
add_matlabfun;

%--read grid info
[ne,np,lx,ly,dp,i34,elnode]=read_schism_hgrid('hgrid.gr3');
[zcor,nvrt]=read_schism_vgrid('vgrid.in',np,1:np,dp,0,0.1,1);

%--read ICM data
ntr=25;
for r1=1:ntr
    fname=['./CB5.1_',num2str(r1),'.txt'];
    Ti=load(fname);
    Doy=Ti(:,1)+DataStartT;
    Data(r1,:,:)=Ti(:,2:end);
end

%-------------------------------
fp=Doy>=StartT&Doy<=EndT;
ti=Doy(fp);yi=Data(:,fp,:); %(tr,time,nvrt)

tic;
%--calulcate and store the interpolation index
nnode=np;
indR=nan(nvrt,nnode,3);
for r1=1:nnode
    disp(['finished:', num2str(r1/nnode*100),'%']); 
    for r2=1:nvrt
        dpi=zcor(r1,r2);
        if isnan(dpi)
           disp('wrong,dpi=nan');
        end

        r=nan;
        for r3=1:length(Out_h)-1
            if dpi>=Out_h(1)                
                ind1=1;
                ind2=2;
                r=1;
                break;
            elseif dpi<=Out_h(end)
                ind1=length(Out_h)-1;
                ind2=length(Out_h);
                r=0;
                break;
            elseif dpi<=Out_h(r3)&dpi>=Out_h(r3+1)
                ind1=r3;ind2=r3;
                r=(dpi-Out_h(r3+1))/(Out_h(r3)-Out_h(r3+1));   
            end                       
        end 
        if isnan(r)|r>1|r<0
            disp('wrong');
            return;
        end
        indR(r2,r1,:)=[ind1,ind2,r];
    end    
end
toc;

tic;
%---do interpolation
ind1=indR(:,:,1); ind2=indR(:,:,2); r=indR(:,:,3);
ind1=ind1(:); ind2=ind2(:); r=r(:);

fid=fopen('ICM_nu.in.1','w+');
fwrite(fid,[ntr,nvrt,nnode,length(ti)],'int32');

tr=nan(ntr,nvrt,nnode);
for r1=1:length(ti)
  disp(datestr(ti(r1)));
  for r2=1:ntr
    yii=squeeze(yi(r2,r1,ind1)).*r+squeeze(yi(r2,r1,ind2)).*(1-r);
    tr(r2,:,:)=reshape(yii,nvrt,nnode); 
  end
  fwrite(fid,[86400*(ti(r1)-StartT);tr(:)],'float32');
end
fclose(fid);
toc;

%end

