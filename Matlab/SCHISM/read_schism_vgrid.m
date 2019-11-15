function [zcor,nvrt]=read_schism_vgrid(fname,np,node,dp,eta,h0,iflag)
%read SCHISM vertical grid information
%usage:
%   1) [zcor,nvrt]=read_schism_vgrid('vgrid.in',np,node,dp,eta,h0); %get zcor for all nodes
%   2) [zcor,nvrt]=read_schism_vgrid('vgrid.in',np,bndnode{1},dp(bndnode{1}),eta,h0); get zcor for a subset of nodes
%   3) [zcor,nvrt]=read_schism_vgrid('vgrid.in',np,node,dp,eta,h0,1); %get zcor for all nodes, replace NaN in zcor by the bottom depth


fid=fopen(fname);
Ti=textscan(fgetl(fid),'%d');
ivcor=Ti{1};

if ivcor==2 %no sanity check, assume a good vgrid.in
    %ivcor=2, nvrt,kz,h_s,ztot,nsig,h_c,theta_b,theta_f,sigma
    Ti=textscan(fgetl(fid),'%d%d%f');
    nvrt=Ti{1};kz=Ti{2};h_s=Ti{3};
    fgetl(fid);
    for r1=1:kz
        Ti=str2num(fgetl(fid));
        ztot(r1)=Ti(2);
    end
    nsig=nvrt-kz+1;
    fgetl(fid);
    Ti=textscan(fgetl(fid),'%f%f%f');
    h_c=Ti{1};theta_b=Ti{2};theta_f=Ti{3};
    for r1=1:nsig
        Ti=str2num(fgetl(fid));
        sigma(r1)=Ti(2);
    end
    zcor=Zcor_SZ(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma);
elseif ivcor==1
    Ti=textscan(fgetl(fid),'%d');
    nvrt=Ti{1};
    sigma=nan(nvrt,np);
    for r1=1:np
        Ti=str2num(fgetl(fid));
        sigma(Ti(2):nvrt,r1)=Ti(3:end);
    end
    zcor=sigma(:,node).*repmat(dp',nvrt,1);
else
    disp('wrong 2');
    return;
end
fclose(fid);

zcor=zcor';

if nargin==7&iflag==1
    zi0=zcor(:,nvrt);
    for r1=(nvrt-1):-1:1
        zi=zcor(:,r1);
        fp=isnan(zi);zi(fp)=zi0(fp);
        zcor(:,r1)=zi;
        zi0=zi;        
    end
end

end



function zcor=Zcor_SZ(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma)
%get z coodinate zcor, without sanity check
nsig=nvrt-kz+1;
cs=(1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f)+...
    theta_b*(tanh(theta_f*(sigma+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5);
hmod=min(dp,h_s);
zcor=nan(nvrt,length(dp));

%use fp to replace the loop, for sigma layer
fp1=hmod<=h_c;
zcor(kz:nvrt,fp1)=sigma'*(hmod(fp1)'+eta)+eta;
fp2=eta<=(-h_c-(hmod-h_c)*theta_f/sinh(theta_f));
if sum(fp2)>0
    disp('wrong');
    return;
else
    fp3=~fp1;
    zcor(kz:nvrt,fp3)=repmat(eta*(1+sigma')+h_c*sigma',1,sum(fp3))+cs'*(hmod(fp3)'-h_c);
end

% for z layer
kbp=zeros(length(dp),1);
kbp(dp<=h_s)=kz;
fp4=dp>h_s;ind=find(fp4);
for r1=ind
    if isempty(r1)
        break;
    end
    for k=1,kz-1
        if(-dp(r1)>=ztot(k)&-dp(r1)<=ztot(k+1))
            kbp(r1)=k;
            break;
        end
    end
    if kbp(r1)==0
        disp('can not find a bottom level for node');
        return;
    else kbp(r1)<1|kbp(r1)>=kz
        disp('imposisble kbp,kz');
        return;
    end
    zcor(kbp(r1),r1)=-dp;
    for k=kbp(r1)+1:kz-1
        zcor(k,r1)=ztot(k);
    end
end
end
