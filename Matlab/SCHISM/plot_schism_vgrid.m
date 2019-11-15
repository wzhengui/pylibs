function hf=plot_schism_vgrid(xi,yi,zcor,flag)
%plot vertical profile given a schism vertical profile
%usage:
%1) dispaly vertial grid profile in 3D space,x and y are the coordinates of nodes
%   plot_schism_vgrid(x,y,zcor,1);
%2) display vertial grid profile in 2D space, 
%    plot_schism_vgrid(x,y,zcor,2);      %diff(xi)=abs(pt(i+1)-pt(i));
%    plot_schism_vgrid(x,'',zcor,2);     %xi=dist;
%    plot_schism_vgrid('','',zcor,2);    %xi=1:size(zcor,1);

figure;

np=size(zcor,1); nvrt=size(zcor,2);

%---plot vertial grid profile in 3D space
if(flag==1)
    for r1=1:nvrt
        plot3(xi,yi,zcor(:,r1),'.-'); hold on;
    end
    for r1=1:np
        xii=repmat(xi(r1),1,nvrt);
        yii=repmat(yi(r1),1,nvrt);
        plot3(xii,yii,zcor(r1,:),'.-'); hold on;
    end
    hold off;
    
end

%---plot vertial grid profile in 2D space
if(flag==2)
    if(~isempty(xi)&&~isempty(yi))
        dist(1)=0;
        for r1=1:length(xi)-1
            dist(r1+1)=dist(r1)+sqrt((xi(r1+1)-xi(r1))^2+(yi(r1+1)-yi(r1))^2);
        end                
    elseif(~isempty(xi)&&isempty(yi))
        dist=xi;
    elseif(isempty(xi)&&isempty(yi))
        dist=0:np-1;
    end
    
    %----plot vertical profile
    for r1=1:size(zcor,2)
        plot(dist,zcor(:,r1),'.-');hold on;
    end
    
    for r1=1:size(zcor,1)
        yi=zcor(r1,:);ti=ones(size(yi))*dist(r1);
        plot(ti,yi,'.-');hold on;
    end
    hold off;
end

set(gcf,'position',[566   478   963   490]);
hf=gcf;
end