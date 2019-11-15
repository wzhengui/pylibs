function [xtick,xticklabel]=plot_xtick(xi,format,iflag)
% %return xtick and xticklabel
% usage: [xtick,xticklabel]=plot_xtick(datenum(2001,1,1):datenum(2017,1,1),'yyyy');
% if flag is user defined format
% %---code for generate xtick----------------------------------
% flag=0;
% for r1=2012:2014
%     flag=flag+1;
%     xtick(flag)=datenum(r1,1,1);
%     xticklabel{flag}=datestr(xtick(flag),'yyyy');
% end
% set(gca,'xtick',xtick,'xticklabel',xticklabel);
% set(gca,'xlim',[datenum(2012,1,1),datenum(2015,1,1)]);

%-----general format---
if nargin==2
    for r1=1:length(xi)
        xticklabel{r1}=datestr(xi(r1),format);
    end
    xtick=xi;
end

%---user self-defined format
if nargin==3&iflag==1
    flag=0;
    for r1=xi
        flag=flag+1;
        xtick(flag)=datenum(r1,1,1);
        xticklabel{flag}=datestr(xtick(flag),'yy');
    end
end


%---user self-defined format
if nargin==3&iflag==2
    flag=0;
    for r1=2002:2016
        for r2=1:1:12
            flag=flag+1;
            if r2==1
                xtick(flag)=datenum(r1,r2,1);
%                 xticklabel{flag}=datestr(xtick(flag),'yyyy');
                xticklabel{flag}=datestr(xtick(flag),'mmm');
            else
                xtick(flag)=datenum(r1,r2,1);
                xticklabel{flag}=datestr(xtick(flag),'mmm');
            end
        end
    end    
end

%---user self-defined format
if nargin==3&iflag==3
    flag=0;
    for r1=1:13
        flag=flag+1;
        xtick(flag)=datenum(2005,r1,1);
        xticklabel{flag}=datestr(xtick(flag),'mmm');
        xtick(flag)=xtick(flag)-datenum(2005,1,1);
    end
end

if nargin==3&iflag==4
    flag=0;
    for r1=1900:2100
        flag=flag+1;
        xtick(flag)=datenum(r1,1,1);
        xticklabel{flag}=datestr(xtick(flag),format);
    end
%     set(gca,'xtick',xtick,'xticklabel',xticklabel);
%     set(gca,'xlim',[datenum(2012,1,1),datenum(2015,1,1)]);
end

if nargin==3&iflag==5
    flag=0;
    for r1=1:13
        flag=flag+1;
        
        xtick(flag)=datenum(2012,r1,1);
        if r1==1|r1==13
            xticklabel{flag}=datestr(xtick(flag),'yyyy');
        else
            xticklabel{flag}=datestr(xtick(flag),'mmm');
        end
        
    end
end

return;
end