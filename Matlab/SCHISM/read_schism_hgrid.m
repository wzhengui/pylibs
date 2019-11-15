%%%------------------------------------------------------------------------
%%read SCHISM hgrid.gr3
%%usage:
%%1)no boundary info: [ne,np,lx,ly,dp,i34,elnode]=read_schism_hgrid(fname)
%%2)with boundary info: [ne,np,lx,ly,dp,i34,elnode,bndinfo]=read_schism_hgrid(fname,1)
%%%------------------------------------------------------------------------
function [ne,np,lx,ly,dp,i34,elnode,bndinfo]=read_schism_hgrid(fname,ibnd)

%another way to read text; 
% [node(:,1),xx,yy,node(:,4)]=textread('hgrid.gr3', '%d %f %f %f',np, 'headerlines', 2);    

fid=fopen(fname);
fgetl(fid);
Ti=textscan(fgetl(fid),'%f%f');

ne=Ti{1}; np=Ti{2};
for r1=1:np
    Ti=str2num(fgetl(fid));
    lx(r1,1)=Ti(2);
    ly(r1,1)=Ti(3);
    dp(r1,1)=Ti(4);
end
elnode=nan(ne,4);
for r1=1:ne
     Ti=str2num(fgetl(fid));
     i34(r1,1)=Ti(2);
     elnode(r1,1:i34(r1))=Ti(3:end);
end

%------read boundary information
if nargin==2&ibnd==1
    
    %---open bnd info----
    Ti=textscan(fgetl(fid),'%f');bndinfo.nob=Ti{1};
    fgetl(fid);
    for r1=1:bndinfo.nob
        Ti=textscan(fgetl(fid),'%f');bndinfo.nobn(r1,1)=Ti{1};
        for r2=1:bndinfo.nobn(r1)
            Ti=textscan(fgetl(fid),'%f');bndinfo.iobn{r1}(r2,1)=Ti{1};
        end
    end
    
    %---land bnd info----    
    Ti=textscan(fgetl(fid),'%f');bndinfo.nlb=Ti{1};
    fgetl(fid);
    for r1=1:bndinfo.nlb
        linetext=fgetl(fid);
        Ti=textscan(linetext,'%f');bndinfo.nlbn(r1,1)=Ti{1}(1);
        if strfind(linetext,'island')
            bndinfo.island(r1,1)=0;
        else
            bndinfo.island(r1,1)=1;
        end
        for r2=1:bndinfo.nlbn(r1)
            Ti=textscan(fgetl(fid),'%f');bndinfo.ilbn{r1}(r2,1)=Ti{1};
        end
    end
    bndinfo.tag={'nob: number of open boundaries'; ...
                 'nobn: number of nodes for each open boundary'; ...
                 'iobn: nodes for each open boundary';...
                 'nlb: number of land boundaries'; ...
                 'nlbn: number of nodes for each land boundary'; ...
                 'ilbn: nodes for each land boundary';};
    
end
fclose(fid);

end