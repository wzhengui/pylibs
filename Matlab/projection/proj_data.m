function [xout,yout]=proj_data(xin,yin,proj,method)
%usage
%  S=proj_xy(): output utm_zone list;
%  [xout,yout]=proj_xy(xin,yin,proj,method)
%   xin: x coordinates or Lon
%   yin: y coordinates or Lat
%   proj: projection stcture data or projection names
%        1st eg: proj=defaultm('utm');
%        2nd eg: proj='epsg:26918';
%   method: (1: lon&lat to xy; 2: xy to lon&lat;


filepath='D:\OneDrive\Matlab\projection\';
if nargin==0
    fname=[filepath,'utm_zones.txt'];
    fid=fopen(fname);
    flag=0;
    while(~feof(fid))
        flag=flag+1;
        xout{flag,1}=fgetl(fid);
    end
    fclose(fid);
    return;
end

proj0=proj;
if isstr(proj0)
    proj=GetProjection(proj0,filepath);
end

if method==1
    [xout,yout]=mfwdtran(proj,yin,xin); %from LL to XY    
else
    [yout,xout]=minvtran(proj,xin,yin);    %from XY to LL
end

end


function proj=GetProjection(proj_str,filepath0)
%-----prepare proj structure data-----
filepath=[filepath0,'prj_files\'];
Ti=strsplit(proj_str,':');
fname=[filepath,Ti{1},'.',Ti{2},'.prj'];
[projdata,VarName,Var]=parse_prj_file(fname);
% [projdata,VarName,Var]=parse_prj_file('D:Work/ChesBay/Loading/P6/WatershedSegment/nps.prj');

%---get projection id----------
projname=projdata.PROJECTION.Name;
if strcmp(projname,'Transverse_Mercator')
    projname='Universal Transverse Mercator (UTM)';
elseif strcmp(projname,'Lambert_Conformal_Conic')
    projname='Lambert Conformal Conic';
elseif strcmp(projname,'Albers')
    projname='Equal Area Conic (Albers)';
end

%---save maplist---------------------------------------------
mlistdata=maplist;
for r1=1:length(mlistdata)
    mlist.ID{r1,1}=mlistdata(r1).IdString;
    mlist.Name{r1,1}=mlistdata(r1).Name;
end

ind=find(strcmp(mlist.Name,projname));
proj_id=mlist.ID{ind};
proj=defaultm(proj_id);

%--get zone number---
PN=replace(replace(projdata.Name,'_',' '),'Zone','zone');
Ti=strsplit(PN,' '); ind=find(strcmp(Ti,'zone'));
if length(ind)==1
    zone=Ti{ind+1};
    proj.zone=zone;
end

%----------Modify Var List-------------
SVar={'GEOGCS.DATUM.SPHEROID.Value','false_easting','false_northing','scale_factor',...
    'latitude_of_origin','central_meridian','standard_parallel_1','standard_parallel_1',...
    };
DVar={'geoid','falseeasting','falsenorthing','scalefactor',...
    'origin(1)','origin(2)','mapparallels(1)','mapparallels(2)',...
    };

for r1=1:length(SVar)   
    ind=strfind(SVar{r1},'.');
    if length(ind)>0
        vn=['projdata.',SVar{r1}(1:ind(end)-1)];
        SVari=SVar{r1}(ind(end)+1:end);
    else
        vn='projdata';
        SVari=SVar{r1};
    end
    eval(['vi=',vn,';']);
    
      
    if isfield(vi,SVari)
        eval(['vii=vi.',SVari,';']);
        if strcmp(SVar{r1},'GEOGCS.DATUM.SPHEROID.Value')
            viii=1/vii(2);
            vii(2)=sqrt(2*viii-viii*viii);
        end
        eval(['proj.',DVar{r1},'=vii;']);
        if strcmp(SVari,'standard_parallel_1')
            proj.origin(3)=0;
        end
    end
end
proj=defaultm(proj);
%----------------------------------------------


end

