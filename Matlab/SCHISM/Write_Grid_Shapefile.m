function Write_Grid_Shapefile(gridname,fname,prj)
%-------read grid info and write out shapefile for ArcGis------
%written by Zhengui Wang on 5/13/2016
%availabel proj: utm10N,utm18N,utm19N,ll
%eg: Write_Grid_Shapefile('hgrid.gr3','York','utm18N'); %utm projection
%    Write_Grid_Shapefile('hgrid.gr3','York','ll'); %ll projection
%------------------------------------------------------------

%------read grid info----
% fname='hgrid';
% fid=fopen('./hgrid.gr3');
% fname='York';
% fid=fopen('./York.gr3');

if nargin==3
    fid=fopen(gridname);
    %---write projection file----------
    if strcmp(prj,'utm10N')
        proj_str='PROJCS["NAD_1983_UTM_Zone_10N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-123.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0],AUTHORITY["EPSG",26910]]';
    elseif strcmp(prj,'utm18N')
        proj_str='PROJCS["NAD_1983_UTM_Zone_18N",GEOGCS["GCS_North_American_1983",DATUM["D_NORTH_AMERICAN_1983",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["Meter",1]]';
    elseif strcmp(prj,'utm19N')       
        proj_str='PROJCS["NAD_1983_UTM_Zone_19N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-69.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0],AUTHORITY["EPSG",26919]]';
    elseif strcmp(prj,'ll')
        proj_str='GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]';
    else
        disp('wrong prj');
        return;
    end
else
    disp('need 3 arguments');
    return;
end

fgetl(fid);
Ti=textscan(fgetl(fid),'%f%f'); ne=Ti{1}; np=Ti{2};
lx=nan(np,1); ly=nan(np,1); dp=nan(np,1);
for r1=1:np
    Ti=str2num(fgetl(fid));
    lx(r1,1)=Ti(2); ly(r1,1)=Ti(3); dp(r1,1)=Ti(4); 
end

i34=nan(ne,1);elnode=nan(ne,4);
for r1=1:ne
    Ti=str2num(fgetl(fid));
    i34(r1,1)=Ti(2);
    elnode(r1,1:i34(r1,1))=Ti(3:end);
end
fclose(fid);

%---write shapefile for elements--------
for r1=1:ne
    if mod(r1,1e4)==0
        disp(['element #: ',num2str(r1)]);
    end
    lxi=lx(elnode(r1,1:i34(r1)));
    lyi=ly(elnode(r1,1:i34(r1)));    
    BoundingBox=[min(lxi),min(lyi); max(lxi),max(lyi)];
    X=flipud(lxi); X=[X;X(1);nan];
    Y=flipud(lyi); Y=[Y;Y(1);nan];           
    E(r1,1)=struct('Geometry','Polygon','BoundingBox',BoundingBox,'X',X','Y',Y','ElementNumber',r1);
    if r1==1
        E=repmat(E,ne,1);
    end
end

%---write shapefile for nodes-----------
for r1=1:np
    if mod(r1,1e4)==0
        disp(['node #: ',num2str(r1)]);
    end
    N(r1,1)=struct('Geometry','Point','X',lx(r1),'Y',ly(r1),'NodeNumber',r1);
    if r1==1
        N=repmat(N,np,1);
    end
end


%--------------------------------------------------------------------------
disp(['writing node shapefile']);
fid=fopen([fname,'_node.prj'],'w+'); fprintf(fid,'%s',proj_str); fclose(fid);
disp(['writing element shapefile']);
fid=fopen([fname,'_element.prj'],'w+'); fprintf(fid,'%s',proj_str); fclose(fid);

shapewrite(N,[fname,'_node.shp']);
shapewrite(E,[fname,'_element.shp']);
end