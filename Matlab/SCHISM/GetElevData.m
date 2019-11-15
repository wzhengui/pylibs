function F=GetElevData(StaName,VarName,StartT,EndT,flag)
%Write Elev datebase for (Stations, Variabel, Time)
%
%format: GetElevData(StaName,VarName,StartT,EndT,flag)
%  StaName: accept station number, cell of station name, or string of stations.
%         eg. [8573927,8574680], or 'Chesapeake City', or {'Chesapeake, City','Baltimore'};
%  VarName: same format for paramters. eg. ,or 'Elev'
%  StartT and EndT: Time range
%  flag==1: save database in .mat format
%  eg.
%  Outputs=GetElevData()  %output station information
%  Outputs=GetElevData([8573927,8574680]) %outputs all parameters (1980-2015)
%  Outputs=GetElevData([8573927,8574680],{'Elev'}) %default time range (1980-2015)
%  Outputs=GetElevData(8573927,'Elev',datenum(2003,1,1),datenum(2015,12,31),1)
%  Outputs=GetElevData(' ',' ',datenum(2001,1,1),datenum(2015,12,31),1) %all stations and all variables

DataBase='D:\Work\Database\Elevation\Data.mat\';
StaInfo=load([DataBase,'StaInfo.mat']);

if nargin==0
    format compact;
    fp=StaInfo.DataAvailability==1;
    disp(StaInfo.Station(fp));
    F.Lat=StaInfo.Lat(fp);
    F.Lon=StaInfo.Lon(fp);
    F.Var=StaInfo.Var(fp); 
    F.Station=StaInfo.Station(fp);
    F.StaName=StaInfo.StaName(fp);
    return;
end

%---------pre-proc.----
if isstr(StaName)
    if strcmp(StaName,' ')
        fp=StaInfo.DataAvailability==1;
        StaName=StaInfo.StaName(fp);
    else
        StaName={StaName};
    end
end

if ~isnumeric(StaName)
    StaNamei=StaName; StaName=[];
    for r1=1:length(StaNamei)
        StaName(r1)=StaInfo.Station(find(strcmp(StaInfo.StaName,StaNamei{r1})));
    end
end


if exist('VarName','var')&isstr(VarName)
    VarName={VarName};
end

%-----
if nargin<2|strcmp(VarName,' ')
    %Variable Name
%     VarName=StaInfo.Var;%output all parameter;
    VarName={'Elev'};
end

if nargin<4|(strcmp(StartT,' ')&strcmp(EndT,' '))
    %time range;    
    StartT=datenum('1980-01-01');EndT=datenum('2016-12-31');
end

if nargin<5
    flag=0;
end

%---get data---
F=OutputDataBase(StaName,VarName,StartT,EndT,flag,DataBase);

end


function F=OutputDataBase(StaName,VarName,StartT,EndT,flag,DataBase)
%----output databased based StaName,VarName,StarT,EndT,flag
StaInfo=load([DataBase,'StaInfo.mat']);

%-------Initial Variables---
F.Doy=[];F.Station=[];F.Lat=[];F.Lon=[];F.StaName={};
for r1=1:length(VarName)
    eval(['F.',VarName{r1},'=[];']);
end

for r1=1:length(StaName)
    Staind=find(StaInfo.Station==StaName(r1));    
    SN=num2str(StaName(r1)); 
    fname=[DataBase,SN,'.mat'];
    D=load(fname);
    
    %---filter date--
    fp=D.Doy>=StartT&D.Doy<=EndT;    
    F.Doy=[F.Doy;D.Doy(fp)]; 
    F.Station=[F.Station;repmat(StaName(r1),sum(fp),1)];
    F.StaName=[F.StaName;repmat(StaInfo.StaName(Staind),sum(fp),1)];
    F.Lat=[F.Lat;D.Lat(fp)];
    F.Lon=[F.Lon;D.Lon(fp)];
    for r2=1:length(VarName)
        eval(['F.',VarName{r2},'=[F.',VarName{r2},';D.',VarName{r2},'(fp)];']);
    end    
end

if flag==1
    savestr='save ElevData Doy Station StaName Lat Lon';
    Doy=F.Doy; Station=F.Station; Lat=F.Lat; Lon=F.Lon;
    for r1=1:length(VarName)
        savestr=[savestr,' ',VarName{r1}];
        eval([VarName{r1},'=F.',VarName{r1},';']);
    end
    savestr=[savestr,' -v7.3'];
    eval(savestr);
end

return;
end