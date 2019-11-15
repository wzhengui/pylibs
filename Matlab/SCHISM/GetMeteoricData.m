function F=GetMeteoricData(StaName,VarName,StartT,EndT,flag)
%Write WQ datebase for (Stations, Variabel, Time)
%
%format: GetMeteoricData(StaName,VarName,StartT,EndT,flag)
%  StaName: accept both cell and string of stations. eg. ,or 'bltm2'
%  VarName: same format for paramters. eg. ,or 'WSPD'
%  StartT and EndT: Time range
%  flag==1: save database in .mat format
%  eg.
%  Outputs=GetMeteoricData()  %output station information
%  Outputs=GetMeteoricData({'bltm2','44063'}) %outputs all parameters (1990-2015)
%  Outputs=GetMeteoricData({'bltm2','44063'},{'WDIR','WSPD','ATMP'}) %default time range (1990-2015)
%  Outputs=GetMeteoricData('bltm2','WDIR',datenum(2001,1,1),datenum(2001,12,31),1)
%  Outputs=GetMeteoricData(' ',' ',datenum(2001,1,1),datenum(2001,12,31),1) %all stations and all variables

DataBase='D:\Work\Database\Meteorology\Data.mat\';

StaInfo=load([DataBase,'StaInfo.mat']);
if nargin==0
    format compact;
    disp(StaInfo.Station(StaInfo.DataAvailability)');
    
    F.Station=StaInfo.Station(StaInfo.DataAvailability);
    F.Lat=StaInfo.Lat(StaInfo.DataAvailability);
    F.Lon=StaInfo.Lon(StaInfo.DataAvailability);
    F.StationNote=StaInfo.StationNote(StaInfo.DataAvailability);
    F.Var=StaInfo.Var;  
    return;
end

%---------pre-proc.----
if isstr(StaName)
    if strcmp(StaName,' ')
        StaName=StaInfo.Station;
    else
        StaName={StaName};
    end
end
if exist('VarName','var')&isstr(VarName)
    VarName={VarName};
end

%-----
if nargin<2|strcmp(VarName,' ')
    %Variable Name
    VarName=StaInfo.Var;%output all parameter;
end

if nargin<4|(strcmp(StartT,' ')&strcmp(EndT,' '))
    %time range;    
    StartT=datenum('1990-01-01');EndT=datenum('2015-12-31');
end

if nargin<5
    flag=0;
end

%---delete station without data----
StaName=setdiff(StaName,StaInfo.Station(~StaInfo.DataAvailability));

%---get data---
if length(StaName)==0
    F=[];
else
    F=OutputDataBase(StaName,VarName,StartT,EndT,flag,DataBase);
end

end

function F=OutputDataBase(StaName,VarName,StartT,EndT,flag,DataBase)
%----output databased based StaName,VarName,StarT,EndT,flag
StaInfo=load([DataBase,'StaInfo.mat']);

%-------Initial Variables---
F.Doy=[];F.Station={};F.Lat=[];F.Lon=[];
for r1=1:length(VarName)
    eval(['F.',VarName{r1},'=[];']);
end

for r1=1:length(StaName)
    SN=StaName{r1}; Staind=find(strcmp(StaInfo.Station,SN));
    fname=[DataBase,SN,'.mat'];
    D=load(fname);
    
    %---filter date--
    fp=D.Doy>=StartT&D.Doy<=EndT;    
    F.Doy=[F.Doy;D.Doy(fp)]; 
    F.Station=[F.Station;repmat({StaName{r1}},sum(fp),1)];
    F.Lat=[F.Lat;ones(sum(fp),1)*StaInfo.Lat(Staind)];
    F.Lon=[F.Lon;ones(sum(fp),1)*StaInfo.Lon(Staind)];
    for r2=1:length(VarName)
        eval(['F.',VarName{r2},'=[F.',VarName{r2},';D.',VarName{r2},'(fp)];']);
    end    
end

if flag==1
    savestr='save MeteoricData Doy Station Lat Lon';
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