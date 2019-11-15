function F=GetCurrentData(StaName,VarName,StartT,EndT,flag)
%Write Current datebase for (Stations, Variabel, Time)
%
%format: GetCurrentData(StaName,VarName,StartT,EndT,flag)
%default variables included: Doy,Lat,Lon,CDir,CSpd,Depth,Orientation,
%available optional variables:  bin, nbin, platform_orientation,platform_pitch_angle, platform_roll_angle, SensorDepth, WTEMP
%
%  StaName: accept both cell and string of stations. eg. ,or 'cb0102'
%  VarName: same format for paramters. eg. ,or 'WTEMP', or bin
%  StartT and EndT: Time range
%  flag==1: save database in .mat format
%  eg.
%  Outputs=GetCurrentData()  %output station information
%  Outputs=GetCurrentData({'cb0201','cb0102'}) %outputs all parameters (2009-2014)
%  Outputs=GetCurrentData({'cb0201','cb0102'},{'WTEMP','bin'}) %default time range (2009-2014)
%  Outputs=GetCurrentData({'cb0201','cb0102'},' ') %default time range (2009-2014), default variables
%  Outputs=GetCurrentData({'cb0201','cb0102'},'All') %default time range (2009-2014), All variables
%  Outputs=GetCurrentData('cb0201','WTEMP',datenum(2010,1,1),datenum(2011,12,31),1)
%  Outputs=GetCurrentData(' ',' ',datenum(2009,1,1),datenum(2015,12,31),1) %all stations and default variables

DataBase='D:\Work\Database\Current\Data.mat\';
StaInfo=load([DataBase,'StaInfo.mat']);

DVar={'Doy','Lat','Lon','CDir','CSpd','Depth','Orientation'};
OVar={'bin','nbin','platform_orientation','platform_pitch_angle','platform_roll_angle','SensorDepth','WTEMP'};
%available optional variables:  bin, nbin, platform_orientation,platform_pitch_angle, platform_roll_angle, SensorDepth, WTEMP

if nargin==0
    format compact;
    disp('available stations:');
    disp(StaInfo.Station');
    disp('optional variables:');
    disp(OVar);
%     F.Lat=StaInfo.Lat;
%     F.Lon=StaInfo.Long;
%     F.Var=StaInfo.Var; 
%     F.Station=StaInfo.Station;
    F=load([DataBase,'StaInfo.mat']);
    F.VarName=F.Var; F.VarName0=F.Var0;
    F.Var=OVar;
    F=rmfield(F,'Var0');
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

%---for variables
if exist('VarName','var')&isstr(VarName)&strcmp(VarName,'All')
    VarName=OVar;
end

if exist('VarName','var')&isstr(VarName)
    VarName={VarName};
end

%-----
if nargin<2|strcmp(VarName,' ')
    %Variable Name
    VarName={};%output all parameter;
end

VarName=[DVar,setdiff(VarName,DVar)];


if nargin<4|(strcmp(StartT,' ')&strcmp(EndT,' '))
    %time range;    
    StartT=datenum('2009-01-01');EndT=datenum('2016-12-31');
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
F.Station={};
for r1=1:length(VarName)
    if strcmp(VarName{r1},'Orientation');
        eval(['F.',VarName{r1},'={};']);
    else
        eval(['F.',VarName{r1},'=[];']);
    end
end

for r1=1:length(StaName)
    SN=StaName{r1}; %Staind=find(strcmp(StaInfo.Station,SN));
    fname=[DataBase,SN,'.mat'];
    D=load(fname);
    
    %---filter date--
    fp=D.Doy>=StartT&D.Doy<=EndT;
    F.Station=[F.Station;repmat({StaName{r1}},sum(fp),1)];
    for r2=1:length(VarName)
        if strcmp(VarName{r2},'Orientation');
            F.Orientation=[F.Orientation;repmat({D.Orientation},sum(fp),1)];            
        else
            eval(['F.',VarName{r2},'=[F.',VarName{r2},';D.',VarName{r2},'(fp)];']);
        end
    end
end

if flag==1
    savestr='save CurrentData Station ';
    Station=F.Station;
    for r1=1:length(VarName)        
%         if strcmp(VarName{r1},'Orientation');
%             F.Orientation=[F.Orientation;repmat({D.Orientation},sum(fp),1)];
%             eval([VarName{r1},'=F.',VarName{r1},';']);
%         else
%             eval(['F.',VarName{r2},'=[F.',VarName{r2},';D.',VarName{r2},'(fp)];']);
%         end        
        savestr=[savestr,' ',VarName{r1}];
        eval([VarName{r1},'=F.',VarName{r1},';']);
    end
    savestr=[savestr,' -v7.3'];
    eval(savestr);
end

return;
end