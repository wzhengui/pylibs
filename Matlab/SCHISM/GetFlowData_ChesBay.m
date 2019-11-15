function C=GetFlowData_ChesBay(StaName,StartT,EndT,flag)
%get Chespeake Bay River flow data
%format: C=GetFlowData(StaName,StartT,EndT,flag)
%eg.
%  1) GetFlowData(); %show all available stations
%  2) GetFlowData('ChoptankRiver') or GetFlowData({'ChoptankRiver','JameRiver_James'}); % Get all available Flow Data
%  3) GetFlowData(StaName,datenum('2001-01-01'),datenum('2015-01-01')) % Get Flow Data in certaion period
%  4) GetFlowData(StaName,datenum('2001-01-01'),datenum('2015-01-01'),1) save Flow Data in Matlab format
%     GetFlowData(StaName,' ',' ',1) save Flow Data in Matlab format

%----Database
DataBase='D:\Work\Database\River_flow_ChesBay\Data.mat\';
StaInfo=load([DataBase,'StaInfo.mat']);

%-show Station Info
if nargin==0
    disp(StaInfo.StaName(:));
    return;
end

if isstr(StaName)
    StaName={StaName};
end

Doy=[]; Flow=[]; StationNum=[]; StationName=[];
for r1=1:length(StaName)
    Data=load([DataBase,StaName{r1}]);
    if nargin==3&~strcmp(StartT,' ')&~strcmp(EndT,' ')
        fp=Data.Doy<StartT|Data.Doy>EndT;
        Data.Doy(fp)=[];
        Data.Flow(fp)=[];
    end  
    
    ind=strcmp(StaInfo.StaName,StaName{r1});
    StaNumi=StaInfo.StaNum(ind);
    StaNamei=StaInfo.StaName(ind);
    
    Doy=[Doy;Data.Doy];
    Flow=[Flow;Data.Flow];
    StationNum=[StationNum; repmat(StaNumi,length(Data.Doy),1)];
    StationName=[StationName; repmat(StaNamei,length(Data.Doy),1)];    
end

%--save data for future use
if nargin==4&flag==1
    StaNum=StationNum;
    StaName=StationName;
    save RiverFlowData StaNum StaName Doy Flow;
end

C.Doy=Doy;
C.Flow=Flow;
C.StaNum=StationNum;
C.StaName=StationName;
return;
end