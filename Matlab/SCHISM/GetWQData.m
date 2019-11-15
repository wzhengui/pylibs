function F=GetWQData(StaName,VarName,StartT,EndT,flag)
%Write WQ datebase for (Stations, Variabel, Time)
%
%format: GetWQData(StaName,VarName,StartT,EndT,flag)
%  StaName: accept both cell and string of stations. eg. ,or 'CB3.2'
%  VarName: same format for paramters. eg. ,or 'CHLA'
%  StartT and EndT: Time range
%  flag==1: save database in .mat format
%  eg.
%  Outputs=GetWQData({'CB3.2','CB3.3C'}) %outputs all parameters (1984-2015)
%  Outputs=GetWQData({'CB3.2','CB3.3C'},{'WTEMP','SALINITY','CHLA','DO'}) %default time range (1984-2015)
%  Outputs=GetWQData('CB3.2','DO',datenum(2001,1,1),datenum(2001,12,31),1)

%---------pre-proc.----
if isstr(StaName)
    StaName={StaName};
end
if exist('VarName','var')&isstr(VarName)
    VarName={VarName};
end

if nargin==0
    %-----------------input for database output--------------------------------
    %Station Name
    StaName={'CB3.2','CB3.3C'};
end

if nargin<2|strcmp(VarName,' ')
    %Variable Name
    % VarName={'DO','SALINITY','WTEMP','TSS'};
    VarName={};%output all parameter;
end

if nargin<4|(strcmp(StartT,' ')&strcmp(EndT,' '))
    %time range;
    % StartT=datenum('2000-01-01');EndT=datenum('2010-01-01');
    StartT=datenum('1984-01-01');EndT=datenum('2015-12-31');
end

if nargin<5
    flag=0;
end
%get main stations for ChesBay---------------------------------------------
% filepath='E:\Work\WQ_Database\Data.mat\';
% load([filepath,'StaInfo'],'StationNum','StationName');
% flag=0;
% for r1=1:length(StationNum)
%     SN=StationName{r1};
%     if ~(sum(strfind(SN,'CB'))|sum(strfind(SN,'ET'))|sum(strfind(SN,'WT')))            
%         continue;
%     end
%     if sum(strfind(SN,'RET'))
%         continue;    
%     end
% %     disp(SN);
%     flag=flag+1;
%     StaName{flag}=SN;
% end
% --------------------------------------------------------------------------


F=OutputDataBase(StaName,VarName,StartT,EndT,flag);

end



function F=OutputDataBase(StaName,VarName,StartT,EndT,Wflag)
%-------------output database based on (StaName,VarName,StartT,EndT)-------
filepath='D:\Work\Database\WQ\Data.mat\'; %Umaine_Laptop
% filepath='E:\Work\Database\WQ\Data.mat\'; %UMaine_Desktop
load([filepath,'StaInfo'],'StationNum','StationName');
format compact;
disp('output database for these stations: ');
disp(StaName);
[tmp,ind1,ind2]=intersect(StaName,StationName);
StaNum=StationNum(ind2);

flag=0;
for r1=1:length(StaNum)
    SN=StaName{r1};
    fname=[filepath,'WQData_',num2str(StaNum(r1))];
    
    flag=flag+1;
    D=load(fname);
    
    %-----filter Variables
    if length(VarName)~=0
        for r2=1:length(VarName)
            if r2==1
                fp=strcmp(D.Var,VarName{r2});
            else
                fp=fp|strcmp(D.Var,VarName{r2});
            end
        end
        fp=fp&D.Data~=-9999;
    else
        fp=D.Data~=-9999;
    end
    
    if StartT~=0
        fp=fp&D.Doy>=StartT&D.Doy<=EndT;
    end
    
    
    if flag==1
        Data=D.Data(fp);
        Var=D.Var(fp);
        Station=D.Station(fp);
        Depth=D.Depth(fp);
        Layer=D.Layer(fp);
        Doy=D.Doy(fp);
        Lat=D.Lat(fp);
        Long=D.Long(fp);
        Unit=D.Unit(fp);
        TDepth=D.TDepth(fp);
        LPyc=D.LPyc(fp);
        UPyc=D.UPyc(fp);
    else
        Data=[Data; D.Data(fp)];
        Var=[Var; D.Var(fp)];
        Station=[Station; D.Station(fp)];
        Depth=[Depth; D.Depth(fp)];
        Layer=[Layer; D.Layer(fp)];
        Doy=[Doy; D.Doy(fp)];
        Lat=[Lat; D.Lat(fp)];
        Long=[Long; D.Long(fp)];
        Unit=[Unit; D.Unit(fp)];
        TDepth=[TDepth; D.TDepth(fp)];
        LPyc=[LPyc; D.LPyc(fp)];
        UPyc=[UPyc; D.UPyc(fp)];
    end
end

F.Doy=Doy; F.Data=Data; F.Var=Var; F.Station=Station; F.Depth=Depth; 
F.Layer=Layer; F.Lat=Lat; F.Long=Long; F.TDetph=TDepth; F.Unit=Unit;

% save CBWQdata Doy Data Var Station Depth Layer Lat Long TDepth;
% save WQdata Doy Data Var Station Depth Layer Lat Long TDepth Unit LPyc UPyc;
if Wflag==1
    save WQdata Doy Data Var Station Depth Layer Lat Long TDepth Unit;
end
end


















% function GetStationName()
% % %-----------Get StationName------------------
% clear;clc;close all;
% load StaInfo;
% 
% flag=0;
% for r1=1:length(StaNum)
%     fname=['./Data.mat/WQData_',num2str(StaNum(r1))];
%     if ~exist([fname,'.mat'])
%         disp('worng 1');
%         continue;
%     end
%     eval(['load ',fname,' Station;']);
%     if isempty(Station)
%         continue;
%     end
% 
%     flag=flag+1;
%     StationNum(flag)=StaNum(r1);
%     Usta=unique(Station);
%     if length(Usta)~=1
%         disp('wrong');
%         return;
%     end
%     StationName{flag}=Usta{:};
%     clear Station;
% end
% save StaInfo ParNum StaNum StationNum StationName;
% %-------------------------------------------
% end