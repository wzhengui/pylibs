function F=GetWODData(StaName, VarName, StartT, EndT, Sflag)
%Write WOD database for (Station, Variable, Time)
%
%format: Output=GetWODData(StaName, VarName, StartT, EndT, flag)
%  StaName: accept both cell and string of stations.
%           Examples:  'XBTO1207','OSDO',{'XBTO1207','PFLO1006'}, {'XBTO1207','OSDO'}
%  VarName: array of variable codes: [1,2 3]
%  StartT and EndT: Time range
%  flag=0: not save; flag=1: save database in .mat format not including CC
%       (country code); flag=1: save database in .mat format including CC
%  including all: ''
%Examples: 
%  Outputs=GetWODData({'XBTO1207','PFLO1006'});
%  Outputs=GetWODData('',[1]);
%  Outputs=GetWODData('XBTO1207',[1],datenum(2001,1,1),datenum(2002,1,1));
%  Outputs=GetWODData('XBTO1207',[1],'','',1);

filepath='D:\Work\Database\WOD\mat';
DataInfo=load(fullfile(filepath,'DataInfo.mat'));
%----pre proc--------------

%----StaName------------
if nargin==0
    StaName=DataInfo.filenames;
else 
    if isstr(StaName)
        if strcmp(StaName,'')|strcmp(StaName,' ')
            StaName=DataInfo.filenames;
        else
            StaName={StaName};
        end
    end
    
    flag=0;SN={};
    for r1=1:length(StaName)
        if length(StaName{r1})==4
            flag=flag+1;
            SN{flag,1}=StaName{r1};
            if flag==1
                fp=strcmp(DataInfo.filetypes,SN{flag});
            else
                fp=fp|strcmp(DataInfo.filetypes,SN{flag});
            end
        end
    end
    if length(SN)~=0
        StaName=[setdiff(StaName,SN);DataInfo.filenames(fp)];
    end    
end
StaName=unique(StaName);

%----VarName------------------
if nargin<=1
    VarName='';
end   

if strcmp(VarName,'')|strcmp(VarName,' ')
    VarName=DataInfo.vcode;
end

VarName=unique(VarName);
if isstr(VarName)
    VarName=str2num(VarName);
end

%----StartT and EndT-------------
if nargin<=2
    StartT=''; EndT='';
end
if nargin==3
    disp('too few argument');
    return;
end

if strcmp(StartT,'')|strcmp(StartT,' ')
    StartT=datenum(1000,1,1);
end
if strcmp(EndT,'')|strcmp(EndT,' ')
    EndT=datenum(3000,1,1);
end

%---flag-----------
if nargin<5
    Sflag=0;
end

if nargin>5
    disp('too many argument');
    return;
end

%-------extract results---------------
flag=0;
for r1=1:length(StaName)    
    %check whether VarName exists
    D=load(fullfile(filepath,StaName{r1}),'uvcode');
    cvar=intersect(D.uvcode,VarName);
    if isempty(cvar)
        continue;
    end
    
    D=load(fullfile(filepath,StaName{r1}));
    %---variable filter
    for r2=1:length(VarName)
        if r2==1
            fp=D.vcode==VarName(r2);
        else
            fp=fp|D.vcode==VarName(r2);
        end
    end
    
    fp=fp&D.Doy>=StartT&D.Doy<=EndT;
    
    if sum(fp)==0
        continue;
    end
    
    disp(['reading: ',StaName{r1}]);
    flag=flag+1;
    %Get Results-----
    if flag==1
        if Sflag>=2
            CC=D.CC_value(D.CC_code(fp));
        end
        Cast=D.Cast(fp);
        Data=D.Data(fp);
        Doy=D.Doy(fp);
        Latitude=D.Latitde(fp);
        Longitude=D.Longitde(fp);
        cruise=D.cruise(fp);
        nlevels=D.nlevels(fp);
        z=D.z(fp);
        vcode=D.vcode(fp);
    else
        if Sflag>=2
            CC=[CC;D.CC_value(D.CC_code(fp))];
        end
        Cast=[Cast;D.Cast(fp)];
        Data=[Data;D.Data(fp)];
        Doy=[Doy;D.Doy(fp)];
        Latitude=[Latitude;D.Latitde(fp)];
        Longitude=[Longitude;D.Longitde(fp)];
        cruise=[cruise;D.cruise(fp)];
        nlevels=[nlevels;D.nlevels(fp)];
        z=[z;D.z(fp)];
        vcode=[vcode;D.vcode(fp)];
    end    
end

if ~exist('vcode')
    F={};
    return;
end

uvcode=unique(vcode);

if Sflag>=2; F.CC=CC; end
F.Cast=Cast; F.Data=Data; F.Doy=Doy; F.Latitude=Latitude;
F.Longitude=Longitude; F.cruise=cruise; F.nlevels=nlevels; F.z=z; 
F.vcode=vcode; F.uvcode=uvcode;

if Sflag==1 
    save WODData Cast Data Doy Latitude Longitude cruise nlevels z vcode uvcode;
elseif Sflag==2
    save WODData CC Cast Data Doy Latitude Longitude cruise nlevels z vcode uvcode;
end

return;
end