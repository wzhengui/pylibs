function C=ReadSflux(fname,flag)
%this script is to read narr files
%C=ReadSflux(fname,flag);
%flag==1: sflux_air
%flag==2: sflux_prc
%flag==3: sflux_rad


if flag==1
    ncid=netcdf.open(fname,'NC_NOWRITE'); %open the file
    
    VarName={'time','lon','lat','uwind','vwind','prmsl','stmp','spfh'};
    AttName={{'long_name','standard_name','units','base_date'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'}
        };
    
    for r1=1:3 %read dimension info
        [C.dimname{r1},C.dims(r1)]=netcdf.inqDim(ncid,r1-1);
    end
    
    for r1=1:length(VarName) %read variables
        vid=netcdf.inqVarID(ncid,VarName{r1});
        
        %---------get attributes of variable
        for r2=1:length(AttName{r1})
            attval=netcdf.getAtt(ncid,vid,AttName{r1}{r2}); %read attributes
            eval(['C.',VarName{r1},'.',AttName{r1}{r2},'=attval;']);
        end
        %--------------get values of variables
        val=netcdf.getVar(ncid,vid); % read variable values
        if r1>=1
            val=double(val);
        end
        eval(['C.',VarName{r1},'.val=val;']);
    end
    netcdf.close(ncid);
    
    %----modify time
    C.time.val=C.time.val+datenum([double(C.time.base_date),0,0]);
    
elseif flag==2
    ncid=netcdf.open(fname,'NC_NOWRITE'); %open the file
    
    VarName={'time','lon','lat','prate'};
    AttName={{'long_name','standard_name','units','base_date'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'}
        };
    
    for r1=1:3 %read dimension info
        [C.dimname{r1},C.dims(r1)]=netcdf.inqDim(ncid,r1-1);
    end
    
    for r1=1:length(VarName) %read variables
        vid=netcdf.inqVarID(ncid,VarName{r1});
        
        %---------get attributes of variable
        for r2=1:length(AttName{r1})
            attval=netcdf.getAtt(ncid,vid,AttName{r1}{r2}); %read attributes
            eval(['C.',VarName{r1},'.',AttName{r1}{r2},'=attval;']);
        end
        %--------------get values of variables
        val=netcdf.getVar(ncid,vid); % read variable values
        if r1>=1
            val=double(val);
        end
        eval(['C.',VarName{r1},'.val=val;']);
    end
    netcdf.close(ncid);
    %----modify time
    C.time.val=C.time.val+datenum([double(C.time.base_date),0,0]);
elseif flag==3
    ncid=netcdf.open(fname,'NC_NOWRITE'); %open the file
    
    VarName={'time','lon','lat','dlwrf','dswrf'};
    AttName={{'long_name','standard_name','units','base_date'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'},...
        {'long_name','standard_name','units'}
        };
    
    for r1=1:3 %read dimension info
        [C.dimname{r1},C.dims(r1)]=netcdf.inqDim(ncid,r1-1);
    end
    
    for r1=1:length(VarName) %read variables
        vid=netcdf.inqVarID(ncid,VarName{r1});
        
        %---------get attributes of variable
        for r2=1:length(AttName{r1})
            attval=netcdf.getAtt(ncid,vid,AttName{r1}{r2}); %read attributes
            eval(['C.',VarName{r1},'.',AttName{r1}{r2},'=attval;']);
        end
        %--------------get values of variables
        val=netcdf.getVar(ncid,vid); % read variable values
        if r1>=1
            val=double(val);
        end
        eval(['C.',VarName{r1},'.val=val;']);
    end
    netcdf.close(ncid);    
    %----modify time
    C.time.val=C.time.val+datenum([double(C.time.base_date),0,0]);
end



end