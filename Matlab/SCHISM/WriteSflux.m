function WriteSflux(data,VarName,AttName,fname)
%this script is to write narr files
%C=WriteSflux(data,fname,flag);

ncid=netcdf.create(fname,'CLOBBER');

%difine dimension
for r1=1:length(data.dims);
    dims(r1)=netcdf.defDim(ncid,data.dimname{r1},data.dims(r1));
end

%define variables
for r1=1:length(VarName)
    if strcmp(VarName{r1},'time')
        vid(r1)=netcdf.defVar(ncid,VarName{r1},'float',dims(3));
    elseif strcmp(VarName{r1},'lon')|strcmp(VarName{r1},'lat')
        vid(r1)=netcdf.defVar(ncid,VarName{r1},'float',dims(1:2));
    else
        vid(r1)=netcdf.defVar(ncid,VarName{r1},'float',dims(1:3));
    end
    
end

%--put attribute
for r1=1:length(VarName)
    for r2=1:length(AttName{r1})
        AttVal=eval(['data.',VarName{r1},'.',AttName{r1}{r2}]);
        netcdf.putAtt(ncid,vid(r1),AttName{r1}{r2},AttVal);
    end
end
netcdf.endDef(ncid);

%--put variables---
for r1=1:length(VarName)
    VarVal=eval(['data.',VarName{r1},'.val']);
    netcdf.putVar(ncid,vid(r1),VarVal);
end

netcdf.close(ncid);

end