function WriteNC(data,fname)
%aumatic write netcdf files for structure Data which has a format from
%ReadNC. 

ncid=netcdf.create(fname,'CLOBBER');

%difine dimension
for r1=1:length(data.dims);
    dims(r1)=netcdf.defDim(ncid,data.dimname{r1},data.dims(r1));
end

VarName=data.Var;
for r1=1:length(VarName)    
    datai=eval(['data.',VarName{r1}]);  
    
    %get dim index
    [tmp,ind1,ind2]=intersect(datai.dimname,data.dimname);
    
    %define variables
    vid(r1)=netcdf.defVar(ncid,VarName{r1},'float',dims(ind2)); 
    
    %put attribute
    AttName=fieldnames(datai);
    for r2=1:length(AttName)
        if strcmp(AttName{r2},'dimname')|strcmp(AttName{r2},'dims')|strcmp(AttName{r2},'val')
            continue;
        end
        AttVal=eval(['datai.',AttName{r2}]);
        netcdf.putAtt(ncid,vid(r1),AttName{r2},AttVal);
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