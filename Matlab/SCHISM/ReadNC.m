function C=ReadNC(fname)
%%---automatic read all information in netcdf files-------------
%eg. C=ReadNC('COSINE.nc');

N=ncinfo(fname);
%---read dimension informaton-----
for r1=1:length(N.Dimensions)
    C.dimname{r1}=N.Dimensions(r1).Name;
    C.dims(r1)=N.Dimensions(r1).Length;
end

%---read variables------------
ncid=netcdf.open(fname,'NC_NOWRITE'); %open the file
for r1=1:length(N.Variables)    
    C.Var{r1,1}=N.Variables(r1).Name;
    
    %---read variable dimension and size
    for r2=1:length(N.Variables(r1).Dimensions)
        Ci.dimname{r2}=N.Variables(r1).Dimensions(r2).Name;
        Ci.dims(r2)=N.Variables(r1).Size(r2);
    end
    
    %--read variable attributes--------------------
    for r2=1:length(N.Variables(r1).Attributes)
        namei=N.Variables(r1).Attributes(r2).Name;
        if strcmp(namei(1),'_')
            namei(1)=[];
        end
        valuei=N.Variables(r1).Attributes(r2).Value;
        if strcmp(class(valuei),'char')        
            valuei=replace(valuei,'''','''''');            
            eval(['Ci.',namei,'=''',valuei,''';']);   
        elseif strcmp(class(valuei),'int32')
            eval(['Ci.',namei,'=valuei;']);   
        else
            eval(['Ci.',namei,'=valuei;']);               
        end
    end
    
    %--read variable values----
    vid=netcdf.inqVarID(ncid,C.Var{r1});
    val=netcdf.getVar(ncid,vid);
    
    %----processing variable values---------
    if strcmp(N.Variables(r1).Datatype,'char')
        for r2=1:Ci.dims(1)
            pval{r2,1}=strtrim(val(r2,:));
        end
    elseif strcmp(N.Variables(r1).Datatype,'single')
        pval=double(val);
    elseif strcmp(N.Variables(r1).Datatype,'int16')
        pval=double(val);
    elseif strcmp(N.Variables(r1).Datatype,'double')
            pval=val;
    end
    Ci.val=pval;
    
    if isfield(Ci,'FillValue')
        fp=Ci.val==double(Ci.FillValue);
        Ci.val(fp)=nan;
    end
    if isfield(Ci,'miss_value')
        fp=Ci.val==double(Ci.miss_value);
        Ci.val(fp)=nan;
    end
    if isfield(Ci,'scale_factor')
        Ci.val=Ci.val*double(Ci.scale_factor);
    end
    if isfield(Ci,'add_offset')
        Ci.val=Ci.val+double(Ci.add_offset);
    end
    
    tmp=str2num(C.Var{r1}(1));
    if ~isempty(tmp)&isnumeric(tmp)
        C.Var{r1}=['A_',C.Var{r1}];
    end
    eval(['C.',C.Var{r1},'=Ci;']);
    clear Ci;
    
    %--get USta----------
    if strcmp(C.Var{r1},'Station')
        flag=0;
        for r2=1:length(pval)
            if r2==1
                flag=flag+1;
                USta{flag,1}=pval{r2};
            end
            if ~strcmp(USta{flag,1},pval{r2})
                flag=flag+1;
                USta{flag,1}=pval{r2};
            end                
        end
        C.USta=USta;
    end
end
netcdf.close(ncid);

return;
end
