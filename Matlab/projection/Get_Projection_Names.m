function prj_names=Get_Projection_Names()

fid=fopen('utm_zones.txt');

flag=0;
while(~feof(fid))
    Ti=fgetl(fid);
    if length(Ti)<4||~strcmp(Ti(1:4),'EPSG')
        continue;
    end
    flag=flag+1;
    i1=strfind(Ti,':'); i2=strfind(Ti,'-');
    prj_names.Ref{flag,1}='epsg';
    prj_names.Id(flag,1)=str2num(strtrim(Ti(i1+1:i2-1)));
    prj_names.Name{flag,1}=strtrim(Ti(i2+1:end));    
end
fclose(fid);

end