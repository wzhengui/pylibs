function eflag=convert_matfile_format(fname,fn)

eflag=0;

load(fname);
MatFileInfo=whos(matfile(fname));

savestr=['save -v7 ',fn];
for r1=1:length(MatFileInfo)
    if ~strcmp(MatFileInfo(r1).class,'double')&~strcmp(MatFileInfo(r1).class,'cell')
        eflag=1;
        return;
    end
    savestr=[savestr,' ',MatFileInfo(r1).name];
end

for r1=1:length(MatFileInfo)
    if strcmp(MatFileInfo(r1).class,'cell')
        Vari=MatFileInfo(r1).name;
        eval(['yi=',Vari,';']);
        fp=cellfun(@isstr,yi); 
        if sum(~fp)~=0
            eflag=2;
            return;
        end 
        UVar=unique(eval(MatFileInfo(r1).name));
        for r2=1:length(UVar)
            if ~isstr(UVar{r2})
                eflag=3;
                return;
            end
        end
        
        eval(['fp=strcmp(',Vari,','''');',Vari,'(fp)={'' ''};'])
    end  
end

eval(savestr);
return;

end