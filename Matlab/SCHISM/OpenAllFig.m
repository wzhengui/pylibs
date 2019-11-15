function OpenAllFig(fname,flag)
%open all figures in directory or file
if isdir(fname)
    files=dir([fname,'/*.fig']);
    for r1=1:length(files)
       Fig{r1}=[fname,'/',files(r1).name];
    end   
else
    Fig{1}=fname;
end

%----open all figures---------
for r1=1:length(Fig)
        openfig(Fig{r1});
        if nargin==2&flag==1
            set(gcf,'position',get(0,'ScreenSize'));
        end
end

return; 
end

