%%%------------------------------------------------------------------------
%%read SCHISM station.bp file
%%usage:
%%1) BP=read_schism_bpfile(fname)
%%2) with station name: BP=read_schism_bpfile(fname,1,delimiter)
%%%------------------------------------------------------------------------
function BP=read_schism_bpfile(fname,flag,delimiter)

if nargin~=3
    delimiter='!';
end

%---read bp info-----------
fid=fopen(fname);
fgetl(fid);

BP.npt=str2num(fgetl(fid));
for r1=1:BP.npt
    Ti=fgetl(fid);
    
    %--with station info
    if nargin>=2&flag~=0
        ind=strfind(Ti,delimiter);
        if isempty(ind)
            BP.Station{r1,1}='';
            Ti=str2num(Ti);
        else
            BP.Station{r1,1}=Ti(ind+1:end);
            Ti=str2num(Ti(1:ind-1));
        end
       
    else
        Ti=str2num(Ti);
    end
    BP.id(r1,1)=Ti(1);
    BP.lx(r1,1)=Ti(2);
    BP.ly(r1,1)=Ti(3);
    BP.dp(r1,1)=Ti(4);
end
fclose(fid);
end