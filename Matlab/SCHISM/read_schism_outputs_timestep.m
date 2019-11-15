function [data,ts]=read_schism_outputs_timestep(h,n,fn)
%read time step from schism outputs
%[data,ts]=read_schism_outputs_timestep(h,n,fn)
%
% input 
%      h: header from read_schism_outputs_header
%      n: array of timestep numbers
%      fn: optional argument for reading from filename fn instead of h.fname
% output


n=round(n);
if max(n)>h.nrec
    error('time step out of range');
end

fname=h.fname;
if nargin==3&~isequal(fname,fn)
    fname=fn;
end

fid=fopen(fname);

for r1=1:length(n)
    [datai,tsi]=readTS(fid,h,n(r1));
    data{r1}=datai;
    ts{r1}=tsi;    
end

fclose(fid);
end

function [Data,T]=readTS(fid,h,n)
%read certain time step results
fseek(fid,h.DataStartPos+h.StepSize*(n-1),'bof');
T.time=fread(fid,1,'float32');
T.it=fread(fid,1,'int32');
T.eta=fread(fid,h.hgrid.np,'float32');
datai=fread(fid,[h.ivs,h.GridSize],'float32');

if h.i23d==3
    if h.ivs==2
        Data(1,:,:)=reshape(datai(1,h.index),h.vgrid.nvrt,h.hgrid.np);
        Data(2,:,:)=reshape(datai(2,h.index),h.vgrid.nvrt,h.hgrid.np);
    elseif h.ivs==1
        Data(:,:)=reshape(datai(1,h.index),h.vgrid.nvrt,h.hgrid.np);
    else
        disp('wrong h.ivs');
        return;
    end
elseif h.i23d==2
    if h.ivs==2
        Data(1,:)=reshape(datai(1,h.index),h.hgrid.np,1);
        Data(2,:)=reshape(datai(2,h.index),h.hgrid.np,1);
    elseif h.ivs==1
        Data=reshape(datai(1,h.index),1,h.hgrid.np);
    else
        disp('wrong h.ivs');
        return;
    end
else
    disp('wrong h.i23d');
    return;
end

return;
end

