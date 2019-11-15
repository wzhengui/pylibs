function h=read_schism_outputs_header(fname)
% h=read_schism_outputs_header(fname)
% read the header of schism binary output based on sz_readHeader
% Data Format V5.00
% Written by Zhengui Wang on 08/09/2016
h.fname=fname;
fid=fopen(fname);

%----read header--
h=readHeader(fid,h);

%----vgrid info.----
h=readVgrid(fid,h);

%----hgrid info.----
h=readHgrid(fid,h);

%----stating position----
h.DataStartPos=ftell(fid);

%----compute step size----
h=ComputeStepSize(h);

%----compute index----
h=ComputeIndex(h);


fclose(fid);
end

function h=readHeader(fid,h)

h.DataFormat   = char(fread(fid,48,'uchar')');
h.Version      = char(fread(fid,48,'uchar')');
h.StartTime    = char(fread(fid,48,'uchar')');
h.VarName      = char(fread(fid,48,'uchar')');
h.VarDimension = char(fread(fid,48,'uchar')');
h.nrec         = fread(fid,1,'int32');
h.dt           = fread(fid,1,'float32');
h.skip         = fread(fid,1,'int32');
h.ivs          = fread(fid,1,'int32');
h.i23d         = fread(fid,1,'int32');

end

function h=readHgrid(fid,h)
h.hgrid.StartPos=ftell(fid);
h.hgrid.np=fread(fid,1,'int32');
h.hgrid.ne=fread(fid,1,'int32');

sp=ftell(fid);
Ti=fread(fid,[4,h.hgrid.np],'float32');
h.hgrid.x=Ti(1,:)';
h.hgrid.y=Ti(2,:)';
h.hgrid.dp=Ti(3,:)';

fseek(fid,sp,'bof');
Ti=fread(fid,[4,h.hgrid.np],'int32');
h.hgrid.kbp00=Ti(4,:)';

h.hgrid.i34=nan(h.hgrid.ne,1);
h.hgrid.elnode=nan(h.hgrid.ne,4);
for r1=1:h.hgrid.ne
    i34=fread(fid,1,'int32');
    h.hgrid.i34(r1)=i34;    
    h.hgrid.elnode(r1,1:i34)=fread(fid,i34,'int32');
end
end

function h=readVgrid(fid,h)
h.vgrid.StartPos=ftell(fid);
h.vgrid.nvrt    = fread(fid,1,'int32');
h.vgrid.kz      = fread(fid,1,'int32');
h.vgrid.h0      = fread(fid,1,'float32');
h.vgrid.hs      = fread(fid,1,'float32');
h.vgrid.hc      = fread(fid,1,'float32');
h.vgrid.theta_b = fread(fid,1,'float32');
h.vgrid.theta_f = fread(fid,1,'float32');             
h.vgrid.ztot    = fread(fid,h.vgrid.kz-1,'float32');
h.vgrid.sigma   = fread(fid,h.vgrid.nvrt-h.vgrid.kz+1,'float32');
end


function h=ComputeStepSize(h)
if h.i23d==3
    nlevel=h.vgrid.nvrt-max(1,h.hgrid.kbp00)+1;
    h.GridSize=sum(nlevel);
elseif h.i23d==2
    h.GridSize=h.hgrid.np;
else
    error('wrong i23d');
end
h.StepSize=4*(2+h.hgrid.np+h.GridSize*h.ivs);
end

function h=ComputeIndex(h)
nvrt=h.vgrid.nvrt;
kbp=h.hgrid.kbp00;
np=h.hgrid.np;

fp=kbp==0; kbp(fp)=1;
%---incremental index for every column--
if h.i23d==3
    %-----index base----    
    nlevel=nvrt-kbp+1;
    cum_nlevel=cumsum(nlevel);
    index0=repmat(cum_nlevel',nvrt,1);
    
    index=zeros(nvrt,np);
    indexN=zeros(nvrt,np);
    for r1=1:nvrt
        fp=kbp>r1;
        index(r1,fp)=-nvrt+kbp(fp);
        indexN(r1,fp)=nan;
        
        fp=kbp<=r1;
        index(r1,fp)=r1-nvrt;
    end
    h.index=index0(:)+index(:);
    h.indexN=isnan(indexN(:));
elseif h.i23d==2
    h.index=1:np;
else
    error('unknown i23d');
end
end