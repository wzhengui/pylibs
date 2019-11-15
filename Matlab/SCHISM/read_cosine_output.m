function C=read_cosine_output(fname,Var,sf)
%read CoSiNE model outputs
%eg. C=read_cosine_output('cstation',{'S2'},'cstation.in')

% fname='run2r/cstation.out';
% Var={'S2','ph','T'};
% sf='./run2r/cstation.in';


Tracer_name={'T','S','NO3','SiO4','NH4','S1','S2','Z1','Z2','DN','DSi','PO4','DO','CO2','ALK',...
    'NPS1','RPS1','NPS2','RPS2','MTS1','MTS2','MTZ1','MTZ2','EXZ1','EXZ2','GS1Z1','GS2Z2','GZ1Z2',...
    'GDNZ2','GTZ2','SKS2','SKDN','SKDSi','Nit','MIDN','MIDSi','pnh4s1','pnh4s2','pih1','pih2',...
    'fS1','fS2','bfNO3S1','bfNH4S1','bfNO3S2','bfNH4S2','fNO3S1','fNH4S1','fNO3S2','fNH4S2',...
    'fPO4S1','fPO4S2','fCO2S1','fCO2S2','fSiO4S2','o2flx','co2flx','OXR','ADPT','ph'};
[CVar,ind1,ind2]=intersect(Var,Tracer_name);
if ~isequal(sort(CVar),sort(Var))
    disp('unknown variable');
end

%read station info
BP=read_schism_bpfile(sf,1);

%read data
Data=load(fname);

C.time=Data(:,1);
for r1=1:length(CVar)
    i1=(ind2(r1)-1)*BP.npt+2; 
    i2=i1+BP.npt-1;
    eval(['C.',CVar{r1},'=Data(:,i1:i2);']);
end
C.StaInfo=BP;

return;
end