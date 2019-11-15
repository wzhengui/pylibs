function [Proj,VarName,Var]=parse_prj_file(fname)
%parse projection files (.prj)

%---read prj text---------
fid=fopen(fname);
prj_text=fgetl(fid);
fclose(fid);

%--parse it---------------
[Proj,VarName,Var]=parse_prj_text(prj_text);

end

function [C,VarName,Var]=parse_prj_text(prj_text)
% prj_text,
%----get index for delimiters---------
ind_fwd=strfind(prj_text,'[');
ind_bwd=strfind(prj_text,']');
ind_comma=strfind(prj_text,',');
ind_qmk=strfind(prj_text,'"');

ind_brkt=sort([ind_fwd,ind_bwd]);

%---get variable names--------------
Var=prj_text(1:ind_fwd(1)-1);
ii=search_ind(ind_qmk,ind_fwd(1),2,2); i1=ii(1); i2=ii(2);
C.Name=prj_text(i1+1:i2-1);
if strcmp(Var,'PARAMETER')
    C.Name=lower(C.Name);
end
VarName{1,1}=Var;
VarName{1,2}=C.Name;

%--determine whether to exit loop------
i1=search_ind(ind_comma,ind_fwd(1),2,1);
i2=search_ind(ind_bwd,ind_fwd(1),2,1);
i3=search_ind(ind_bwd,ind_fwd(1),2,2); i3=i3(end);
if isnan(i1(1))
    return;

end

if isnan(i3)
    ti=strtrim(prj_text(i1(1)+1:i2-1));
    vi=str2num(ti);
    if isempty(vi)
        disp('wrong value: not a number');
        return;
    end
    C.Value=vi;  
    VarName{1,3}=ti;
    return;
end

%-----continue search for fileds---------
i0=search_ind(ind_comma,ind_fwd(1),2,1);
while(~isnan(i0))    
    flag=0;
    while(1)
        flag=flag+1;
        i1=search_ind(ind_fwd,i0,2,flag+1);
        i2=search_ind(ind_bwd,i0,2,flag);
        prj_text_next=prj_text(i0+1:i2(end));
        
        if isnan(i1(end))|i1(end)>i2(end)            
            break;
        end
    end
    
    [Ci,VarNamei,Vari]=parse_prj_text(prj_text_next);
    
    if strcmp(Vari,'PARAMETER')
        eval(['C.',Ci.Name,'=Ci.Value;']);
    else
        eval(['C.',Vari,'=Ci;']);
    end
    len=size(VarName,1);
    for r1=1:size(VarNamei,1)
        VarName{r1+len,1}=[VarName{1},'.',VarNamei{r1,1}];
        VarName{r1+len,2}=VarNamei{r1,2};
        if size(VarNamei,2)==3
            VarName{r1+len,3}=VarNamei{r1,3};
        end
    end
    i0=search_ind(ind_comma,i2(end),2,1);
end

end


function indout=search_ind(ind,origin,dir,n)
%search indexes
%   ind: index series
%   origin: point of search start
%   dir: search direction (1: backward, 2:foward)
%   n: number of outputs points

if dir==1
    fp=ind<origin;
elseif dir==2
    fp=ind>origin;
else
    disp('wrong search direction');
    return;
end

ind_left=sort(ind(fp));
len=length(ind_left);

indout=nan(1,n);
if len>=n
    if dir==1
        indout=ind_left(end-n+1:end);
    elseif dir==2
        indout=ind_left(1:n);
    end
else
    indout(1:len)=ind_left;
end
    
end


% clear;clc;close all;
% fname='./prj_files/epsg.26918.prj';
% 
% % function fout=parse_prj_file(fname)
% %---parse prj text and outputs structured data;
% fid=fopen(fname);
% prj_text=fgetl(fid),
% fclose(fid);
% 
% ind_fwd=strfind(prj_text,'[');
% ind_bwd=strfind(prj_text,']');
% ind_comma=strfind(prj_text,',');
% ind_qmk=strfind(prj_text,'"');
% 
% ind_brkt=sort([ind_fwd,ind_bwd]);
% 
% %---get variable names--------------
% for r1=1:length(ind_fwd)
%     i0=ind_fwd(r1);
%     
%     %---get varname
%     if r1==1
%       ii=0;
%     else
%         ii=search_ind(ind_comma,i0,1,1);
%     end
%     S.VarName{r1,1}=strtrim(prj_text(ii+1:i0-1));
%     
%     %---get fullanme
%     fpf=ind_fwd<i0; fpb=ind_bwd<i0;
%     sb=sum(fpf)-sum(fpb);    
%     FullVarNamei=S.VarName{r1,1};
%     for r2=[sb:-1:1]
%         FullVarNamei=[S.VarName{r2,1},'.',FullVarNamei];
%     end    
%     S.FullVarName{r1,1}=FullVarNamei;
%     
%     %---get Name------
%     ii=search_ind(ind_qmk,i0,2,2); i1=ii(1); i2=ii(2);
%     S.Name{r1,1}=strtrim(prj_text(i1+1:i2-1));
%     
%     %---get Value-----
%     i3=search_ind(ind_comma,i2,2,1);
%     if ~isnan(i3)
%         i4=search_ind(ind_bwd,i2,2,1);
%         i5=search_ind(ind_fwd,i2,2,1);
%         
%         if (isnan(i5)&i4>i3)|(~isnan(i5)&i4>i3&i5>i4)
%             
%             ti=strtrim(prj_text(i3+1:i4-1));
%             vi=str2num(ti);
%             if isempty(vi)
%                 disp('wrong value: not a number');
%                 disp(S.FullVarName{r1});
%                 return;
%             end                            
%             S.Value{r1,1}=vi;
%         end
%         
%         if ~isnan(i5)&i3<i5&i4>i5
%               flag=0;
%               while(1)
%                   flag=flag+1;
%                   i6=search_ind(ind_fwd,i3,2,flag);
%                   i7=search_ind(ind_bwd,i3,2,flag);
%                   if i6>i7|isnan(i6(end))
%                       S.Value{r1,1}=strtrim(prj_text(i3+1:i7));
%                       break;
%                   end
%               end
% 
%         end                       
%     end    
% end

% return;
% 
% %----restructure the data------------
% nloop=0;
% for r1=1:length(S.VarName)
%     fps=strfind(S.FullVarName{r1},'.');
%     nloop=max(nloop,sum(fps));
% end
% 
% for r1=nloop:-1:0
%     for r2=1:length(S.VarName)
%         fps=strfind(S.FullVarName{r2},'.');
%         if sum(fps)~=r1
%             continue;
%         end        
%         data.Name=S.Name{r2};
%         
%         if isempty(S.Value{r2})
%             S.Data{r2,1}=data;
%             continue;
%         end
%         if isnumeric(S.Value{r2})
%             data.Value=S.Value{r2};
%             S.Data{r2,1}=data;
%             continue
%         else
%             for r3=1:length(S.VarName)
%                 fpsi=strfind(S.FullVarName{r3},'.');
%                 if sum(fpsi)~=(r1+1)
%                     continue;
%                 end
%                 eval([]);
%             end            
%         end
%         
% 
%     end
% end