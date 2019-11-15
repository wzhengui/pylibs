function write_schism_hgrid(ne,np,lx,ly,dp,i34,elnode,fname,note)
%write schism hgrid.gr3
%format:  write_schism_hgrid(ne,np,lx,ly,dp,i34,elnode,fname,note);
%eg. 
%    1) write_schism_hgrid(ne,np,lx,ly,dp,i34,elnode,'hgrid.gr3')
%    2) write_schism_hgrid(ne,np,lx,ly,dp,i34,elnode,'hgrid.gr3','add annotation')


indp=[1:np]';
lx=lx(:); ly=ly(:);dp=dp(:);i34=i34(:);
%----write hgrid.gr3----
fid=fopen(fname,'w+');
if nargin==8
    fprintf(fid,'%s \n',fname);
elseif nargin==9
    fprintf(fid,'%s \n',note);
else
    disp('incorrect nargin');
    return;    
end

fprintf(fid,'%d %d\n',ne,np);
fprintf(fid,'%d %14.6f %14.6f %14.7e\n',[indp,lx,ly,dp]');
for r1=1:ne
    if i34(r1)==3
        fprintf(fid,'%d %d %d %d %d\n',r1,i34(r1),elnode(r1,1:3));
    elseif i34(r1)==4
        fprintf(fid,'%d %d %d %d %d %d\n',r1,i34(r1),elnode(r1,1:4));
    else
        disp('wrong i34 value');
        return;
    end
end
fclose(fid);

end