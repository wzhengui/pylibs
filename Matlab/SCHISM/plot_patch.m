% %----patch for node----------------------------
% ind1=i34==3; ind2=i34==4;
% patch('Faces',elnode(ind1,1:3),'Vertices',[lx,ly],'FaceVertexCdata',dp,'FaceColor','interp','EdgeColor','none'); hold on;
% patch('Faces',elnode(ind2,1:4),'Vertices',[lx,ly],'FaceVertexCdata',dp,'FaceColor','interp','EdgeColor','none'); hold on;
% 
% %----patch for element----------------------------
% ind1=i34==3; ind2=i34==4;
% elnodet=elnode; fp=i34==3; elnodet(fp,4)=1; 
% dpe=dp(elnodet); dpe(fp,4)=0; dpe=sum(dpe,2)./i34;
% patch('Faces',elnode,'Vertices',[lx,ly],'Cdata',dpe,'Facecolor','flat','EdgeColor','none'); hold on;
%
%%-----plot bndinfo-----------------------------------
% for br1=1:bndinfo.nob
%     plot(lx(bndinfo.iobn{br1}),ly(bndinfo.iobn{br1}),'k'); hold on;
% end
% for br1=1:bndinfo.nlb
%     plot(lx(bndinfo.ilbn{br1}),ly(bndinfo.ilbn{br1}),'k'); hold on;
% end