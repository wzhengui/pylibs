%%code for taylor diagram   
% % axes(hF(r1));
% % 
% % [hp ht hax]=taylordiag(stdi,rmsei,corri,...
% %     'tickRMS',[0.5:0.5:2],'titleRMS',1,'showlabelsRMS',1,...
% %     'widthRMS',1,'colRMS','g','tickRMSangle',115,...
% %     'tickSTD',[0:0.5:2],'limSTD',2,'styleSTD',':','widthSTD',1,'colCOR','b',...
% %     'tickCOR',[0.1:0.1:0.9,0.95,0.99],'showlabelsCOR',1,'titleCOR',1);
% % 
% % set(hax(1).handle,'position',[-0.2,1.15,-1],'FontSize',12,'string','Normalized Standard Deviation');
% % set(hp(1),'color','k','MarkerSize',30);set(ht(1),'string',{});
% % 
% % lstr={' ',' ','Obs'};
% % for r2=1:length(Station)
% %     set(hp(r2+1),'Marker',Marker{r2},'color',color{r2},'MarkerSize',8,'MarkerFacecolor',color{r2});
% %     set(ht(r2+1),'String',{});
% %     lstr{r2+3}=Station{r2};
% % end
% % if r1==1
% %     hl=legend(lstr);
% %     set(hl,'position',[0.6513    0.1    0.0493    0.2808]);
% % end
% % 
% % text(0.7,2.2,Var{r1},'FontSize',14,'FontWeight','bold');