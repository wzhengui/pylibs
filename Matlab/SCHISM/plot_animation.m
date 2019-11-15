% %how to save an animation
% 
% for r1=n1:n2
%     %plot---
%     pause(0.01);
%     
%     
%     %--save-----    
%     frame = getframe(gcf);
%     flag=flag+1;
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);    
%     fname='FlowDist.gif';
%     if flag==1;
%         imwrite(imind,cm,fname,'gif','Loopcount',inf,'delaytime',1);
%     else
%         imwrite(imind,cm,fname,'gif','WriteMode','append','delaytime',1);
%     end
% end