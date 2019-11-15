function hF=GetFigurePosition(hg,PosModal,Space,Structure,flag)
%hF=GetFigurePosition(hg,PosModal,Space,Structure,flag);
%the handle of the figure;
%the size of PosModal: <n,4>
%the size of Space: <VSpace,HSpace>
%the size of Structure: <VNum,HNum);
%if flag==1, automatic assign a PosModel, and Space Based on Structure
%

if nargin==5&flag==1
    dx=1/Structure(2);
    dy=1/Structure(1);
    PosModal=[0.175*dx,1-0.925*dy,0.825*dx,0.8*dy],
    Space=[0.125*dy,0.15*dx],
end

VNum=Structure(1);HNum=Structure(2);
hF=nan(size(PosModal,1),VNum,HNum);
figure(hg);

Modal=nan(1,4);
Modal(1)=min(PosModal(:,1));
Modal(2)=min(PosModal(:,2));
Modal(3)=max(PosModal(:,1)+PosModal(:,3));
Modal(4)=max(PosModal(:,2)+PosModal(:,4));
PosModal(:,1)=PosModal(:,1)-Modal(1);
PosModal(:,2)=PosModal(:,2)-Modal(2);
for r1=1:Structure(1)
    for r2=1:Structure(2)
        VBase=Modal(2)-(r1-1)*(Modal(4)-Modal(2)+Space(1));
        HBase=Modal(1)+(r2-1)*(Modal(3)-Modal(1)+Space(2));
        for r3=1:size(PosModal,1)
            Posi=[HBase+PosModal(r3,1),VBase+PosModal(r3,2),PosModal(r3,3:4)];
            h=axes('position',Posi,'xticklabel',{},'yticklabel',{},...
                'box','on');
            hF(r3,r1,r2)=h;
        end
    end
end
if size(PosModal,1)==1
    hF=squeeze(hF);
end
end
