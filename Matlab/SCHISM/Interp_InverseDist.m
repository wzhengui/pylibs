function f=Interp_InverseDist(X,V,factor,iCheck)
%Input: X(npt,3) %(number of obs, xyz)
%        V(npt,2) %(number of pts, xy)
%        factor: distance factor.
%        iCheck: plot data and result for double check
%
%eg.  f=Interp_InverseDist(X,V,2)   %
%     f=Interp_InverseDist(X,V,1,1)  %plot

lxi=repmat(V(:,1),1,size(X,1));
lyi=repmat(V(:,2),1,size(X,1));
xi=repmat(X(:,1)',size(V,1),1);
yi=repmat(X(:,2)',size(V,1),1);
dist=1./sqrt((lxi-xi).^2+(lyi-yi).^2);
dist=dist.^(factor);
f=dist*X(:,3)./sum(dist,2);

if nargin==4&iCheck==1
    figure;
    subplot(2,1,1);
    scatter(X(:,1),X(:,2),40,X(:,3),'filled');
    subplot(2,1,2);
    scatter(V(:,1),V(:,2),40,f,'filled');
end
end