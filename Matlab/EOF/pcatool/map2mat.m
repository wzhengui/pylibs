% D = MAP2MAT(F,C) Reshaping matrix
%
% => Conversion of a 'map' matrix C(TIME,LON,LAT) into a D(TIME,PT) matrix
% under the mask F(LON,LAT).
% F is a matrix contenaing 1 where you would like to keep the point and
%  0 elsewhere.
%
% Rq: No check is done about the input.
%
% See also: mat2map
%================================================================

% March 2004
% gmaze@univ-brest.fr

function [D] = map2mat(F,C);

% Get dimensions
[tps nolon nolat] = size(C);

% So output matrix will be:
D = zeros(tps,nolon*nolat);

% point indice
ipt = 0;

% 'Un-mapping' :
for iy=1:nolat
  for ix=1:nolon
      if F(ix,iy)>0
         ipt = ipt + 1;
         D(:,ipt)=squeeze(C(:,ix,iy));
      end % if
  end % for
end %for

% OUTPUT:
D = squeeze(D(:,1:ipt));


