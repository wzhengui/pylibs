% C = MAP2MAT(F,D) Reshaping matrix
%
% => Creation of a 'map' matrix C(TIME,LON,LAT) from D(TIME,PT)
% following mask F(LON,LAT).
% F is a matrix contenaing 1 where you would like to keep the point and
%  0 elsewhere (see mat2map).
% 
% Rq: No check is done about the input.
%
% See also: map2mat
%================================================================

% March 2004
% gmaze@univ-brest.fr

function [C] = mat2map(F,D);

% Get dimensions
[nolon nolat] = size(F);
[time npt] = size(D);

% So output 'map' matrix has the form:
C = zeros(time,nolon,nolat);

% Variables
nul = NaN.*ones(time,1);
ipt = 0 ;


% 'mapping' :
  for iy=1:nolat
      for ix=1:nolon
         if F(ix,iy)>0
            ipt = ipt + 1;
            C(:,ix,iy) = D(:,ipt);
         else
            C(:,ix,iy) = nul;
         end %if
      end %for ix
  end %for iy

