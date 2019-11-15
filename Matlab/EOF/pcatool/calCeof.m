% [e,pc,expvar] = calCeof(M,N,METHOD)
%   Compute Complex EOF

    function [e,pc,expvar] = calCeof(M,N,METHOD)


% Compute Hilbert transform Mh and create M as: M = M + i*Mh
M = hilbert(M);

% Compute EOFs
 [e,pc,expvar] = caleof(M,N,METHOD);
