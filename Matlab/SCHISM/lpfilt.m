function fdata = lpfilt(data,delta_t,cutoff_f)
%  fdata = lpfilt(data,delta_t,cutoff_f)
%
%  Peforms low-pass filter by multiplication in frequency domain.
%  Uses three-point taper (in frequency space) between pass-band and
%  stop band.

%  Coefficients suggested by D. Coats, Battelle, Ventura.
%  Even better coefficients for the taper can be found in:
%  Rabiner, L.R., Gold,B., and McGonegal, C.A. (1980).  An Approach to
%  the approximation problem for nonrecursive digital filters.  IEEE Tran.
%  vol. AU-18(2):83-106.

%  Written by Chris Sherwood, Battelle PNL
%  Last revised: 9/03/89
data=data(:);
n = length(data);
%c = polyfit((1:n)',data,1);
%trend = polyval(c,1:n)';
%data = data-trend;
mn = mean(data);
data = data-mn;
P = fft(data);
N = length(P);
filt = ones(N,1);
k = floor(cutoff_f*N*delta_t);
filt(1:k) = filt(1:k)*1;
filt(k+1) = .715;
filt(k+2) = .24;
filt(k+3) = .024;
filt(k+4:N-(k+4)) = filt(k+4:N-(k+4))*0.;
filt(N-(k+3)) = .024;
filt(N-(k+2)) = .24;
filt(N-(k+1)) = .715;
P = P .* filt;
fdata = real(ifft(P));
fdata = fdata(1:n)+mn;
%fdata = fdata(1:n)+trend;
fdata=fdata(:);
end