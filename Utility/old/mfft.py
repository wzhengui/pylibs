#!/usr/bin/evn python3
from pylib import *
from scipy.fftpack import fft

def mfft(xi,dt):
    #input
    #xi: time series
    #dt: interval
    #
    #output
    #perid[period],afx[amplitude],pfx[phase]
    N=xi.size;
    fx=fft(xi);
    afx=abs(fx[1:N//2])*2.0/N;
    pfx=angle(fx[1:N//2]);
    period=dt*N/arange(1,N//2);
    return period,afx,pfx

if __name__=="__main__":
    pass;

#    plt.close('all')
#    T=10; dt=0.01; N=T/dt;
#    x=linspace(0.0, T, N);
#    y=4*cos(2.0*pi*(x-0.3)/0.5)+2*cos(2.0*pi*(x-0.4)/1.0)+4*cos(2.0*pi*(x-0.5)/2.0)
#    f,a,p=mfft(y,dt)
#
#    subplot(2,1,1)
#    plot(x,y,'k-')
#    subplot(2,1,2)
#    plot(f,a,'k.',ms=20)
#    setp(gca(),'xlim',[0,5])
