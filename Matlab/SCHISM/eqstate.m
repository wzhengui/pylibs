function rho=eqstate(T,S,P)
%--------------------------------------------------------------------
%written by Zhengui Wang on Jun 19, 2017
%calculate seawater density (Millero and Poisson 1981, Gill, 1982)
%Input:
%    1)T: temperature in [C]
%    2)S: salinity in [PSU]
%    3)P: pressure [bars] (P=0: at 1 atm. pressure)
%Output:
%    1)rho: seawater density in [Kg/m3]
%--------------------------------------------------------------------
%   implicit none
%   integer,parameter :: rkind=8
%   real(rkind), intent(in) :: T,S,P
%   real(rkind), intent(out) :: rho
%
%   %local
%   real(rkind) :: T2,T3,T4,T5,S05,S15,S2,P2
%   real(rkind) :: rho_pw,rho_st,A,B,C,K_pw,K_st,K_stp

%pre_calculation
T2=T*T;
T3=T^3;
T4=T^4;
T5=T^5;
S05=sqrt(S);
S15=S*S05;
S2=S*S;
P2=P*P;

%pure water S=0,at 1 atm.
rho_pw=999.842594d0+6.793952d-2*T-9.095290d-3*T2+1.001685d-4*T3 ...
    -1.120083d-6*T4+6.536332d-9*T5;

%density with Salinity
A=0.0; B=0.0; C=0.0;
if S~=0.0
    A=8.24493d-1-4.0899d-3*T+7.6438d-5*T2-8.2467d-7*T3+5.3875d-9*T4;
    B=-5.72466d-3+1.0227d-4*T-1.6546d-6*T2;
    C=4.8314d-4;
end %S
rho_st=rho_pw+A*S+B*S15+C*S2;

%pressure not zero
K_stp=0.0;
if P~=0.0
    K_pw=19652.21d0+148.4206d0*T-2.327105d0*T2+1.360477d-2*T3-5.155288d-5*T4;
    K_st=K_pw;
    if S~=0.0
        K_st=K_st+S*(54.6746d0-0.603459d0*T+1.09987d-2*T2-6.1670d-5*T3) ...
            +S15*(7.944d-2+1.6483d-2*T-5.3009d-4*T2);
    end %S
    K_stp=K_st+P*(3.239908d0+1.43713d-3*T+1.16092d-4*T2-5.77905d-7*T3) ...
        +P*S*(2.2838d-3-1.0981d-5*T-1.6078d-6*T2)+1.91075d-4*P*S15 ...
        +P2*(8.50935d-5-6.12293d-6*T+5.2787d-8*T2) ...
        +P2*S*(-9.9348d-7+2.0816d-8*T+9.1697d-10*T2);
    rho=rho_st/(1.d0-P/K_stp);
else
    rho=rho_st;
end %P

return
end

