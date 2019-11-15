function f=daytimelength(Lat,Doy)
%calculate daytime length
%D=DaytimeLength(Lat,Doy)
%Lat: latitude, Doy: (1-365), sunrise=12-D/2, sunset=12+D/2
P=asin(0.39795*cos(0.2163108 + 2*atan(0.9671396*tan(0.00860*(Doy-186)))));
T=(sin(0.8333*pi/180)+sin(Lat*pi/180)*sin(P))/cos(Lat*pi/180)./cos(P);
f=24-(24/pi)*acos(T);
end
