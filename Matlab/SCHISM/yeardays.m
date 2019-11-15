function n=yeardays(Yr)
%return number of days in a year;
%n=yeardays(Year);
n=datenum(Yr+1,1,1)-datenum(Yr,1,1);
end