clear;clc;close all;
prj_names=Get_Projection_Names;

url0='http://spatialreference.org/ref/';

for r1=114:length(prj_names.Id)
    URL=[url0,prj_names.Ref{r1},'/',num2str(prj_names.Id(r1)),'/prj/'];
    fname=['./prj_files/',prj_names.Ref{r1},'.',num2str(prj_names.Id(r1)),'.prj'];
    disp([num2str(r1),':',fname]);
    try
        [filestr,status] = urlwrite(URL,fname);
    catch me
        disp(filestr);       
    end
end
