% %write excel files, format---------------------
% Data=ones(4,3);
% Station={'x1','x2','x3'};
% Time(1:4,1)={'Time'};
% 
% %%--1st------------------
% Table=[{'1st Test'} Station; Time num2cell(Data)];
% xlswrite('my_file',Table,'1st sheet');
% 
% %%--2nd--------------------------------
% xlswrite('my_file',2*Data,'2nd sheet','B2');
% xlswrite('my_file',Station,'2nd sheet','B1');
% xlswrite('my_file',Time,'2nd sheet','A2');
% xlswrite('my_file',{'2nd test'},'2nd sheet','A1');
