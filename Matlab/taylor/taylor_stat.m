function stat=taylor_stat(obs,model)
%%--calculate statistics for taylor diagram-----
%%--stat=taylor_stat(obs,model)--------
%%--stat=[std_obs,std_model,rmsd,corr] with definitions below:
%
%       - The STANDARD DEVIATION is computed as:
%                                 /  sum[ {C-mean(C)} .^2]  \
%                       std = sqrt|  ---------------------  |
%                                 \          N              /
%
%       - The CENTERED ROOT MEAN SQUARE DIFFERENCE is computed as:
%                                  /  sum[  { [C-mean(C)] - [Cr-mean(Cr)] }.^2  ]  \
%                       rmsd = sqrt|  -------------------------------------------  |
%                                  \                      N                        /
%
%       - The CORRELATION is computed as:
%                             sum( [C-mean(C)].*[Cr-mean(Cr)] )
%                       corr = ---------------------------------
%                                     N*STD(C)*STD(Cr)

mobs=mean(obs); mmodel=mean(model);
obsi=obs-mobs; modeli=model-mmodel;
std_obs=sqrt(mean(obsi.^2));
std_model=sqrt(mean(modeli.^2));
rmsei=sqrt(mean((obsi-modeli).^2));
corri=mean(obsi.*modeli)/std_obs/std_model;

stat=[std_obs,std_model,rmsei,corri];
end