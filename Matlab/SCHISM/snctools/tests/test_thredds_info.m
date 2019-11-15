function test_thredds_info()

fprintf('\t\tTesting THREDDS_INFO ...  ' );

run_motherlode_test;
fprintf('OK\n');




%--------------------------------------------------------------------------
function run_motherlode_test (  )
    
% This catalog should have 3 top-level sub-catalogs.
url = 'http://tashtego.marine.rutgers.edu:8080/thredds/roms/latte/catalog.xml';
info = thredds_info(url);
if ~strcmp(info.name, 'ROMS LaTTE')
    error('failed');
end
if numel(info.dataset) ~= 3
    error('failed');
end
if ~strcmp(info.URL,url)
    error('failed');
end
return

