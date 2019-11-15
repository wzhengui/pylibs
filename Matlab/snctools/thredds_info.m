function info = thredds_info(url)
% THREDDS_INFO  Retrieve metadata from a thredds catalog.
%   INFO = THREDDS_INFO(XML_URL) recursively dumps metadata from a thredds
%   catalog URL.  The URL must be for the catalog XML file.
%
% Example:
%   url = 'http://tashtego.marine.rutgers.edu:8080/thredds/roms/latte/catalog.xml';
%   info = thredds_info(url);


import thredds.catalog.*
x = InvCatalogFactory('test',true);
catalog = x.readXML(url);

info = catalog_info(catalog);

info.URL = url;

%--------------------------------------------------------------------------
function info = catalog_info(catalog)

z = catalog.getDatasets();
num_elements = z.size();
if num_elements == 0
    
    info = get_dataset_info(catalog);
    
else
    
    info.name = char(catalog.getName());
    
    for j = 0:num_elements-1
        sub_catalog = z.get(j);
        x = catalog_info(sub_catalog);
        info.dataset(j+1) = x;
    end
end

%--------------------------------------------------------------------------
function info = get_dataset_info(dataset)

urlpath = dataset.getUrlPath();
if isempty(urlpath)
    return;
end

info = struct('name',[],'service',[],'time_coverage',[]);

opendap_url = construct_opendap_url(dataset);

info.name = char(dataset.getFullName());
info.service.opendap = opendap_url;

time_coverage = dataset.getTimeCoverage();
if isempty(time_coverage)
    info.time_coverage.start = [];
    info.time_coverage.stop = [];
else
    start = time_coverage.getStart();
    s = char(start.toString());
    info.time_coverage.start = datenum(s(1:end-1));
    
    stop = time_coverage.getEnd();
    s = char(stop.toString());
    if strcmp(s,'present')
        info.time_coverage.stop = now;
    else
        info.time_coverage.stop = datenum(s(1:end-1));
    end
end
    
    
%--------------------------------------------------------------------------
function url = construct_opendap_url(dataset)
access = dataset.getAccess();
n = access.size();

for j = 0:n-1
    atype = access.get(j);
    service = atype.getService();
    if strcmp(char(service.getServiceType()),'OPENDAP')
        url = atype.getStandardUri();
        url = char(url.toString());
        return
    end
end
