function values = snc_pp_strings(jobj,jdata,shape)
% Post process NC_STRING data into cell arrays.
if isempty(jdata)
    values = {''};
    return;
elseif strcmp(version('-release'),'14')
    % In R14, we must use the 2.2.x release of java.  No access to the
    % "getObject" method.  Assuming a single-valued string.
    values = {char(jobj.getStringValue())};
    return;
end


% Make it a 1D vector until we know how to shape it properly.
values = cell([1 prod(shape)]);


switch class(jdata.getObject(0))
    case 'opendap.dap.DString'
        for j = 1:prod(shape)
            values{j} = char(jdata.getObject(j-1).getValue());
        end
        
    otherwise       
        for j = 1:prod(shape)
            values{j} = jdata.getObject(j-1);
        end
end

% If just a single string, then do not store it as a cell array.
% netcdf-java treats nc_char as just a single-valued case of nc_string, 
% but mexcdf cannot do this.
if numel(values) == 1
    values = values{1};
else
    values = reshape(values, fliplr(shape));
    n = jdata.getShape();
    if n > 1
        values = permute(values, ndims(values):-1:1);
    end
end

