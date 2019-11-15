function createLinks(filterspec,startDirName,recurse,linkToDirName)
% Create Windows shell link (.LNK) shortcuts
% 
% SYNTAX:
% createLink(filterspec,startDirName,recurse,linktoDirName)
%     creates Windows shell link (.lnk) shortcuts to files
%     specified by extensions in filterspec, for files in
%     directory startDirName, including or excluding
%     subdirectories, and writing links to directory
%     linkToDirName.
%
% INPUTS:
% filterspec:    a cell array of extensions to include.
%                Specify filterspec as {'*.ext1';'*.ext2';...;'*.extn'}.
%                (NOTE: See help for UIGETFILE for a
%                discussion of filterspec.)
% startDirName:  the top-level directory from which to start
%                searching for files to link. Default = pwd;
% recurse:       Logical (t/f); do you want to recurse (look
%                sub-directories)? Default = false;
% linkToDirName: The name of the directory in which you want
%                to put the shortcuts. Default = pwd.
%
% EXAMPLES:
% Example 1: To create, in the current directory, links to
%            all JPG or TIFF files in the current directory:
% 
%            createLinks({'*.jpg';'*.jpeg';'*.tif';'*.tiff'})
%
% Example 2: To create links to all files of type *.m in
%            directory 'c:\mfiles', recursively, and place
%            the links in directory 'c:\myMfileLinks':
%
%            createLinks({'*.m'},'c:\mfiles',true,'c:\myMfileLinks')
%
% NOTE: To EXTRACT the target name from the links, file
%       getTargetFromLink might be useful. (It is available
%       on the MATLAB Central File Exchange.)
% 
% Written by Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 
% Thanks to Jiro Doke, PhD. for his assistance.
%
% Copyright 2012 MathWorks, Inc.
%
% See Also: getTargetFromLink

if nargin < 1
    error('CREATELINKS: Requires at a minimum the specification of the filterspec.')
end
if nargin < 4
    linkToDirName = pwd;
end
if ~exist(linkToDirName,'dir')
    mkdir(linkToDirName)
end
if nargin < 3
    recurse = false;
end
if nargin < 2
    startDirName = pwd;
end
%
numspecs=length(filterspec);
filterspec = sprintf('%s;',filterspec{:});
keepFilterspec = filterspec;
% numspecs = length(find(filterspec=='*'));

if recurse
    if strcmp(startDirName(end),'\')
        allPaths = genpath(startDirName);
    else
        allPaths = genpath([startDirName,'\']);
    end
else
    allPaths = [startDirName ';'];
end
nPaths = numel(strfind(allPaths,';'));
files = cell(100000,1);
dirs = files;

nFound = 0;
for jj = 1:nPaths
    [pathToken,allPaths] = strtok(allPaths,';');
    for ii = 1:numspecs
        [filterToken,filterspec] = strtok(filterspec,';');
        if strfind(filterToken,'*')
            tmp = dir([pathToken, filesep, filterToken]);
            for kk = 1:numel(tmp)
                nFound = nFound + 1;
                files{nFound} = tmp(kk).name;
                dirs{nFound} = pathToken;
            end
        else            
            nFound = nFound + 1;
            files{nFound} = filterToken;
            dirs{nFound} = pathToken;
        end

    end
    filterspec = keepFilterspec;
end
files(nFound+1:end) = [];
dirs(nFound+1:end) = [];

asvr = actxserver('WScript.Shell');
for ii = 1:numel(files)
    if strcmp(files(ii),'.')|strcmp(files(ii),'..')
        continue;
    end
    [~,fn] = fileparts(files{ii});
    b = asvr.CreateShortcut(fullfile(linkToDirName, [fn,'.lnk']));
    b.TargetPath = fullfile(dirs{ii}, files{ii});
    b.Save();
end
