% init.m — Add all subdirectories of the current folder to MATLAB path
% Cross-platform: works on Windows, macOS, and Linux
function init()
% INIT  Add all subdirectories of the current folder to MATLAB path
%   Works across Windows, macOS, and Linux.

    % Get full path to the folder where this script is located
    projectRoot = fileparts(mfilename('fullpath'));

    % Generate list of all subfolders recursively
    allFolders = genpath(projectRoot);

    % Add to MATLAB path
    addpath(allFolders);

    % Display confirmation
    fprintf('All folders added to path from:\n%s\n', projectRoot);
end