function ensure_gpml_path(gpml_path)
%ENSURE_GPML_PATH  Add GPML toolbox to path. Errors if not found.
%
%   ensure_gpml_path()
%   ensure_gpml_path(gpml_path)
%
%   If gpml_path given, uses it. Otherwise tries common locations.
    if nargin < 1, gpml_path = []; end

    if ~isempty(gpml_path)
        gpml_root = gpml_path;
    else
        script_dir = fileparts(mfilename('fullpath'));
        cand = {fullfile(matlabroot, 'toolbox', 'gpml');
                'C:\gpml-matlab';
                'C:\Program Files\MATLAB\gpml-matlab';
                fullfile(getenv('USERPROFILE'), 'gpml-matlab');
                fullfile(script_dir, '..', 'gpml-matlab');
                fullfile(script_dir, '..', 'gpml-matlab-master');
                'C:\Users\chorc\GP_and_SINDy\gpml-matlab-master'};
        gpml_root = [];
        for i = 1:length(cand)
            if ~isempty(cand{i}) && exist(fullfile(cand{i}, 'gp.m'), 'file')
                gpml_root = cand{i};
                break;
            end
        end
    end
    if isempty(gpml_root) || ~exist(fullfile(gpml_root, 'gp.m'), 'file')
        error('ensure_gpml_path:not_found', ...
            'GPML not found. Download from https://gitlab.com/hnickisch/gpml-matlab and pass opts.gpml_path.');
    end
    addpath(gpml_root);
    if exist(fullfile(gpml_root, 'startup.m'), 'file')
        run(fullfile(gpml_root, 'startup.m'));
    end
end
