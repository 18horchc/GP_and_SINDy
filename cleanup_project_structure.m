%% Project Structure Cleanup Script
% This script helps clean up the nested GP_and_SINDy/GP_and_SINDy structure
% 
% Current structure:
%   C:\Users\chorc\GP_and_SINDy\          (outer folder)
%     - GP_and_SINDy\                     (inner folder - git repo is here)
%       - .git\
%       - [all project files]
%     - [some duplicate files]
%
% Target structure:
%   C:\Users\chorc\GP_and_SINDy\          (git repo root)
%     - .git\
%     - [all project files]

fprintf('=== Project Structure Cleanup ===\n\n');

% Get current directory
currentDir = pwd;
fprintf('Current directory: %s\n', currentDir);

% Check if we're in the inner folder (git repo)
[status, result] = system('git rev-parse --show-toplevel');
if status == 0
    gitRoot = strtrim(result);
    fprintf('Git root: %s\n', gitRoot);
    
    % Check if git root is the inner nested folder
    if contains(gitRoot, [filesep 'GP_and_SINDy' filesep 'GP_and_SINDy'])
        fprintf('\nDetected nested structure. Git repo is in inner folder.\n');
        fprintf('To clean up, you have two options:\n\n');
        fprintf('OPTION 1: Keep current structure (safer)\n');
        fprintf('  - Git repo stays in inner folder\n');
        fprintf('  - Open MATLAB project from: %s\n', gitRoot);
        fprintf('  - All files are already in the git repo\n\n');
        
        fprintf('OPTION 2: Move git repo to outer folder (cleaner)\n');
        fprintf('  - This requires manual steps (see CLEANUP_PLAN.md)\n');
        fprintf('  - More complex but results in cleaner structure\n\n');
        
        fprintf('RECOMMENDATION: Use Option 1 for now.\n');
        fprintf('Open the MATLAB project file at: %s%sGP_and_SINDy.prj\n', ...
            gitRoot, filesep);
    else
        fprintf('\nGit repo is already at the root level. Structure looks good!\n');
    end
else
    fprintf('Warning: Not in a git repository!\n');
end

fprintf('\n=== Cleanup Complete ===\n');

