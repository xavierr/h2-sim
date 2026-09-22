function startupH2sim

% This startup file set up the MATLAB path
%
%% We use first `MRST <https://github.com/SINTEF-AppliedCompSci/MRST>`_ setup for MRST modules.
% The source code for MRST is synchronized to H2Sim using a git submodule (see ``submodules/mrst``).
%

fprintf('\n  /-\\\n | + | H2Sim\n  \\-/\n\n');
fprintf('Welcome to the H2 Storage Modeling Toolbox (H2sim)!\n');
fprintf('H2sim is based on MRST, which will now be initialized.\n\n');

rootdirname = fileparts(mfilename('fullpath'));

%% Make sure the git submodules (MRST, PhreeqcMatlab) are present, downloading them on demand
%
% If a submodule directory is missing (e.g. the repository was cloned without
% ``--recurse-submodules``), it is fetched automatically here using
% ``git submodule update --init``. This always checks out the exact commit that is pinned for
% that submodule in the H2sim repository, so newer, possibly incompatible, upstream changes to
% MRST/PhreeqcMatlab never affect H2sim until that pin is deliberately updated.

ensureSubmodule(rootdirname, 'mrst');
ensureSubmodule(rootdirname, 'PhreeqcMatlab');

run(fullfile(rootdirname, 'submodules', 'mrst', 'startup'));

dirnames = {'h2-biochem', ...
            'h2-store', ...
            'utils'};

for ind = 1 : numel(dirnames)
    dirname = fullfile(rootdirname, dirnames{ind});
    addpath(genpath(dirname));
end

%% Install Phreeqcs

run(fullfile(rootdirname, 'submodules', 'PhreeqcMatlab', 'startup'));

%% Octave requires some extra functionality
if mrstPlatform('octave')

    % Octave MRST settings
    run(fullfile(rootdirname, 'submodules', 'mrst', 'core', 'utils', 'octave_only', 'startup_octave.m'));

    % Disable warnings
    warning('off', 'Octave:possible-matlab-short-circuit-operator');
    warning('off', 'Octave:data-file-in-path');

    % Install package for json files for older octave
    if compare_versions(version, "6.4", "<=")
        try
            pkg load jsonstuff
        catch
            fprintf('Trying to install jsonstuff...\n');
            pkg install "https://github.com/apjanke/octave-jsonstuff/releases/download/v0.3.3/jsonstuff-0.3.3.tar.gz"
            pkg load jsonstuff
        end
    end

    % For running Julia from Octave, a tcp client such as
    % https://gnu-octave.github.io/packages/instrument-control/ is
    % needed
    try
        pkg load instrument-control
    catch
        fprintf('Trying to install instrument-control...\n');
        pkg install "https://downloads.sourceforge.net/project/octave/Octave%20Forge%20Packages/Individual%20Package%20Releases/instrument-control-0.9.1.tar.gz"
        pkg load instrument-control
    end

end

mrstModule add compositional ad-blackoil ad-core ad-props deckformat

end

function ensureSubmodule(rootdirname, name)
% Clone the git submodule ``name`` (under ``submodules/<name>``) if it is not already present,
% checking out the commit that is pinned for it in the H2sim repository (see ``.gitmodules``).

    submodulepath = fullfile(rootdirname, 'submodules', name);

    if exist(fullfile(submodulepath, 'startup.m'), 'file')
        % Already downloaded
        return
    end

    if ~isfolder(fullfile(rootdirname, '.git'))
        error('h2sim:missingGitRepo', ['Cannot auto-download the ''%s'' submodule because\n  %s\n' ...
                            'is not a git checkout. Please install H2sim using ' ...
                            '''git clone --recurse-submodules'' as described in the documentation ' ...
                            '(doc/installation.rst), or download ''%s'' manually into ''%s''.'], ...
              name, rootdirname, name, submodulepath);
    end

    fprintf('Submodule ''%s'' was not found, downloading it now (this may take a while)...\n', name);

    relpath = strrep(fullfile('submodules', name), filesep, '/');
    cmd = sprintf('git -C "%s" submodule update --init --recursive -- "%s"', rootdirname, relpath);
    [status, cmdout] = system(cmd);

    if status ~= 0 || ~exist(fullfile(submodulepath, 'startup.m'), 'file')
        error('h2sim:submoduleSetupFailed', ['Failed to download the ''%s'' submodule automatically.\n' ...
                            'Command: %s\nOutput:\n%s\n\n' ...
                            'You can also do it manually by running from the H2sim root:\n' ...
                            '  git submodule update --init --recursive -- %s'], ...
              name, cmd, cmdout, relpath);
    end

    fprintf('Submodule ''%s'' downloaded (pinned commit, matching the tested H2sim version).\n\n', name);

end

