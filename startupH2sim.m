function startupH2sim

% This startup file set up the MATLAB path
%
%% We use first `MRST <https://github.com/SINTEF-AppliedCompSci/MRST>`_ setup for MRST modules.
% The source code for MRST is synchronized to H2Sim using git-submodule mechanisms (In the MRST directory in BattMo, you
% should find the subdirectories given by the ``names`` cell array below)
%

fprintf('\n  /-\\\n | + | H2Sim\n  \\-/\n\n');
fprintf('Welcome to the H2 Storage Modeling Toolbox (H2sim)!\n');
fprintf('H2sim is based on MRST, which will now be initialized.\n\n');

rootdirname = fileparts(mfilename('fullpath'));

run(fullfile(rootdirname, 'submodules', 'mrst', 'startup'));

dirnames = {'h2-biochem', ...
            'h2-store', ...
            'utils'};

for ind = 1 : numel(dirnames)
    dirname = fullfile(rootdirname, dirnames{ind});
    addpath(genpath(dirname));
end

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

%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
