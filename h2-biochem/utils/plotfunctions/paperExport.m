function files = paperExport(fig, name, varargin)
% Export a paper/presentation figure to EPS and PNG with one call.
%
% SYNOPSIS:
%   paperExport(fig, 'h2_loss_over_time')
%   files = paperExport(fig, 'h2_loss_over_time', 'dir', myDir)
%
% PARAMETERS:
%   fig  - Figure handle (e.g. the return value of paperFigure).
%   name - Base file name without extension. Non-alphanumeric characters are
%          replaced by underscores.
%
% OPTIONAL PARAMETERS:
%   'dir'        - Output directory. Default: a 'figures' sub-folder of the
%                  calling script's directory, or of pwd if that cannot be
%                  resolved. Created if missing.
%   'formats'    - Cell array of extensions to write. Default {'eps','png'}.
%   'resolution' - Raster resolution in DPI for PNG. Default 300.
%   'style'      - If true, run styleAxes on every axes of the figure before
%                  export. Default false.
%
% RETURNS:
%   files - Cell array of the full paths written.
%
% SEE ALSO:
%   paperFigure, styleAxes, paperColors

    opt = struct('dir', '', 'formats', {{'eps', 'png'}}, ...
        'resolution', 300, 'style', false);
    opt = merge_options(opt, varargin{:});

    assert(isgraphics(fig, 'figure'), 'fig must be a figure handle.');

    if isempty(opt.dir)
        opt.dir = defaultExportDir();
    end
    if ~isfolder(opt.dir)
        [ok, msg] = mkdir(opt.dir);
        assert(ok, 'Could not create export directory %s: %s', opt.dir, msg);
    end

    if opt.style
        axList = findall(fig, 'Type', 'axes');
        for k = 1:numel(axList)
            styleAxes(axList(k));
        end
    end

    base = regexprep(name, '[^A-Za-z0-9]+', '_');
    base = regexprep(base, '^_+|_+$', '');

    files = cell(1, numel(opt.formats));
    for k = 1:numel(opt.formats)
        ext = lower(strrep(opt.formats{k}, '.', ''));
        target = fullfile(opt.dir, [base, '.', ext]);
        writeOneFormat(fig, target, ext, opt.resolution);
        files{k} = target;
    end
    fprintf('Exported %s -> %s\n', base, opt.dir);
end

function writeOneFormat(fig, target, ext, resolution)
    switch ext
        case {'eps', 'epsc'}
            if useExportgraphics()
                exportgraphics(fig, target, 'ContentType', 'vector', ...
                    'BackgroundColor', 'white');
            else
                print(fig, target, '-depsc2', '-painters');
            end
        case 'pdf'
            if useExportgraphics()
                exportgraphics(fig, target, 'ContentType', 'vector', ...
                    'BackgroundColor', 'white');
            else
                print(fig, target, '-dpdf', '-painters');
            end
        case 'png'
            if useExportgraphics()
                exportgraphics(fig, target, 'Resolution', resolution, ...
                    'BackgroundColor', 'white');
            else
                print(fig, target, '-dpng', sprintf('-r%d', resolution));
            end
        otherwise
            error('paperExport:unsupportedFormat', ...
                'Unsupported export format "%s".', ext);
    end
end

function tf = useExportgraphics()
    tf = exist('exportgraphics', 'file') == 2 || ...
        exist('exportgraphics', 'builtin') == 5;
end

function outDir = defaultExportDir()
    stack = dbstack('-completenames');
    outDir = '';
    for k = 1:numel(stack)
        [folder, fname] = fileparts(stack(k).file);
        if ~strcmp(fname, mfilename) && ~isempty(folder)
            outDir = fullfile(folder, 'figures');
            return;
        end
    end
    outDir = fullfile(pwd, 'figures');
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
