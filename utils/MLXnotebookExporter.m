classdef MLXnotebookExporter
    methods


        function setupMfromMlx(mne, filename, varargin)
        % Setup M file from mlx

            opt = struct('outputDirectory', []);
            opt = merge_options(opt, varargin{:});

            assert(exist(filename, 'file') == 2, 'File %s not found.', filename); 
            inputfilename = which(filename);
            
            if isempty(opt.outputDirectory)
                % if the file is in a directory denoted notebooks, we move the output in the directory above
                [dirpath, filename, ext] = fileparts(inputfilename);
                dirpaths = split(dirpath, filesep);
                lastdir = dirpaths{end};
                if strcmp(lastdir, 'notebooks')
                    outputDirectory = strrep(dirpath, [filesep, 'notebooks'], '');
                elseif strcmp(ext, '.m')
                    inputfilename = fullfile(dirpath, 'notebooks', [filename, '.mlx']);
                    outputDirectory = dirpath;
                else
                    outputDirectory = dirpath;
                end
            else
                outputDirectory = opt.outputDirectory;
            end

            outputfilename = fullfile(outputDirectory, [filename, '.m']);

            export(inputfilename, outputfilename, 'format', 'm');            
            
        end
        
        
        function setupIpynbFromMlx(mne, filename, varargin)

            opt = struct('outputDirectory', [], ...
                         'removeSolverOutput', true);
            opt = merge_options(opt, varargin{:});
            
            assert(exist(filename, 'file') == 2, 'File %s not found.', filename);
            
            fullfilename = which(filename);
            [dirpath, filename, ext] = fileparts(fullfilename);
                                   
            if strcmp(ext, '.mlx')
                % same directory
                inputfilename   = fullfilename;
                outputDirectory = dirpath;
            else
                outputDirectory = fullfile(dirpath, 'notebooks');
                inputfilename = fullfile(outputDirectory, [filename, '.mlx']);
            end            

            if ~isempty(opt.outputDirectory)
                outputDirectory = opt.outputDirectory;
            end

            outputfilename = fullfile(outputDirectory, [filename, '.ipynb']);
            
            export(inputfilename, outputfilename, 'format', 'ipynb');

            if opt.removeSolverOutput
                % Remove the solver output
                fid = fopen(outputfilename, 'r+');
                txt = fread(fid, '*char')';
                fclose(fid);

                txt = mne.cleanup(txt);

                fid = fopen(outputfilename, 'w+');
                fwrite(fid, txt);
                fclose(fid);
            end

            pyfilename = fullfile(battmoDir(), 'Utilities', 'Various', 'setupIpynbForH2sim.py');

            pyrunfile([pyfilename ' ' outputfilename]);
            
        end

        
        function setupMlxFromM(mne, filename, varargin)

            opt = struct('run'            , false, ...
                         'outputDirectory', []);
            opt = merge_options(opt, varargin{:});

            [inputfile, outputfile] = getIOfiles(mne, filename, ...
                                                 'outputDirectory', opt.outputDirectory, ...
                                                 'outputFormat', 'mlx');
            
            matlab.internal.liveeditor.openAndSave(inputfile, outputfile);

            if opt.run
                matlab.internal.liveeditor.executeAndSave(outputfile);
            end
            
        end
        
        function setupIpynbFromM(mne, filename, varargin)

        % From a M-file, generate the notebooks (mlx and ipynb)
        % The default directory for the generated notebooks is a sub-directory where the M-file is located called 'notebooks'

            opt = struct('run'            , false, ...
                         'outputDirectory', []   , ...
                         'generateIpynb'  , true);
            opt = merge_options(opt, varargin{:});

            mne.setupMlxFromM(filename, ...
                              'run'            , opt.run, ...
                              'outputDirectory', opt.outputDirectory);

            if opt.generateIpynb
                mne.setupIpynbFromMlx(filename, 'outputDirectory', opt.outputDirectory);
            end
            
        end

        function [inputfile, outputfile] = getIOfiles(mne, filename, varargin)
            
            opt = struct('outputDirectory', [], ...
                         'outputFormat', 'mlx');
            
            opt = merge_options(opt, varargin{:});

            assert(exist(filename, 'file') == 2, 'File %s not found.', filename); 
            fullfilename = which(filename);
            [dirpath, filename, ext] = fileparts(fullfilename);

            if isempty(opt.outputDirectory)
                outputDirectory = fullfile(dirpath, 'notebooks');
                if ~exist(outputDirectory)
                    % Create the directory if it does not exist
                    mkdir(outputDirectory);
                end
            else
                outputDirectory = opt.outputDirectory
            end

            switch opt.outputFormat
              case 'mlx'
                inputfile  = fullfile(dirpath, [filename, '.m']);
                outputfile = fullfile(outputDirectory, [filename, '.mlx']);
              case 'ipynb'
                inputfile  = fullfile(outputDirectory, [filename, '.mlx']);
                outputfile = fullfile(outputDirectory, [filename, '.ipynb']);
              otherwise
                error('outputFormat not recognized');
            end
            
        end

        function runMlxAndSave(mne, filename)

        % To run and update the mlx notebook programmatically, it is possible to use:
        % matlab.internal.liveeditor.executeAndSave('fullpathnameto.mlx')
            
            matlab.internal.liveeditor.executeAndSave(filename);

        end

    end

    methods (Static)

        function txt = cleanup(txt)

            [txt, found_one] = MLXnotebookExporter.cleanup_one(txt);

            while found_one
                [txt, found_one] = MLXnotebookExporter.cleanup_one(txt);
            end

        end
        
        function [txt, found_one] = cleanup_one(txt)
            

            pos = strfind(txt, "Solving timestep");

            if isempty(pos)
                found_one = false;
                return
            else
                found_one = true;
            end
            
            pos = pos(1);

            %  find the first occurence of 'output' before 'Solving timestep'
            outputpos = strfind(txt, 'output');
            outputpos = outputpos(outputpos < pos);
            outputpos = outputpos(end);

            %% We now find the position of the matching square brackets after the occurence of 'output'
            
            bleft  = strfind(txt, '[');
            bright = strfind(txt, ']');

            bleft  = [bleft', ones(numel(bleft), 1)];
            bright = [bright', -ones(numel(bright), 1)];

            bracket = [bleft; bright];

            [~, ind] = sort(bracket(:, 1));

            bracket = bracket(ind, :);

            bracket(:, 2) = cumsum(bracket(:, 2));

            ind = find(bracket(:, 1) > outputpos, 1);

            bracket = bracket(ind : end, :);
            
            % Position after the opening square bracket
            startpos = bracket(1, 1) + 1;

            ind = find(bracket(:, 2) == bracket(1, 2) - 1, 1);

            % Position before the closing square bracket
            endpos = bracket(ind, 1) - 1;

            %% We remove the captured text
            txt = [txt(1 : startpos), ...
                   txt(endpos + 1: end)];

        end

    end

end

