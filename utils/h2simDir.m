function dir = h2simDir()
    dir = fileparts(mfilename('fullpath'));
    dir = fullfile(dir, '..');
end
