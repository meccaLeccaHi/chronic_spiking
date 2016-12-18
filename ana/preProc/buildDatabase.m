function buildDatabase(header,paths,varargin)
%
% Usage: buildDatabase(header,paths,varargin)
% HEADER file-specific and project-specific information
% PATHS project path structure containing at minimum paths.rare field
% Batch script for parsing stimscreen database collected on TDT
% system.
%
% last modified 6-14-16
% apj

% step through each recording day in header
for n = 1:length(header)
    
    % name file
    filename            = [header{n}.filename header{n}.src.block header{n}.suffix.save '2'];
    savename            = fullfile(paths.mas,filename);
    
    if n==1||n==length(header)
        fprintf('%s\n',['Building: ' savename])
    else
        fprintf('%s\n',['Building: ' filename])
    end
    
    % add session to header
    header{n}.session   = n;
    
    % parse
    dat                 = parseData(header{n},paths,varargin{1});

    keyboard
    
    if ~isempty(dat.s)
        % save file
        save(savename,'dat');
        if n==1||n==length(header)
            fprintf('%s\n',['Saved: ' savename])
        else
            fprintf('%s\n',['Saved: ' filename])
        end
    else
        fprintf('%s\n',['NOT saved: ' filename ' ... .dat empty'])
    end
    
end