function [paths, monkeys] = setFilePath(proj)
%
% Script for setting up appropriate paths, depending on the project specified
% and the machine being used.
% Usage: [paths, monkeys] = setFilePath(proj)
% PROJ is the project identifier number
% PATHS is the output structure fields
%
% Example output:
% paths.data       /data/mcmahond/faceLearning
% paths.raw        /rawdata/mcmahond/
% paths.rare       /data/mcmahond/faceLearning/rare
% paths.mas        /data/mcmahond/faceLearning/mas
% paths.ofs        ~/ana/OfflineSDK
% paths.han        ~/ana/han06
% paths.ana        ~/ana/faceLearning
% paths.results    ~/results/faceLearning
% paths.stim       ~/sim/faceLearning
%
% last modified 5-10-16
% apj

% get hostname
hostname = [];
while isempty(hostname)
[~, hostname]                       = system('hostname');
if isunix
    if regexp(hostname, 'adam-P55A-UD3')
        proc_drive                  = '/home/adam/Cloud3';
        paths.raid              = '/mnt/raid1';
        paths.scratch           = '/mnt/scratch/';
    elseif regexp(hostname, 'lab-All')

        %%  LAB MACHINE
        proc_drive                  = '/home/lab/Cloud3/';
        paths.raid              = [];
        paths.scratch           = '/mnt/scratch/';
    end
else
    if regexp(hostname, 'DESKTOP-GTV1LPS')
        keyboard
        %% need to fix these values to be correct- LAPTOP
        proc_drive                  = 'E:\';
        paths.raid              = [];
        paths.scratch           = [];
    elseif regexp(hostname, 'Adam-PC')
        proc_drive                  = '\\srvr1\LeopoldShare\';
        paths.raid              = '\\srvr1\Raid1\';
        paths.scratch           = '\\srvr1\Scratch\';
    end
end
end

%% project
if proj==1
    proj_title                      = 'categorySelectivity';
    proj_code                       = 'cs'; % ?
elseif proj==2
    proj_title                      = 'faceLearning';
    proj_code                       = 'fl'; % ?
elseif proj==3
    proj_title                      = 'radTan';
    proj_code                       = 'rt'; % ?
end

% monkey names
monkeys{1}                          = 'toroid';
monkeyCodes{1}                      = 'tor';
monkeys{2}                          = 'rhombus';
monkeyCodes{2}                      = 'rom';
monkeys{3}                          = 'matcha';
monkeyCodes{3}                      = 'mat';
monkeys{4}                          = 'spice';
monkeyCodes{4}                      = 'spc';

paths.proj                     = proj_title;
% is this necessary?

paths.code                          = proj_code;
paths.drive                         = proc_drive;

if ~isempty(paths.raid)
    paths.raw          = fullfile(paths.raid, 'raw_data');
end

paths.data         = fullfile(proc_drive, 'proc_data', proj_title);
paths.rare         = fullfile(paths.data, 'rare');
paths.sortFigs     = fullfile(paths.data, 'spikeSort_figs');
paths.mas          = fullfile(paths.data, 'mas');
paths.projects     = fullfile(proc_drive, 'projects');
paths.physio       = fullfile(paths.projects, 'Physio');
paths.ana          = fullfile(paths.projects, 'ana');
paths.anaProj      = fullfile(paths.projects, 'ana', proj_title);
paths.analysis     = fullfile(paths.projects, 'analysis', proj_title);
paths.preProc      = fullfile(paths.ana, 'preProc');
paths.han          = fullfile(paths.ana, 'han06');
paths.ofs          = fullfile(paths.preProc, 'OfflineSDK');
paths.wave         = fullfile(paths.preProc, 'wave_clus_2.0wb', 'Wave_clus');
paths.ref          = fullfile(paths.projects, 'reference_code');
paths.results      = fullfile(paths.projects, 'results', proj_title);
paths.stim         = fullfile(paths.projects, 'stimuli', proj_title);
paths.stimCode     = fullfile(paths.projects, 'stimuli', 'stimCode', proj_title);

addpath(paths.preProc);
addpath(paths.anaProj);
% addpath(genpath(paths.wave));
addpath(paths.analysis);
addpath(paths.ofs);


