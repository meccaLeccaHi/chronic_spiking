function header = initFaceLearningHeader(expt)
%
% Function for initializing a blank header structure for faceLearning
% project.
%
% last modified 6-17-16
% apj

% Specify variables
if isequal(expt,'stimscreen')
    vv.SESSION = 1;
    vv.DATE = 2;
    vv.TR = 3;
    vv.STIM = 4;
    vv.TYPE = 5; % HUMAN, MONKEY, or REFERENCE
    vv.FACE = 6; % which individual; 0:12 humans, 100:112 monkeys
    vv.STEP = 7; % morph step, [0:100]
    vv.TREL = 8; % time relative to start of first trial
    vv.T_PRE = 9;
    vv.T_STIMON = 10;
    vv.T_STIMOFF = 11;
    vv.T_POST = 12;
    vv.SPIKE_DIMS = 80;
    vv.HUMAN = 1;
    vv.MONKEY = 2;
    vv.REFERENCE = 3;
elseif isequal(expt,'faceTraining')
    vv.SESSION = 1;
    vv.DATE = 2;
    vv.TR = 3;
    vv.STIM = 4;
    vv.TYPE = 5; % HUMAN, MONKEY, or REFERENCE
    vv.FACE = 6; % which individual; 0:12 humans, 100:112 monkeys
    vv.STEP = 7; % morph step, [0:100]
    vv.RESP = 8;
    vv.OUTCOME = 9;
    vv.TREL = 10; % time relative to start of first trial
    vv.T_PRE = 11;
    vv.T_STIMON = 12;
    vv.T_STIMOFF = 13;
    vv.T_SAC = 14;
    vv.T_REWARD = 15;
    vv.T_FIXBREAK = 16;
    vv.T_POST = 17;
    vv.SPIKE_DIMS = 200;
    vv.HUMAN = 1;
    vv.MONKEY = 2;
    vv.REFERENCE = 3;
    % %% OUTCOME codes
    %vv.FIXBREAK = 666;% fixbreak before T_STIMOFF
    %vv.NOSTIM = 0; % monkey never acquires fixation
    vv.CORRECT = 1;
    vv.WRONG = 0;
    vv.INVALID = -999;
end

% training history: First stimscreen file collected after learning.
training.startDate = 735375; %2013-may-21
training.training1Date = 735395; % 2013-jun-10; first day of saccade training
training.training2Date = 735434; % 2013-jul-19; first day on saccade tasks with adaptive threshold
training.training3Date = 735771; % update training3Date when possible - date switched to rhombus

training.filename = [];
training.session = [];

stim.set = 'faceLearning';
stim.all = [1:206];
stim.morphSteps = [1 2 3 6 12 25 50 100];
stim.img = 'imgFaceLearning.mat';
% stim.set = [];
% stim.all = [];
% stim.morphSteps = [];
% stim.img = [];

% raw data source files
src.tank = [];
src.block = [];
src.spike = []; % needed for spike files with multiple blocks.

% init header fields
header.filename = [];
header.expt = expt;
header.stim = stim;
header.session = [];
header.date = [];
header.ver = [];
header.parseDate = [];
header.vv = vv;
header.timestamps = [vv.T_PRE : vv.T_POST];
header.chan = [];
header.snames = [];
header.iso = [];
header.training = training;
header.src = src;
header.bands = [];
%header.lfp_channels = [];
%header.Taligned = [];


