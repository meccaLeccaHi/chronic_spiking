function header = initRadTanHeader(expt)
%
% Function for initializing a blank header structure for radTan
% project.
%
% last modified 6-17-16
% apj


% Specify variables
if isequal(expt,'stimscreen')
    vv.SESSION = 1;
    vv. DATE = 2;
    vv.TR = 3;
    vv.STIM = 4;
    vv.TYPE = 5; % HUMAN, MONKEY, or REFERENCE
    vv.FACE = 6; % which individual; 0:12 humans, 100:112 monkeys
    vv.STEP = 7; % percentage morph
    vv.DIRECTION = 8; % radial v. tangential morph direction
    vv.TREL = 9; % time relative to start of first trial
    vv.T_PRE = 10;
    vv.T_STIMON = 11;
    vv.T_STIMOFF = 12;
    vv.T_POST = 13;
    vv.SPIKE_DIMS = 80;
    vv.HUMAN = 1;
    vv.MONKEY = 2;
    vv.REFERENCE = 3;
    vv.RADIAL = 1;
    vv.TANGENTIAL = 2;
end

stim.set = 'radTan';
stim.all = [1:536];
stim.radSteps = [-100 -50 -25 -10 10 25 50 75 100 125 150 200];
stim.tanSteps = [10 25 50 75 90];
stim.img = 'imgRadTan.mat';

% raw data source files
src.tank = [];
src.block = [];
src.plx = []; % needed for plx files with multiple blocks.

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
header.src = src;
header.bands = [];
%header.lfp_channels = [];
%header.Taligned = [];


