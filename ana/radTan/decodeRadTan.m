function [type, face, step, direction] = decodeRadTan(stim,h)
%
% [type face step] = decodeRadTan(stim,h);
%
% Decoder function for radTan, passive screening stimuli.
%
% Calling the function with no input argument returns the lookup table:
% lut = decodeRadTan();
%
%look up table is organized as follows:
% LUT = [STIMID TYPE FACE STEP DIRECTION]
%
% last modified 5-15-16
% apj

% radLevels = [-100 -50 -25 -10 10 25 50 75 100 125 150 200];
% tanLevels = [10 25 50 75 90];

HUMAN = h.vv.HUMAN;
MONKEY = h.vv.MONKEY;
RADIAL = h.vv.RADIAL;
TANGENTIAL = h.vv.TANGENTIAL;

if isfield(h.vv,'REFERENCE')
    REFERENCE = h.vv.REFERENCE;
end
RADSTEPS = h.stim.radSteps;
TANSTEPS = h.stim.tanSteps;

lut = nan(536,5);
lut(:,1) = [1:536];

% for type
lut(1:12*(length(RADSTEPS)*2+length(TANSTEPS)),2) = HUMAN;
lut((12*(length(RADSTEPS)*2+length(TANSTEPS)))+1:18*(length(RADSTEPS)*2+length(TANSTEPS)),2) = MONKEY;
lut(523:534,2) = REFERENCE;
lut(535,2) = HUMAN;
lut(536,2) = MONKEY;

% for face
for i = 1:12
    lut([1:12]+((i-1)*29),3) = i+(i-1);
    lut([13:24]+((i-1)*29),3) = i+1+(i-1);
    lut([25:29]+((i-1)*29),3) = i+.5+(i-1);
end

% for step
for i = 1:24
    temp = length(TANSTEPS)+2*length(RADSTEPS);
    lut([1:2*length(RADSTEPS)]+temp*(i-1),4) = [RADSTEPS RADSTEPS];
    lut([25:25+length(TANSTEPS)-1]+temp*(i-1),4) = TANSTEPS;
end

% for direction
for i = 1:24
    temp = length(TANSTEPS)+2*length(RADSTEPS);
    lut([1:2*length(RADSTEPS)]+temp*(i-1),5) = RADIAL;
    lut([25:25+length(TANSTEPS)-1]+temp*(i-1),5) = TANGENTIAL;
end

getrows = find(ismember(lut(:,1),stim));
type = lut(getrows,2);
face = lut(getrows,3);
step = lut(getrows,4);
direction = lut(getrows,5);




