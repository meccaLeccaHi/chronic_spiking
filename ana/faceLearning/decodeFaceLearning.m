function [type, face, step] = decodeFaceLearning(stim,h)
%
% [type face step] = decodeFaceLearning(stim,h);
%
% Decoder function for faceLearning, passive screening stimuli.
%
% Calling the function with no input argument returns the lookup table:
% lut = decodeFaceLearning();
%
%look up table is organized as follows:
% LUT = [STIMID TYPE FACE STEP]
%
% UPDATED 2013-JUN-17
% Returns the appropriate lookup table for the date of data collection.
% PHASE I : 2013-may-21 to jun-09
% PHASE II : 2013-jun-10 to 2013-jul-14
% Phase III : 2013-jul-19 to 2014-jun-28
%
% 2013-jul-16 Fixed bug that was preventing return of lut2
% (Needed for post 2013-jun-10 training phase)
% 2013-jul-26 Added phase 3 condition in ph switch - apj
%
% last modified 2013-aug-10
% apj


if strcmp(h.expt,'stimscreen')
        if (round(datenum(h.date(1,:))) < h.training.training1Date)|| ...
                 (round(datenum(h.date(1,:))) > datenum('01-Jan-2014')) % length(h.stim.all)==206; 
            ph = 1;
        elseif round(datenum(h.date(1,:))) >= h.training.training1Date
            ph = 2;
        end
elseif strcmp(h.expt,'faceTraining')
    if round(datenum(h.date(1,:))) < h.training.training2Date
        ph = 3;
    elseif round(datenum(h.date(1,:))) >= h.training.training2Date
        ph = 4;
    end
end


HUMAN = h.vv.HUMAN;
MONKEY = h.vv.MONKEY;
if isfield(h.vv,'REFERENCE')
    REFERENCE = h.vv.REFERENCE;
end
STEPS = h.stim.morphSteps;

if ph==1
    lut = nan(206,4);
    lut(:,1) = [1:206];
    lut(1:96,2) = HUMAN;
    lut(97:192,2) = MONKEY;
    lut(193,2) = HUMAN;
    lut(194,2) = MONKEY;
    lut(195:206,2) = REFERENCE;
    FACE = [];
    for n=1:24
        if n<=12
            FACE = [FACE ; n*ones(8,1)];
        else
            FACE = [FACE ; 100+((n-12)*ones(8,1))];
        end
    end
    lut(1:192,3) = FACE;
    lut(193,3) = 0;
    lut(194,3) = 100;
    lut(195:206,3) = 201:212;
    nsteps = length(STEPS);
    for s=1:nsteps
        index = s:nsteps:192;
        lut(index,4) = STEPS(s);
    end
    lut([193 194],4) = 0;
    
elseif ph==2
    controlHumans = [1 5 7 12];
    controlMonkeys = [1 2 4 5];
    testMonkeys = [6 7 8 10];
    lut = nan(110,4);
    lut(:,1) = [1:110];
    lut(1:32,2) = HUMAN;
    lut(33:96,2) = MONKEY;
    lut(97,2) = HUMAN;
    lut(98,2) = MONKEY;
    lut(99:110,2) = REFERENCE;
    FACE = [];
    for n=1:12
        if ismember(n,[1:4]);
            FACE = [FACE ; controlHumans(n)*ones(8,1)];
        elseif ismember(n,[5:8])
            FACE = [FACE ; 100+controlMonkeys(n-4)*ones(8,1)];
        elseif ismember(n,[9:12])
            FACE = [FACE ; 100+testMonkeys(n-8)*ones(8,1)];
        end
    end
    lut(1:96,3) = FACE;
    lut(97,3) = 0;
    lut(98,3) = 100;
    lut(99:110,3) = 201:212;
    nsteps = length(STEPS);
    for s=1:nsteps
        index = s:nsteps:96;
        lut(index,4) = STEPS(s);
    end
    lut([97 98],4) = 0;
    
elseif ph==3
    controlMonkeys = [1 2 4 5];
    lut = nan(4,4);
    lut(:,1) = 1:4;
    lut(:,2) = MONKEY;
    lut(:,3) = controlMonkeys;
    lut(:,4) = ones(1,4).*100;
    
elseif ph == 4
    controlMonkeys = [1 2 4 5];
    lut = nan(401,4);
    lut(:,1) = 1:401;
    lut(:,2) = MONKEY;
    FACE = [];
    for n=1:4
        FACE = [FACE ; controlMonkeys(n)*ones(100,1)];
    end
    lut(1:end-1,3) = FACE;
    
    nsteps = length(STEPS);
    for s = 1:nsteps
        index = s:nsteps:400;
        lut(index,4) = STEPS(s);
    end
    lut(end,3) = 100;
    lut(end,4) = 0;
end

getrows = find(ismember(lut(:,1),stim));
type = lut(getrows,2);
face = lut(getrows,3);
step = lut(getrows,4);


