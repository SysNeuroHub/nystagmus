function pursuit2D(subject,varargin)
% pursuit task in two dimensions
%
% Call example 
% >>pursuit2D('TST','stimType',2,'nRepPerCond',1)
% Will run subject named TST
% 
% After paradigm starts:
% - Press Enter to see your eye - adjust focus of camera
% - Press 'a' for automatic luminance adjustment, then `c` to calibrate. 
% - Press 'space' when you first fixate and 'Enter' when
%   you've finished looking at the 9 spots. 
% - Press 'v' to validate your calibration. (Again, space to start, Enter at
% end)
% - Press Esc (twice) when you're ready to start protocol. 
% - Press F8 during task to recalibrate eyetracker
% - Press 'a' or 'z' to start a new trial
%
%   *********** Press "Esc" twice to exit a running experiment ************
%
%
% Trial structure:
% - acquire and hold fixation for a fixed (+ random) duration
% - pursuit of moving stimulus
% - keypress to start next trial
% Trials that are not completed (due to eye-tracking errors) get randomly
% shuffled back into sequence. 
%
% Implementation
% - initial static fixation spot is one "stimulus"
% - moving spot is drawn by a second "stimulus", with position controlled by
% a third stimulus function - targXY.m - whose only role is to update the
% XY position. 
% - movement starts at pre-defined time after fixation is acquired
%
% STIM-TYPE SELECTION
% 1 - sinusoids with 3 frequencies (horizontal, vertical, diagonal and circular)
% 2 - sum-of-sinusoids - 15 unique conditions in 2D
%
% Also 
% - restart Matlab before you collect data in a new session. 
% - close all potentially demanding programs like Firefox, Slack, GitKraken etc. 
%       File manager and terminals are fine.
% - turn off EyeLink monitor or push it to edge of table (after you know it's working)
% - close door to hallway
% - have lights in anteroom on, in testing room off
%
% Auto backup recent data:
% - Ctrl-Alt-T (pop a terminal)
% - syncit (mount syncitium) It may ask "password for nic", which is the computer password - (Hint: monkey)
%   it will also ask "Enter syncitium username", which is your standard Monash username and then password.
% - rsyncData (push data)
%
% Check frame drops:
% - fd=get(c.prms.frameDrop,'struct',true);
%
% 
%

%% PARAMETER DEFINITIONS

if ~exist('subject','var')
  error('No subject name provided. Type ''help facecal'' for usage information.');
end

validateattributes(subject,{'char'},{'nonempty'},'','subject',1);

% parse arguments...
p = inputParser();
p.KeepUnmatched = true;
p.addRequired('subject',@(x) validateattributes(x,{'char'},{'nonempty'}));
p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('stimType',1,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'})); % 0=step-ramp; 1=single sinusoids; 2=sum-of-sinusoids
p.addParameter('tDur',8000,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));  % trial duration (ms)
p.addParameter('nRepPerCond',3,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));  % number of repeats of each condition
p.addParameter('tolerance',6,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));  % (deg) eye tolerance - radius

%for gabor patch
p.addParameter('contrast',1);
p.addParameter('orientation',0);
p.addParameter('frequency',1);
p.addParameter('spd',1);
p.addParameter('sigma',3);

p.addParameter('flickerFrequency',1);%not necessary?
p.addParameter('color',[0 0 .5]); 
p.addParameter('colorPolarity',[1 1 1]); 

p.parse(subject,varargin{:});

args = p.Results;


import neurostim.*
commandwindow;

%% ========= Specify rig configuration  =========

%Create a Command and Intelligence Centre object (the central controller for everything). Here a cic is returned with some default settings for this computer, if it is recognized.
if ~args.debug
    c = marmolab.rigcfg('debug',args.debug, p.Unmatched); % set to false to save githash at start of each experiment!
else
    c = dsOffice('smallWindow',true);
end
c.paradigm = 'pursuit2D';

% if ~args.debug % log git hash
%   hash = marmolab.getGitHash(fileparts(mfilename('fullpath')));
%   c.githash('pursuit.git') = hash;
% end

%Make sure there is an eye tracker (or at least a virtual one)
if isempty(c.pluginsByClass('eyetracker')) 
    e = neurostim.plugins.eyetracker(c);      %Eye tracker plugin not yet added, so use the virtual one. Mouse is used to control gaze position (click)
    e.useMouse = true;   
 end


%% ============== Add stimuli ==================

%Fixation dot
f = stimuli.fixation(c,'fix');    % Add a fixation stimulus object (named "fix") to the cic. It is born with default values for all parameters.
f.shape = 'CIRC';               %The seemingly local variable "f" is actually a handle to the stimulus in CIC, so can alter the internal stimulus by modifying "f".               
f.size = 0.25; %0.25; % units? 
f.color = [1 1 1];
f.on=0;                         % What time should the stimulus come on? (all times are in ms)
f.duration = '@f1.startTime.fixating+600'; % Show spot briefly after fixation acquired
f.X = '@traj.X';
f.Y = '@traj.Y';

% draw moving fixation target
tDur = args.tDur; % (ms) applies to fixation target and trajectory
% fm = stimuli.fixation(c,'targ');
% fm.X = '@traj.X';
% fm.Y = '@traj.Y';
% fm.on = '@traj.startTime'; % was .on
% fm.duration = tDur;
% fm.shape = f.shape; %'CIRC';
% fm.size = f.size; %0.25; %0.25; % units? 
% fm.color = f.color; %[1 1 1];

% draw gabor patch
fm = tuning.cgabor(c,'gabor'); % neurostim.stimuli.gabor
fm.color = args.color;
fm.colorPolarity = args.colorPolarity;
fm.width = 4*max(args.sigma);
fm.height = fm.width;
fm.X = '@traj.X';
fm.Y = '@traj.Y';
fm.on = '@traj.startTime'; % was .on
fm.duration = tDur;
%  tMoveWin1 = args.tPreBlank + args.tPreStat;
%     fm.addProperty('tMoveWin1',tMoveWin1);
%     tMoveWin2 = tMoveWin1 + args.tMove;
%     fm.addProperty('tMoveWin2',tMoveWin2);    
%         fm.phaseSpeed = '@iff(cic.trialTime>=gabor.tMoveWin1 & cic.trialTime<gabor.tMoveWin2, gabor.spd,0)';
fm.frequency = args.frequency;
fm.orientation = args.orientation;
fm.contrast = args.contrast;
fm.spd = args.spd;
fm.flickerMode = 'none';
fm.flickerFrequency = args.flickerFrequency;
fm.square = true;


%% trajectory definition
% there are no internal drawing commands - this stimulus just updates position .X/Y
t = targXY(c,'traj');
t.on = '@fix.stopTime'; % was .off';
t.duration = tDur; % ms - pursuit duration


% define properties for targXY trajectory control
% see also experimental design - factorial setup further below
c.traj.addProperty('stimType',args.stimType);
switch c.traj.stimType
    case 0 
        %% STEPRAMP
        t.sigFun = @(x) targXY.stepRamp(x);
        t.sigParams = [];
        c.traj.addProperty('dir',15);
        c.traj.addProperty('sp',6);
        c.traj.addProperty('step', -1);
        t.initPos = [5 -5];
        t.offset = [5 -5];
    
    case 1 
        %% SINGLE SINUSOIDS (H/V/diag/circle)
        t.sigFun = @(x) targXY.sumSines(x);
        % t.sigParams = [];
        
        c.traj.addProperty('fund',0.25); % (Hz) fundamental frequency
        c.traj.addProperty('freqMult',1); % component multipliers for each sinusoid
        c.traj.addProperty('xAmp',5); % separate X/Y component amplitudes and phases
        c.traj.addProperty('xPh',0);
        c.traj.addProperty('yAmp',0);
        c.traj.addProperty('yPh',0);
        
        t.initPos = []; % leave empty to start at first frame of trajectory
        t.offset = [0 0];
        
    case {2,3} % 2D sum-of-sinusoids 
        %% SUM SINES
        t.sigFun = @(x) targXY.sumSines(x);
        
        c.traj.addProperty('fund',0.25); % (Hz) fundamental frequency
        c.traj.addProperty('freqMult',[1 2]); % component multipliers for each sinusoid
        c.traj.addProperty('xAmp',[2.5 2.5]); % separate X/Y component amplitudes and phases
        c.traj.addProperty('xPh',[0 0]);
        c.traj.addProperty('yAmp',[2.5 2.5]);
        c.traj.addProperty('yPh',[3*pi/4 pi/4]);
        
        t.initPos = []; % leave empty to start at first frame of trajectory
        t.offset = [0 0];
end

%% NOISE FUNCTION
% t.noiseFun = @(x) targXY.gauss(x{:});
% t.nNoise = 1;
% t.noiseParams = {'sigma',0.5};


%% ========== Add required behaviours =========
%Subject's 2AFC response
k = behaviors.keyResponse(c,'choice');
k.from = '@traj.off'; % end of pursuit
k.maximumRT= Inf;                   %Allow inf time for a response
k.keys = {'a','z'};
k.required = false; %   setting false means that even if this behavior is not successful (i.e. the wrong answer is given), the trial will not be repeated.

%Maintain gaze on the fixation point until the dots disappear
g = behaviors.fixate(c,'f1');
g.from = 5000; % If fixation has not started at this time, move to the next trial
g.to = '@traj.off'; % stop tracking when trajectory ends
g.X = '@traj.X';
g.Y = '@traj.Y';
g.tolerance = args.tolerance; % (deg) allowed eye position error - should be aiming to get this as small as possible
g.required = true; % This is a required behavior. Any trial in which fixation is not maintained throughout will be retried. (See myDesign.retry below)
g.failEndsTrial = true; 


%% Turn off logging
stopLog(c.fix.prms.X);
stopLog(c.fix.prms.Y);
stopLog(c.targ.prms.X);
stopLog(c.targ.prms.Y);
stopLog(c.traj.prms.X);
stopLog(c.traj.prms.Y);
stopLog(c.f1.prms.X);
stopLog(c.f1.prms.Y);
c.traj.setChangesInTrial('X')
c.traj.setChangesInTrial('Y')


%% ========== Specify feedback/rewards ========= 
% Play a correct/incorrect sound for the 2AFC task
plugins.sound(c);           %Use the sound plugin

% Add correct/incorrect feedback
s= plugins.soundFeedback(c,'soundFeedback');
s.add('waveform','correct.wav','when','afterTrial','criterion','@choice.correct');
s.add('waveform','incorrect.wav','when','afterTrial','criterion','@ ~choice.correct');

%% Experimental design
c.trialDuration = Inf; %'@choice.stopTime';       %End the trial as soon as the 2AFC response is made.
k.failEndsTrial = true; 
k.successEndsTrial  = true; %false;

%  Specify experimental conditions
% For threshold estimation, we'd just vary speed
myDesign=design('myFac');                      %Type "help neurostim/design" for more options.
switch args.stimType
    case 0 
        %% STEPRAMP
    
    case 1 
        %% SINGLE SINUSOIDS (H/V/diag/circle)
        % For each component position is A.sin(2.pi.f.t)
        % Therefore peak speed is 2.pi.f.A. 
        % So to match peak speed as we vary frequency, 
        % we need to scale amplitude by 1/f
        xAmp = 6*[1 0 1 1]; % (deg)
        yAmp = 6*[0 1 1 1];
        freqMult = [1 2 3]; % note that fundamental frequency is defined above
        fundScale = 1./freqMult; %(fund/fund(1));
        
        myDesign.fac1.traj.xAmp = repmat(xAmp,[1 length(freqMult)]) .* repelem(fundScale, length(xAmp)); % scale amplitude by inverse of frequency multipler
        myDesign.fac1.traj.yAmp = repmat(yAmp,[1 length(freqMult)]) .* repelem(fundScale, length(yAmp));
        myDesign.fac1.traj.xPh = repmat([0 0 0 pi/2], [1 length(freqMult)]); % this will give H/V/diag/circle
        myDesign.fac1.traj.freqMult = repelem(freqMult,length(xAmp)); % frequency multiplier
        
%         % Without peak-speed matching
%         myDesign.fac1.traj.xAmp = 5*[1 0 1 1];
%         myDesign.fac1.traj.yAmp = 5*[0 1 1 1];
%         myDesign.fac1.traj.xPh = [0 0 0 pi/2];
%         myDesign.fac2.traj.fund = [0.25 0.5 0.75]; % frequency (Hz)
        
    case 2 
        %% 2D SUM SINES
        % We have 15 unique conditions. Simplest approach seems to be to
        % define them into a big cell array
        amp = 6;
        % define 3 circles with different frequencies
        for ix = 1:3
            freqMult{ix} = ix; %#ok<*AGROW>
            xAmp{ix} = amp/ix;    yAmp{ix} = amp/ix;
            xPh{ix} = pi/2;        yPh{ix} = 0;
        end
        
        % DEFINE 4 H2V3/H3V2 conditions
        ix = 4;
        freqMult{ix} = [2 3];
        xAmp{ix} = [amp/freqMult{ix}(1) 0]; yAmp{ix} = [0 amp/freqMult{ix}(2)];
        xPh{ix} = [0 0]; yPh{ix} = [0 0];

        % copy and slightly tweak next conditions
        ix = 5;
        freqMult{ix} = [3 2];
        xAmp{ix} = [amp/freqMult{ix}(1) 0]; yAmp{ix} = [0 amp/freqMult{ix}(2)];
        xPh{ix} = [0 0]; yPh{ix} = [0 0];
        
        ix = 6;
        freqMult{ix} = [2 3];
        xAmp{ix} = [amp/freqMult{ix}(1) 0]; yAmp{ix} = [0 amp/freqMult{ix}(2)];
        xPh{ix} = [0 0]; yPh{ix} = [0 pi/3];

        % copy and slightly tweak next conditions
        ix = 7;
        freqMult{ix} = [3 2];
        xAmp{ix} = [amp/freqMult{ix}(1) 0]; yAmp{ix} = [0 amp/freqMult{ix}(2)];
        xPh{ix} = [0 pi/3]; yPh{ix} = [0 0];
        
        % DEFINE 4 H1H2V1V2 conditions - they only vary in phase
        for ix = 8:11
            freqMult{ix} = [1 2];
            xAmp{ix} = [amp amp]/2; yAmp{ix} = [amp amp]/2;
        end
        ix = 8;  xPh{ix} = [pi/2 0]; yPh{ix} = [0 0];
        ix = 9;  xPh{ix} = [pi/2 0]; yPh{ix} = [0 pi/2];
        ix = 10; xPh{ix} = [pi/5 pi/3]; yPh{ix} = [pi/2 pi/7];
        ix = 11; xPh{ix} = [3*pi/4 2*pi/3]; yPh{ix} = [pi/2 pi/7];
        
        % DEFINE 4 speed-balanced H2H3V2V3 conditions
        for ix = 12:15
            freqMult{ix} = [2 3];
            xAmp{ix} = amp./freqMult{ix}; yAmp{ix} = amp./freqMult{ix};
        end
        ix = 12; xPh{ix} = [pi/2 0]; yPh{ix} = [0 0];
        ix = 13; xPh{ix} = [pi/2 0]; yPh{ix} = [0 pi/2];
        ix = 14; xPh{ix} = [pi/2 0]; yPh{ix} = [0 2*pi/3];
        ix = 15; xPh{ix} = [3*pi/4 2*pi/3]; yPh{ix} = [pi/2 pi/7];
        
        myDesign.fac1.traj.freqMult = freqMult;
        myDesign.fac1.traj.xAmp = xAmp;
        myDesign.fac1.traj.yAmp = yAmp;
        myDesign.fac1.traj.xPh = xPh;
        myDesign.fac1.traj.yPh = yPh;
        
    case 3
        %% SINGLE SINUSOIDS * 3 frequencies (H/V/circle)
        % For each component position is A.sin(2.pi.f.t)
        % Therefore peak speed is 2.pi.f.A. 
        % So to match peak speed as we vary frequency, 
        % we need to scale amplitude by 1/f

        amp = 6;
        % H sinusoids
        for ix = 1:3
            freqMult{ix} = ix; %#ok<*AGROW>
            xAmp{ix} = amp/ix;    yAmp{ix} = 0;
            xPh{ix} = 0;        yPh{ix} = 0;
        end
        
        % V sinusoid
        for a = 1:3
            ix = a+3;
            freqMult{ix} = a; %#ok<*AGROW>
            xAmp{ix} = 0;    yAmp{ix} = amp/a;
            xPh{ix} = 0;        yPh{ix} = 0;
        end        
        
        % define 3 circles with different frequencies
        for a = 1:3
            ix = a+6;
            freqMult{ix} = a; %#ok<*AGROW>
            xAmp{ix} = amp/a;    yAmp{ix} = amp/a;
            xPh{ix} = pi/2;        yPh{ix} = 0;
        end

        % DEFINE 2 H2V3/H3V2 conditions
        ix = 10;
        freqMult{ix} = [2 3];
        xAmp{ix} = [amp/freqMult{ix}(1) 0]; yAmp{ix} = [0 amp/freqMult{ix}(2)];
        xPh{ix} = [0 0]; yPh{ix} = [0 0];

        ix = 11;
        freqMult{ix} = [3 2];
        xAmp{ix} = [amp/freqMult{ix}(1) 0]; yAmp{ix} = [0 amp/freqMult{ix}(2)];
        xPh{ix} = [0 pi/3]; yPh{ix} = [0 0];
        
        % DEFINE 4 H1H2V1V2 conditions - they only vary in phase
        for ix = 12:15
            freqMult{ix} = [1 2];
            xAmp{ix} = [amp amp]/2; yAmp{ix} = [amp amp]/2;
        end
        ix = 12;  xPh{ix} = [pi/2 0]; yPh{ix} = [0 0];
        ix = 13;  xPh{ix} = [pi/2 0]; yPh{ix} = [0 pi/2];
        ix = 14; xPh{ix} = [pi/5 pi/3]; yPh{ix} = [pi/2 pi/7];
        ix = 15; xPh{ix} = [3*pi/4 2*pi/3]; yPh{ix} = [pi/2 pi/7];
        
        % DEFINE 3 speed-balanced H2H3V2V3 conditions
        for ix = 16:18
            freqMult{ix} = [2 3];
            xAmp{ix} = amp./freqMult{ix}; yAmp{ix} = amp./freqMult{ix};
        end
        ix = 16; xPh{ix} = [pi/2 0]; yPh{ix} = [0 0];
        ix = 17; xPh{ix} = [pi/2 0]; yPh{ix} = [0 pi/2];
        ix = 18; xPh{ix} = [3*pi/4 2*pi/3]; yPh{ix} = [pi/2 pi/7];
        
        myDesign.fac1.traj.freqMult = freqMult;
        myDesign.fac1.traj.xAmp = xAmp;
        myDesign.fac1.traj.yAmp = yAmp;
        myDesign.fac1.traj.xPh = xPh;
        myDesign.fac1.traj.yPh = yPh;
        
end

myDesign.retry = 'RANDOM'; %'IMMEDIATE' or 'IGNORE';
myDesign.maxRetry = 10;  % Each condition will be retried up to this many times. 
 
a=1;
myBlk{a} = block('myBlock',myDesign);
myBlk{a}.nrRepeats = args.nRepPerCond; %params.nRepPerCond; %nRepeatsPerBlock;

if ~args.debug
    c.eye.doTrackerSetupEachBlock = true;
else
    c.eye.doTrackerSetupEachBlock = false;
    c.cursor = 'arrow';
end

%% Run the experiment.
c.subject = args.subject; %params.subj; %'NP';
c.run(myBlk{1})


% Quick check of framedrops
%
% load last file
% fd=get(c.prms.frameDrop,'struct',true);
% -OR-
% d = marmodata.mdbase('file',f,'loadArgs',{'loadEye',true});
% for a = 1:max(d.meta.cic.trial.data)
%   frDrop{a} = d.meta.cic.frameDrop('trial',a).data; 
% end
% nDrop = cellfun('length',frDrop)
