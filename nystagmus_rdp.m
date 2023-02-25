function nystagmus_rdp(subject,varargin)
% two eyes receiving the same direction, specified as 'ori1List'
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
%p.addParameter('stimType',1,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'})); % 0=step-ramp; 1=single sinusoids; 2=sum-of-sinusoids
p.addParameter('tDur',8000,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));  % trial duration (ms)
p.addParameter('nRepPerCond',3,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));  % number of repeats of each condition
p.addParameter('tolerance',6,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));  % (deg) eye tolerance - radius

%for gabor patch
p.addParameter('dir1List_r',[0]); %direction(s) of red dots [deg] 0: left to right, 90: bottom to top
p.addParameter('dir1List_b',[0]); %direction(s) of blue dots [deg]
p.addParameter('speed',4); %[deg]
p.addParameter('radius',5); %aperture size [pix]
p.addParameter('dotSize',5); %dot size [pix]
p.addParameter('nrDots',30); %number of dots
p.addParameter('coherence',1);
p.addParameter('tPreBlank',0);

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
f.X = 0;%'@traj.X';
f.Y = 0;%'@traj.Y';

% draw moving fixation target
tDur = args.tDur; % (ms) applies to fixation target and trajectory

%c.addProperty('tPreBlank',args.tPreBlank);

s = RandStream('mt19937ar');

nrConds = 2;
fm = cell(nrConds,1);
for ii = 1:nrConds
    
    reset(s, 1);%args.rngSeed);

    stimName = ['patch' num2str(ii)];
    %patch1: blue dots, patch2: red dots
    
    % draw rdp patch
    fm{ii} = neurostim.stimuli.rdp(c,stimName);
    
    %common parameters across stim
    fm{ii}.X = 0;
    fm{ii}.Y = 0;

    
    %rdp specific parameters
    fm{ii}.size = args.dotSize; %dot size [px]
    fm{ii}.type = 0; %square dot
    fm{ii}.maxRadius =  args.radius;%maximum radius of aperture (px)
    fm{ii}.speed = args.speed; %dot speed (deg
    fm{ii}.direction = 0;
    fm{ii}.nrDots = args.nrDots;
    fm{ii}.coherence = 0.7; %dot coherence [0-1]
    fm{ii}.motionMode = 1; %linear
    fm{ii}.lifetime = 50;%lifetime of dots (in frames)
    fm{ii}.dwellTime = 1;
    fm{ii}.coordSystem = 0; %polar coordinates
    fm{ii}.noiseMode = 0; %proportion
    fm{ii}.noiseDist = 1; %uniform 
    fm{ii}.rngSeed = 1;
    %fm{ii}.noiseWidth not used when noiseDist=1
    %fm{ii}.truncateGauss ??
    
end
fm{1}.color = [0 0 1 .5]; %blue
fm{2}.color = [1 0 0 .5]; %red
    
%% ========== Add required behaviours =========
%Subject's 2AFC response
k = behaviors.keyResponse(c,'choice');
k.from = '@patch1.off'; % end of pursuit
k.maximumRT= Inf;                   %Allow inf time for a response
k.keys = {'z'};
k.required = false; %   setting false means that even if this behavior is not successful (i.e. the wrong answer is given), the trial will not be repeated.

%Maintain gaze on the fixation point until the dots disappear
g = behaviors.fixate(c,'f1');
g.from = 5000; % If fixation has not started at this time, move to the next trial
g.to = '@patch1.off'; % stop tracking when trajectory ends
g.X = '@patch1.X';
g.Y = '@patch1.Y';
g.tolerance = args.tolerance; % (deg) allowed eye position error - should be aiming to get this as small as possible
g.required = true; % This is a required behavior. Any trial in which fixation is not maintained throughout will be retried. (See myDesign.retry below)
g.failEndsTrial = true; 


%% Turn off logging
stopLog(c.fix.prms.X);
stopLog(c.fix.prms.Y);
stopLog(c.patch1.prms.X);
stopLog(c.patch1.prms.Y);
stopLog(c.patch2.prms.X);
stopLog(c.patch2.prms.Y);
stopLog(c.f1.prms.X);
stopLog(c.f1.prms.Y);
%c.traj.setChangesInTrial('X');%?
%c.traj.setChangesInTrial('Y');%?


%% ========== Specify feedback/rewards ========= 
% Play a correct/incorrect sound for the 2AFC task
plugins.sound(c);           %Use the sound plugin

% Add correct/incorrect feedback
s= plugins.soundFeedback(c,'soundFeedback');
s.add('waveform','ding.wav','when','afterTrial','criterion','@ choice.correct');

%% Experimental design
c.trialDuration = Inf; %'@choice.stopTime';       %End the trial as soon as the 2AFC response is made.
k.failEndsTrial = true; 
k.successEndsTrial  = true; %false;

%  Specify experimental conditions
% For threshold estimation, we'd just vary speed
myDesign=design('myFac');                      %Type "help neurostim/design" for more options.


facOutList = {'direction'}; % frequency = spatial frequency
facInList = {'dir1List_b'};
for a = 1:length(facInList)
    myDesign.(sprintf('fac%d',a)).patch1.(facOutList{a}) = args.(facInList{a});
end
facInList = {'dir1List_r'};
for a = 1:length(facInList)
    myDesign.(sprintf('fac%d',a)).patch2.(facOutList{a}) = args.(facInList{a});
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
    c.eye.continuous = true;
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
