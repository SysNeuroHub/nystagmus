function [fileOut, msg] = pursuit2D_Merge_nr(files, fileOut, varargin)
%
% load neurostim data and eye-movement data from multiple files using marmodb
%
% saves a combined data file providing simpler access to perceptual
% responses, stimulus trajectory and eye movement traces from complete
% trials
%
% If given a list of files, it will collate them (but do not
% similarity-checking, i.e. it just assumes the two files had identical
% stimulus properties
%
% Note that optional input (...,'align',true') performs destructive
% alignment of eye positions, based on (...,'alignBuffer',[tPre tPost])
% Eye positions from tPre-before to tPost-after the motion trajectory are
% saved.
% 
% % sample usage
% subj = 'NP'; para = 'percPurs_Discrim'; modTwi = 250; 
% path = 'C:\data\Neurostim\';
% fileOut = [subj '.' para num2str(modTwi) '.mat'];
% pathOut = 'c:\data\percPurs';
%  cicprms = {'nrTrials',[],{'ge',45}; {'traj','modTwi'},{'trial',1,'atTrialTime',inf},{'eq',modTwi}};
%  [files, paths]=percPurs_FileList(path, 'paradigm',para,'subject',subj,'cicprms',cicprms)
% % see separate help for setup of percPurs_FileList
% msg = percPurs_Merge(files, fileOut, 'paths', paths, 'pathOut', pathOut, 'align', true,'alignBuffer', [-200 200]); 
% 
% See also percPurs_FileList

p = inputParser;
addRequired(p,'files'); %,@ischar);
addRequired(p,'fileOut', @ischar);
addOptional(p,'paths',[]);
addOptional(p,'pathOut',[]);
addParameter(p,'align',false,@islogical); % if true, re-aligns all eye data to a consistent length, based on the alignBuffer relative to the trajectory
addOptional(p,'alignBuffer',[-200 200],@isnumeric);

p.KeepUnmatched = true;
parse(p,files,fileOut, varargin{:})
params = p.Results;

d = marmodata.mdbase('path',params.paths,'file',params.files,'loadArgs',{'loadEye',true});

%% collate behaviour
[time,trial,frame,keyTmp] = d.meta.choice.keyIx('time',Inf);
key = cell2mat(keyTmp);
ignoreTrial = isnan(key);
keepInd = find(~ignoreTrial); 
psychData.time = time(~ignoreTrial);
psychData.trial = trial(~ignoreTrial);
psychData.frame = frame(~ignoreTrial);
psychData.key = key(~ignoreTrial);
psychData.resp = 2-psychData.key; % 1=faster; 0=slower
nTrial = length(psychData.resp);

%% collate stimulus

k=d.meta.cic.screen.data;
stim.frRate = k{1}.frameRate;

% ## check carefully if we want traj or targ for this
% There's usually a frame difference, but sometimes much more
% Targ is later, so we'll go with that
stim.tStTraj = 1000*(d.meta.traj.startTime('trial',keepInd).time - d.meta.cic.firstFrame('trial',keepInd).time); % start time in ms
stim.tSt = 1000*(d.meta.targ.startTime('trial',keepInd).time - d.meta.cic.firstFrame('trial',keepInd).time); % start time in ms

stim.tDur = d.meta.traj.duration('time',Inf','trial',1).data;
%stim.tMod = d.meta.traj.modTwi('time',Inf','trial',1).data;
%stim.modGain=d.meta.traj.modGain('time',Inf, 'trial', keepInd).data; %#ok<*FNDSB>
 
%% Reconstruct stimulus
%stim.modTraj = d.meta.traj.modTraj('time',Inf,'trial',keepInd).data; 
% cell array if rows are different length :( 
% if iscell(stim.modTraj)
%     len=cellfun(@length,stim.modTraj);
%     stim.modTraj = cell2nanmat(stim.modTraj);
%     if range(len)<=1 % single-frame rounding issue - no effect in modulation at middle of sequence
%         stim.modTraj(:,end) = []; % remove last entry for each trial
%     end
% end

% stim.posStartX = d.meta.traj.X('time',0,'trial',1).data; 
% stim.posStartY = d.meta.traj.Y('time',0,'trial',1).data;
% %stim.rampSp = d.meta.traj.rampSp('trial',1).data;
% 
% stim.posTraj = d.meta.traj.sigTraj('time',Inf, 'trial', 1).data;
% stim.posTraj = reshape(stim.posTraj, 2, length(stim.posTraj)); %temporary solution to posTraj fix array 

stim.posTraj = d.meta.traj.sigTraj('time',Inf, 'trial', psychData.trial).data; % 25Apr - changed this so that we get position trajectory for every trial

%stim.spTraj = targXY.sumSinesSpeed();
nFrame = length(stim.posTraj);
% stim.tTraj = (1:(1+nFrame)) / (stim.frRate/1000);
stim.tTraj = (1:nFrame) / (stim.frRate/1000); % 25Apr2022 - removed extra index

%% get stimulus parameters for each trial
stim.fundList = d.meta.traj.fund('trial',keepInd,'time',Inf).data;
stim.freqMultList = d.meta.traj.freqMult('trial',keepInd,'time',Inf).data;
stim.xAmpList = d.meta.traj.xAmp('trial',keepInd,'time',Inf).data;
stim.yAmpList = d.meta.traj.yAmp('trial',keepInd,'time',Inf).data;
stim.xPhList = d.meta.traj.xPh('trial',keepInd,'time',Inf).data;
stim.yPhList = d.meta.traj.yPh('trial',keepInd,'time',Inf).data;
stim.condIds = d.condIds(keepInd);
stim.numConds = d.numConds;


%% Reconstruct stimulus separately for each trial, taking into account framedrops
ind = 0;
for a = psychData.trial % loop through trials
    ind = ind+1; % index for saving "complete" trials
    stim.posTrajFrDrop{ind} = squeeze(stim.posTraj(ind,:,:));
    frDrop = d.meta.cic.frameDrop('trial',a).data; % 2 columns - [frame, dropDuration (ms)]
    
    % where there are framedrops, we replicate the stimulus position for
    % the relevant number of frames. 
    for b = 1:size(frDrop,1)
       startFr = frDrop(b,1);
       % ignore framedrops after trajectory ends (is this OK - we don't reference to when trajectory starts!?)
       if isnan(startFr) || startFr>length(stim.posTrajFrDrop{ind}), continue, end 
       nFrDrop = round(frDrop(b,2) / (1/stim.frRate)); % how many frames were dropped?
       stim.posTrajFrDrop{ind} = [stim.posTrajFrDrop{ind}(:,1:startFr) stim.posTrajFrDrop{ind}(:,startFr)*ones(1,nFrDrop) stim.posTrajFrDrop{ind}(:,(startFr+1):end)];
    end
end

% for a = psychData.trial % loop through trials
%     ind = ind+1;
%     stim.posTrajFrDrop{a} = squeeze(stim.posTraj(ind,:,:));
%     frDrop = d.meta.cic.frameDrop('trial',a).data; % 2 columns - [frame, dropDuration (ms)]
%     
%     % where there are framedrops, we replicate the stimulus position for
%     % the relevant number of frames. 
%     for b = 1:size(frDrop,1)
%        startFr = frDrop(b,1);
%        % ignore framedrops after trajectory ends (is this OK - we don't reference to when trajectory starts!?)
%        if isnan(startFr) || startFr>length(stim.posTrajFrDrop{a}), continue, end 
%        nFrDrop = round(frDrop(b,2) / (1/stim.frRate)); % how many frames were dropped?
%        stim.posTrajFrDrop{a} = [stim.posTrajFrDrop{a}(:,1:startFr) stim.posTrajFrDrop{a}(:,startFr)*ones(1,nFrDrop) stim.posTrajFrDrop{a}(:,(startFr+1):end)];
%     end
% end

%% Re-align eye (destructive)
ey = d.eye(keepInd);

if params.align
    eyeData.t = params.alignBuffer(1):(stim.tDur+params.alignBuffer(2));
    for a = 1:nTrial % everything is binned here. We might have up to 1 ms of slop, but that is acceptable for eye positions
        tt = eyeData.t - round(1000*ey(a).toffset) + round(stim.tSt(a));

%         ey(a).toffset = 0;
        % re-align x,y,pwdth,phght
        % Note that parea is updated automatically
        for str = {'x','y','pwdth','phght'}
            try
            eyeData.(str{:})(a,:) = ey(a).(str{:})(tt);
            catch,keyboard,end
        end
    end
    
else
    for a = 1:nTrial
        for str = {'t','x','y','pwdth','phght'}
            eyeData.(str{:}){a} = ey(a).(str{:});
        end
    end
end
%% Speed of stimulus
stim.spTraj = findSpeed_sumSines(stim)
%% Putting other variables in order of condition
stim = pursuit2D_condList(stim)

%% Export
if ~isempty(params.pathOut)
    fileOut = [params.pathOut filesep params.fileOut];
else
    fileOut = params.fileOut;
end
save(fileOut, 'stim', 'psychData','eyeData','params')
msg = sprintf('Saved %d trials to %s\n', nTrial, fileOut);
