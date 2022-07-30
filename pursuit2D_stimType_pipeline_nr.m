function pursuit2D_stimType_pipeline(varargin)
%
% Shows how to:
% 1. Define a list of all files (& their paths) associated with
%   - a single user
%   - a single paradigm
%   - a set of stimulus conditions
% 2. Loads files, merges data 
% 3. Calculates eye velocity
% 4. Detects saccades
% At each stage, processed data is saved
%
% TODO - modify output filename so that it includes unique date
%
% NP - 25Apr2022 - Based on pursuit2d_pipeline with modifications to cope with blocking according to stimulus-type

%% DEFAULT PARAMETERS
default.subj = 'JY'; % initials of subject
default.para = 'pursuit2D'; % what paradigm?
default.pathIn = '/home/marmolab/data/2022/07/01'; % where is the neurostim data?
default.pathOut = '/home/marmolab/data/pursuit'; % where will the merged data get saved

default.stimType = 3; % blocking according to stimulus type
default.afterDate = '2022-07-01'; % only include files collected after this date %JY-changing to the date stim type 3 was created
default.minTrials = 1; % minimum number of trials to include a file

default.alignBuffer = [-200 200]; % (ms) when trials are aligned, we grab this much time before and after the trial

% Savitzky-Golay filter parameters for calculating eye velocity
default.ord = 3; 
default.fl = 51; % previously we were using the equivalent of ord=5; fl=51;

% parameters for saccade detection
default.saccRejectFlag = true;
default.spSaccThresh = 16; % (deg/s) speed rejection threshold
% pa.tSaccRejectPrePost = [-200 200]; % (ms) time window either side of speed modulation in which saccades are detected for trial rejection


%% parse arguments...
p = inputParser();
p.KeepUnmatched = true;
p.addParameter('subj',default.subj, @(x) validateattributes(x,{'char'},{'nonempty'}));
p.addParameter('para',default.para);
p.addParameter('pathIn',default.pathIn);
p.addParameter('pathOut',default.pathOut);
p.addParameter('stimType',default.stimType,@(x) validateattributes(x,{'numeric'},{'scalar','nonempty'})); % 0=step-ramp; 1=single sinusoids; 2=sum-of-sinusoids
p.addParameter('afterDate',default.afterDate);
p.addParameter('minTrials',default.minTrials);
p.addParameter('alignBuffer',default.alignBuffer);
p.addParameter('ord',default.ord);
p.addParameter('fl',default.fl);
p.addParameter('saccRejectFlag',default.saccRejectFlag);
p.addParameter('spSaccThresh',default.spSaccThresh);
p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.parse(varargin{:});

pa = p.Results;

%% define some additional parameters
%
dbprms = {'after',pa.afterDate};
% assemble input 'cicprms' for file searching
% cicprms = {'nrTrials',[],{'ge',pa.minTrials}; }; %{'traj','duration'},{'trial',1,'atTrialTime',inf},{'eq',pa.tDur}};
cicprms = {'nrTrials',[],{'ge',pa.minTrials}; {'traj','stimType'},{'trial',1,'atTrialTime',inf},{'eq',pa.stimType}}; %


% output path-file
% fileOut =[pa.subj '.' pa.para '.tDur' num2str(pa.tDur) '.mat'];

% filedate = extractAfter(pa.pathIn,'data/')
% fileOut =[pa.subj '.' pa.para '.stimType' num2str(pa.stimType) '.' datestr(now,30) '.mat'];


%% Go through pipeline

% get list of files
%[files, paths]=percPurs_FileList(pa.pathIn, 'paradigm',pa.para,'subject',pa.subj,'cicprms',cicprms,'dbprms',dbprms);
[files, paths]=percPurs_FileList(pa.pathIn, 'paradigm',pa.para,'subject',pa.subj,'dbprms',dbprms);

filedate = extractAfter(paths{1},'data/');
filedate = strrep(filedate,filesep,'');
fileOut =[pa.subj '.' pa.para '.stimType' num2str(pa.stimType) '.' filedate '.mat'];

% merge files, extract stim/psychData/eyeData and save them into a new data file
% if strcmp(pa.para,'pursuit2D')
    [fullfileOut,msg] = pursuit2D_Merge_nr(files, fileOut, 'paths', paths, 'pathOut', pa.pathOut, 'align', true,'alignBuffer', pa.alignBuffer);
    disp(msg)
% end

% load the new data file, differentiate and append data
load(fullfileOut,'eyeData','stim');
eyeData = percPurs_Diff(eyeData,'order',pa.ord,'framelen',pa.fl);

%% Reject saccades
if pa.saccRejectFlag && strcmp(pa.para,'pursuit2D')
    eyeData = pursuit2D_SaccReject(eyeData, stim, pa.spSaccThresh);
% elseif pa.saccRejectFlag & strcmp(pa.para,'percPurs_Detect')
%     eyeData = percPurs_detectSaccReject(eyeData,stim, pa.tSaccRejectPrePost,true);
end

% spTraj = targXY.sumSinesSpeed(c.traj);

%% Fitting sin to stimulus speed + adding variables fitting sin to stimulus speed

nTrials = length(stim.xAmpList);

for a = 1:nTrials
    noSacc.indRemove = isnan(eyeData.noSaccSpX(a,:)) | isnan(eyeData.noSaccSpY(a,:)); % logical index identifying locations of NaNs in this trial
    noSacc.spX{a} = eyeData.noSaccSpX(a, ~noSacc.indRemove); % save NaN-removed data from this trial into a cell because trials are now different lengths
    noSacc.spY{a} = eyeData.noSaccSpY(a, ~noSacc.indRemove);
    noSacc.t{a} = eyeData.t(~noSacc.indRemove);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    % constrain frequencies (b1/b2 = variables 2 and 5 in list
    if pa.stimType == 1 %a variable in pipeline
%         ft = fittype( 'sin1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [-Inf 0 -Inf];
% opts.StartPoint = [7.3054 0.00317292529083681 1.45965962462462];
        ft = fittype('sin1');
        fr = (stim.fundList(a)*stim.freqMultList(a) * 2*pi) / 1000; % account for 2*pi and milliseconds
        opts.Lower = [0 fr(1) -Inf];
        opts.Upper = [Inf fr(1) Inf];
        [noSacc.fitH{a},noSacc.gofH{a}] = fit(noSacc.t{a}(200:end)', noSacc.spX{a}(200:end)',ft, opts); %stim type 1
        [noSacc.fitV{a},noSacc.gofV{a}] = fit(noSacc.t{a}(200:end)', noSacc.spY{a}(200:end)',ft, opts); 
        noSacc.stimspeedH{a} = stim.xAmpList(a)*stim.freqMultList(a)*stim.fundList(a)*2*pi;
        noSacc.stimspeedV{a} = stim.yAmpList(a)*stim.freqMultList(a)*stim.fundList(a)*2*pi;
        noSacc.gainH{a} = noSacc.fitH{a}.a1./noSacc.stimspeedH{a}; % GAIN
        noSacc.gainV{a} = noSacc.fitV{a}.a1./noSacc.stimspeedV{a};
    elseif pa.stimType == 2
        ft = fittype('sin2');
        fr = (stim.fundList(a)*stim.freqMultList{a} * 2*pi) / 1000; % account for 2*pi and milliseconds
        if length(fr)>1
            opts.Lower = [0 fr(1) -Inf 0 fr(2) -Inf];
            opts.Upper = [Inf fr(1) Inf Inf fr(2) Inf];
        else
            opts.Lower = [0 fr(1) -Inf];
            opts.Upper = [Inf fr(1) Inf];
        end
%         opts.StartPoint = [7.3054 0.00317292529083681 1.45965962462462];
%         [newVer,noSacc.gofH{a}] = fit(noSacc.t{a}(200:end)', noSacc.spX{a}(200:end)',ft,opts); %stim type 2
%         [oldVer,noSacc.gofH{a}] = fit(noSacc.t{a}(200:end)', noSacc.spX{a}(200:end)','sin2'); %stim type 2
        [noSacc.fitH{a},noSacc.gofH{a}] = fit(noSacc.t{a}(200:end)', noSacc.spX{a}(200:end)',ft, opts); %stim type 2
        [noSacc.fitV{a},noSacc.gofV{a}] = fit(noSacc.t{a}(200:end)', noSacc.spY{a}(200:end)',ft, opts);
        noSacc.stimspeedH{a} = stim.xAmpList{a}.*stim.freqMultList{a}.*stim.fundList(a)*2*pi;
        noSacc.stimspeedV{a} = stim.yAmpList{a}.*stim.freqMultList{a}.*stim.fundList(a)*2*pi;
        noSacc.fitgainH{a}(1) = noSacc.fitH{a}.a1./noSacc.stimspeedH{a}(1);
        noSacc.fitgainV{a}(1) = noSacc.fitV{a}.a1./noSacc.stimspeedV{a}(1);
        if length(noSacc.stimspeedH{a})>1 % cope with multiple frequencies
            noSacc.fitgainH{a}(2) = noSacc.fitH{a}.a2./noSacc.stimspeedH{a}(2);
        end
        if length(noSacc.stimspeedV{a})>1
            noSacc.fitgainV{a}(2) = noSacc.fitV{a}.a2./noSacc.stimspeedV{a}(2);
        end
    elseif pa.stimType == 3
        singlesin = [1:9];
        if ismember(stim.condIds(a), singlesin) == 1
            ft = fittype('sin1');
            fr = (stim.fundList(a)*stim.freqMultList{a} * 2*pi) / 1000; % account for 2*pi and milliseconds
            opts.Lower = [0 fr(1) -Inf];
            opts.Upper = [Inf fr(1) Inf];
            [noSacc.fitH{a},noSacc.gofH{a}] = fit(noSacc.t{a}(200:end)', noSacc.spX{a}(200:end)',ft, opts); %stim type 1
            [noSacc.fitV{a},noSacc.gofV{a}] = fit(noSacc.t{a}(200:end)', noSacc.spY{a}(200:end)',ft, opts); 
            noSacc.stimspeedH{a} = stim.xAmpList{a}*stim.freqMultList{a}*stim.fundList(a)*2*pi;
            noSacc.stimspeedV{a} = stim.yAmpList{a}*stim.freqMultList{a}*stim.fundList(a)*2*pi;
            noSacc.fitgainH{a} = noSacc.fitH{a}.a1./noSacc.stimspeedH{a}; % GAIN
            noSacc.fitgainV{a} = noSacc.fitV{a}.a1./noSacc.stimspeedV{a};
        else
            ft = fittype('sin2');
            fr = (stim.fundList(a)*stim.freqMultList{a} * 2*pi) / 1000; % account for 2*pi and milliseconds
            if length(fr)>1
                opts.Lower = [0 fr(1) -Inf 0 fr(2) -Inf];
                opts.Upper = [Inf fr(1) Inf Inf fr(2) Inf];
            else
                opts.Lower = [0 fr(1) -Inf];
                opts.Upper = [Inf fr(1) Inf];
            end
            [noSacc.fitH{a},noSacc.gofH{a}] = fit(noSacc.t{a}(200:end)', noSacc.spX{a}(200:end)',ft, opts);
            [noSacc.fitV{a},noSacc.gofV{a}] = fit(noSacc.t{a}(200:end)', noSacc.spY{a}(200:end)',ft, opts);
            noSacc.stimspeedH{a} = stim.xAmpList{a}.*stim.freqMultList{a}.*stim.fundList(a)*2*pi;
            noSacc.stimspeedV{a} = stim.yAmpList{a}.*stim.freqMultList{a}.*stim.fundList(a)*2*pi;
            noSacc.fitgainH{a}(1) = noSacc.fitH{a}.a1./noSacc.stimspeedH{a}(1);
            noSacc.fitgainV{a}(1) = noSacc.fitV{a}.a1./noSacc.stimspeedV{a}(1);
            if length(noSacc.stimspeedH{a})>1 % cope with multiple frequencies
                noSacc.fitgainH{a}(2) = noSacc.fitH{a}.a2./noSacc.stimspeedH{a}(2);
            end
            if length(noSacc.stimspeedV{a})>1
                noSacc.fitgainV{a}(2) = noSacc.fitV{a}.a2./noSacc.stimspeedV{a}(2);
            end
        end
    end
    % noSacc.fitHa{a} = noSacc.fitH{a}.a1; % parameters of fitted sin - readd if needed
    % noSacc.fitHb{a} = noSacc.fitH{a}.b1;
    % noSacc.fitHc{a} = noSacc.fitH{a}.c1;
    % noSacc.fitVa{a} = noSacc.fitV{a}.a1;
    % noSacc.fitVb{a} = noSacc.fitV{a}.b1;
    % noSacc.fitVc{a} = noSacc.fitV{a}.c1;

end
%% Put trials in cell arrays based on condition

noSacc.fitgainHbyCond = cell(1,stim.numConds); % initialise
for a = 1:nTrials
    vari.condID = stim.condIds(a); % index of stim condition for this trial
    noSacc.fitgainHbyCond{vari.condID} = [noSacc.fitgainHbyCond{vari.condID}; noSacc.fitgainH{a}]; % append this trial to cell for thisCond
    
end

noSacc.fitgainVbyCond = cell(1,stim.numConds); 
for a = 1:nTrials
    vari.condID = stim.condIds(a);
    noSacc.fitgainVbyCond{vari.condID} = [noSacc.fitgainVbyCond{vari.condID}; noSacc.fitgainV{a}];
end

%% Taking the mean of cell arrays

noSacc.AvgfitgainH = cellfun(@mean, noSacc.fitgainHbyCond,'UniformOutput',false); % take mean across rows (within a column)
% tmp = cellfun(@(x) mean(x,2), noSacc.fitgainHbyCond,'UniformOutput',false); % take mean across columns
noSacc.AvgfitgainV = cellfun(@mean, noSacc.fitgainVbyCond,'UniformOutput',false); 

%% Creating array of frequencies based on condition
for a = 1:stim.numConds
    vari.condID2 = find(stim.condIds==a,1,'first');
    if iscell(stim.freqMultList) % will be a cell if we have multiple frequencies
        vari.freq{a} = cell2mat(stim.freqMultList(vari.condID2)).*stim.fundList(vari.condID2)';
    else % single frequency.
        % Need to decide if output var.freq should be a cell array or a
        % vector
        vari.freq(a) = stim.freqMultList(vari.condID2).*stim.fundList(vari.condID2)';
    end
end
%%
save(fullfileOut,'eyeData','pa','noSacc','nTrials','vari','-append')
disp('Appended velocity to eyeData') 
