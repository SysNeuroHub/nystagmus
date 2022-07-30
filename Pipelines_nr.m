%% Pipeline 1 
% processing raw data
% input matching subject name and date of raw data

pathBase0 = '/home/marmolab/data/2022/'; %pipeline 1 inputs
%pathBase1 = '/home/marmolab/data/pursuit/pipeline1/'; %where pipeline 1 saves out
pathBase1 = '/home/marmolab/data/nystagmus/pipeline1/'; %where pipeline 1 saves out


subjList = {'sameDir_glasses'};
pathInList = {'07/30'};
%weekList = [1 2 1 2 1 2 1 2 1 2]; 

for a = 1:length(subjList)
    disp(['Loading ' subjList{a} ' ' pathInList{a}])
    pursuit2D_stimType_pipeline_nr('subj',subjList{a},'para','pursuit2D',...
        'pathIn',[pathBase0 pathInList{a}],'pathOut',[pathBase1]);
end

%% Pipeline 1.2
% renaming files to add week IDs

% use second newName to return file names to original output from pipeline 1
% useful after accidentally rerunning this section and producing "week1.week1.mat"

fileFolder = dir(fullfile([pathBase1 ,'*.mat']));

for a = 1:numel(fileFolder)
    oldName = [fileFolder(a).folder filesep fileFolder(a).name];
    [~,f,ext] = fileparts(fileFolder(a).name);
    
    newName = strcat(f,'.Week',num2str(weekList(a)),ext);
%     newName = strcat(extractBefore(fileFolder(a).name,'.Week'),ext);
    
    movefile(oldName,[fileFolder(a).folder filesep newName]);
end

disp(['Processed ' num2str(numel(fileFolder)) ' files'])

%% Pipeline 2
% performing data analysis, adding metrics and variables
% loads all files in pathBase2, processes, and saves in pathBase 3

% NOTE: does NOT work if pathBase1 + 2 are the same and it's already been
% run before! need to delete the file previously created

pathBase1 = '/home/marmolab/data/pursuit/pipeline1/'; %where pipeline 1 saves
pathBase2 = '/home/marmolab/data/pursuit/'; %where pipeline 2 saves out

fileFolder2 = dir(fullfile(pathBase1,'*.mat')); 

for a = 1:numel(fileFolder2)
    fileName = fileFolder2(a).name;
    load([pathBase1 fileName])
    disp(['Loaded file ' num2str(a) ' = '  fileName])
    
    [~,f,ext] = fileparts(fileFolder2(1).name);
    
    [eyeData,eyeResamp,position,saccade,speed] = analysis_pipeline(stim,eyeData,noSacc);
    
    save([pathBase2 f '.P2' ext])
    disp(['Processed ' fileName])
end

disp(['Processed ' num2str(numel(fileFolder2)) ' files'])
%% Pipeline 3
% collating all data into two big files
% requires manually sorting files from Pipeline 2 into week 1 and 2





