params.paths = '/home/marmolab/data/2022/07/30/';


%% unidirectional stimuli. without glasses
% nystagmus_gabor('sameDir'); %phaseSpeed=4
%params.files = 'sameDir.pursuit2D.100848.mat';

%% unidirectional downward stimuli. red/blue glasses 
%nystagmus_gabor('sameDir_glasses','phaseSpeed',2)
 params.files = 'sameDir_glasses.pursuit2D.102145.mat';

%% unidirectional upward stimuli. red/blue glasses 
%nystagmus_gabor('sameDir_glasses','phaseSpeed',8, 'orientation',[180 180])
% params.files = 'sameDir_glasses.pursuit2D.123458.mat';

%% bidirectional stimuli red/blue glasses "stare" condition
%nystagmus_gabor('sameDir_glasses','phaseSpeed',8, 'orientation',[0 180])
% params.files = 'biDir.pursuit2D.125154.mat';
%< no nystagmus

%% bidirectional stimuli red/blue glasses. circle aperture 
%nystagmus_gabor('sameDir_glasses','phaseSpeed',8, 'orientation',[0 180])
% params.files = 'biDir.pursuit2D.125842.mat';
%< yet no nystagmus

%% bidirectional stimuli red/blue glasses. larger circle aperture 
%nystagmus_gabor('sameDir_glasses','phaseSpeed',8, 'orientation',[0 180],'sigma',10)
% params.files = 'biDir.pursuit2D.130429.mat';
%< sigma parameter did not work as intended

%% go back to 125842 bidirectional stimuli red/blue glasses "look" condition
% square aperture again
%nystagmus_gabor('sameDir_glasses','phaseSpeed',8, 'orientation',[0 180])
% params.files = 'biDir.pursuit2D.130957.mat';
%< switch in 2nd trial? but gain is very small compared to unidirectional
%condition
% felt like the stimulus was split between upward and downward

%% bidirectional condition
% nystagmus_gabor('biDir','phaseSpeed',4, 'orientation',[0 180], 'sigma',3, 'frequency',.5);
% params.files = 'biDir.pursuit2D.141606.mat';
%< no nystagmus

%% bidirectional condition smaller aperture
% nystagmus_gabor('biDir','phaseSpeed',4, 'orientation',[0 180], 'sigma',2, 'frequency',.5);
%params.files = 'biDir.pursuit2D.142422.mat';
%< maybe 1st trial?


%% bidirectional condition
%nystagmus_gabor('biDir','phaseSpeed',13, 'orientation',[0 180], 'sigma',2, 'frequency',.5);
% params.files = 'biDir.pursuit2D.144030.mat';
%< maybe the 2nd trial?


%% bidirectional condition (same as above) occluded right eye
%nystagmus_gabor('biDir','phaseSpeed',13, 'orientation',[0 180], 'sigma',2, 'frequency',.5);
% params.files = 'biDir-oneEye.pursuit2D.144733.mat';
%nystagmus happen but quite small gain

%% bidirectional condition, larger aperture again
% nystagmus_gabor('biDir','phaseSpeed',13, 'orientation',[0 180], 'sigma',4, 'frequency',.5);
% params.files = 'biDir.pursuit2D.145235.mat';

%% bidirectional condition, parameter from Logothetis 1990, longer duration
%nystagmus_gabor('biDir','phaseSpeed',10, 'orientation',[0 180], 'sigma',2, 'frequency',.5,'tDur',2e4);
% params.files = 'biDir.pursuit2D.150033.mat';
%< somewhat convincing though the last trial is without nystagmus

%% bidirection condition, right eye closed 
%nystagmus_gabor('biDir','phaseSpeed',10, 'orientation',[0 180], 'sigma',2, 'frequency',.5,'tDur',2e4);
%params.files = 'biDir-oneEye.pursuit2D.150802.mat';
%eye tracking is not great as I closed my right eye by pulling the eyelid

%% bidirection condition, right eye closed 
%nystagmus_gabor('biDir','phaseSpeed',10, 'orientation',[0 180], 'sigma',2, 'frequency',.5,'tDur',2e4);
params.files = 'biDir-oneEye.pursuit2D.151834.mat';


d = marmodata.mdbase('path',params.paths,'file',params.files,'loadArgs',{'loadEye',true});


tiledlayout('flow');
for itr = 1:d.numTrials
    nexttile;
    tidx = d.eye(itr).t< 20; 
    plot(d.eye(itr).t(tidx), d.eye(itr).x(tidx), d.eye(itr).t(tidx), d.eye(itr).y(tidx));
end
xlabel('time [s]')
ylabel('eye position [deg]')