%% NOTE FOR ALL EXPS:
% to start a next trial, press 'z' key


%% simplest experiment function
%a script to test nystagmus_gabor.m
nystagmus_gabor('TST','nRepPerCond',3,'debug',true,'tolerance',10,...
    'tDur',5000,'orientation',[0 180]);

%TODO
% - something is not right about the circle mask (square mask works)
% DONE - superimpose red and blue gratings. multiGaborsN seems not working
% - control timing of appearance (flash suppression)
% - add control condition
% FIXED - bug in useMouse=true mode
%       eye position is not visible on screen
%           >     c.cursor = 'arrow';
%       cannot proceed to the next exp
%       eye position not registered until clicking
% FIXED - correct/incorrect sound feedback is not audible
%       > correct.wav > ding.wav
%       > incorrect.wav > cough.wav
% - provide sound feedback when right after the stimulus offset, not the
% start of trial stimulus onset


%% more proper function for experiment
deg = [90 180];
% two eyes receiving the same direction, specified as 'ori1List'
nystagmus_directions_twoEyes('TST','phaseSpeed',10, 'ori1List',deg, ...
    'nRepPerCond',1, 'sigma',2, 'frequency',.5,'tDur',2000,'debug',true);

deg = [90 180];
% one eye receiving (blue) the directions specified as 'ori1List'
% the other eye (red) is fixed as 0deg
nystagmus_directions_oneEye('TST','phaseSpeed',10, 'ori1List',deg, ...
    'nRepPerCond',1, 'sigma',2, 'frequency',.5,'tDur',2000,'debug',true);

