%a script to test nystagmus_gabor.m
%to start a next trial, press 'z' key
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