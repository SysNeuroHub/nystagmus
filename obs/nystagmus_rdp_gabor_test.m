
deg = [0];% 0:15:180;

nystagmus_rdp('rdp','dir1List_b',deg, 'dir1List_r',[0 90], 'speed',5, 'radius',5, ...
    'nrDots',20, 'coherence',1, 'dotSize', 3, ...
    'nRepPerCond',1,'tDur',1e4);

nystagmus_gabor('gabor','dirList_b',0, 'dirList_r',[90], 'speed',5, 'radius',5, ...
    'nRepPerCond',1,'tDur',1e4);
