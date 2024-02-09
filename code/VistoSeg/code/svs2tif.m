% Add the path to the BioFormats MATLAB Toolbox
addpath(genpath('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/VistoSeg/code/bfmatlab'))

% Specify the path to your .svs file
dt = '/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/images/';
ot = '/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Images/VistoSeg/';
fname = 'V13Y24-346_40x.svs';
% Open the .svs file using BioFormats
reader = bfGetReader(fullfile(dt,fname));

% Get image size
sizeX = reader.getSizeX();
sizeY = reader.getSizeY();
numPlanes = reader.getImageCount();