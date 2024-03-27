% Add the path to the BioFormats MATLAB Toolbox
cd '/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/VistoSeg/'
addpath(genpath('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/VistoSeg/code/bfmatlab'))

% Specify the path to your .svs file
dt = '/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/images/';
ot = '/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Images/VistoSeg/';
fname = 'V13Y24-345_40X.svs';
%fname = 'V13Y24-346_40x.svs';
% Open the .svs file using BioFormats
reader = bfGetReader(fullfile(dt,fname));

% Get image size
sizeX = reader.getSizeX();
sizeY = reader.getSizeY();
numPlanes = reader.getImageCount();
array = {'A1','B1','C1','D1'};
% Loop through chunks and read the image

A=1;
for startX = 1:sizeX/4:sizeX
    endX = min(startX + sizeX/4 - 1, sizeX);
    width = endX - startX + 1;
        img = [];
       for startY = 1:4000:sizeY
        endY = min(startY + 4000 - 1, sizeY);
        height = endY - startY + 1;
        imageData1 = bfGetPlane(reader, 1, startX, startY, width, height);
        imageData2 = bfGetPlane(reader, 2, startX, startY, width, height);
        imageData3 = bfGetPlane(reader, 3, startX, startY, width, height);
        img(startY:endY, 1:width,  :) = cat(3,imageData1,imageData2,imageData3);
       end
          %  Read the chunk
            % imageData1 = bfGetPlane(reader, 1, startX, 1, width, sizeY);
            % imageData2 = bfGetPlane(reader, 2, startX, 1, width, sizeY);
            % imageData3 = bfGetPlane(reader, 3, startX, 1, width, sizeY);
            % Concatenate the chunk to the output image
            % img = uint8(cat(3,imageData1,imageData2,imageData3));

        img = imresize(uint8(img), 0.7);
        save(fullfile(ot,[fname(1:end-7),array{A},'.mat']),'img','-v7.3'); 
        imwrite(img, fullfile(ot,[fname(1:end-7),array{A},'.tif']))
        disp([array{A}, ' done'])
       A=A+1;
end
