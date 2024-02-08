% Add the path to the BioFormats MATLAB Toolbox
addpath(genpath('bfmatlab'))
addpath(genpath('/Users/madhavitippani/Downloads/bfmatlab'))

% Specify the path to your .svs file
svsFilePath = '/Volumes/Neural_Plasticity/Molecular_Profiling/Amygdala/Human/Visium/Br2743/H_and_E/V13Y24-345_40X.svs';

% Open the .svs file using BioFormats
reader = bfGetReader(svsFilePath);

% Get image size
sizeX = reader.getSizeX();
sizeX = sizeX/4;
sizeY = reader.getSizeY();
numPlanes = reader.getImageCount();

% Set chunk size (adjust as needed)
chunkSizeX = 4000;
chunkSizeY = 4000;

% Initialize output image
outputImage = [];

% Loop through chunks and read the image
    for startY = 1:chunkSizeY:sizeY
        endY = min(startY + chunkSizeY - 1, sizeY);
        height = endY - startY + 1;
        
        for startX = 1:chunkSizeX:sizeX
            endX = min(startX + chunkSizeX - 1, sizeX);
            width = endX - startX + 1;
            
            % Read the chunk
            imageData1 = bfGetPlane(reader, 1, startX, startY, width, height);
            imageData2 = bfGetPlane(reader, 2, startX, startY, width, height);
            imageData3 = bfGetPlane(reader, 3, startX, startY, width, height);
            % Concatenate the chunk to the output image
            % outputImage(startY:endY, startX:endX, :) = cat(3,imageData1,imageData2,imageData3);
            outputImage(startY-15000:endY-15000, startX:endX, :) = cat(3,imageData1,imageData2,imageData3);
        end
        disp(startY)
    end

% Save the output image as .tif
save('final.mat','outputImage')
