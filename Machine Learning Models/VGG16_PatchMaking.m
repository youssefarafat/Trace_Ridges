% Define the directory containing the core images
baseDir1 = '/MATLAB/TMA/EligibleData';
% Define the directory to save patches
saveDir = '/MATLAB/TMA/Fibres_Patches';
% ResNet input size
patchSize = [224, 224]; 

% Get list of all files in the directory
fileList = dir(fullfile(baseDir1, '*.png')); % Change the extension if needed

% Initialize patch counters
normalCounter = 1;
tumorCounter = 1;
%%
% Loop through each file
for k = 1:length(fileList)
    disp(k);
    % Read the image
    InputCore = imread(fullfile(baseDir1, fileList(k).name));
    
    % Extract fibers using your CoreReading function
    CoreOutput = CoreReading(InputCore);
    blueFibres = CoreOutput.blueFibres;

     % Ensure blueFibres is in RGB format
    if size(blueFibres, 3) ~= 3
        blueFibres = repmat(blueFibres, [1, 1, 3]); % Convert to 3 channels if necessary
    end
    
    % Calculate the number of patches
    [rows, cols, ~] = size(blueFibres);
    numPatchesRow = floor(rows / patchSize(1));
    numPatchesCol = floor(cols / patchSize(2));
    
    % Loop through and create patches
    for i = 1:numPatchesRow
        for j = 1:numPatchesCol
            % Extract the patch
           patch = blueFibres((i-1)*patchSize(1)+1:i*patchSize(1), (j-1)*patchSize(2)+1:j*patchSize(2), :);
            
            % Determine the type (normal or tumor) and save the patch
            if contains(fileList(k).name, 'Normal', 'IgnoreCase', true)
                patchName = sprintf('BrNormal_patch_%04d.jpg', normalCounter);
                normalCounter = normalCounter + 1;
            elseif contains(fileList(k).name, 'TUM', 'IgnoreCase', true)
                patchName = sprintf('BrTUM_patch_%04d.jpg', tumorCounter);
                tumorCounter = tumorCounter + 1;
            else
                continue; % Skip if the file name does not contain 'normal' or 'TUM'
            end
            
            % Save the patch
            imwrite(patch, fullfile(saveDir, patchName));
        end
    end
end

disp('Patches created and saved successfully.');
