% Define image directory and filename
imageDir = '/MATLAB/UNET/';
imageName = 'A+_N.tif';
imagePath = fullfile(imageDir, imageName);

% Read the image
image = imread(imagePath);
%image = rgb2gray(image); %colour to gray image
% Define patch size
patchSize = [32, 32];

% Create folder to save patches
saveDir = fullfile(imageDir, 'AN_patches');
mkdir(saveDir);

% Extract and save patches
patchIndex = 1;
for r = 1:patchSize(1):size(image, 1)-patchSize(1)+1
    for c = 1:patchSize(2):size(image, 2)-patchSize(2)+1
        % Print loop bounds for debugging
        fprintf('r: %d - %d, c: %d - %d\n', r, r+patchSize(1)-1, c, c+patchSize(2)-1);
        
        % Calculate end row and column indices for the patch
        endRow = r + patchSize(1) - 1;
        endCol = c + patchSize(2) - 1;
        
        % Extract patch
        patch = image(r:endRow, c:endCol, :);
        
        % Save patch
        patchFileName = sprintf('AN_patch_%05d.png', patchIndex);
        imwrite(patch, fullfile(saveDir, patchFileName));
        
        % Increment patch index
        patchIndex = patchIndex + 1;
    end
end
%%
%Define GT
groundTruthFile = 'groundTruth_A+_N_image.mat';

% Load ground truth data
load(groundTruthFile); % This will load the variable 'groundTruth' from the MAT file
groundTruth = max(groundTruth(:,:,:),[],3);
groundTruth = ((groundTruth>0==1)+1);
groundTruth = imdilate(groundTruth,ones(5));
% Create folder to save labels
labelsDir = fullfile(imageDir, 'AN_Dilated_labels');
mkdir(labelsDir);

% Extract and save label patches
patchIndex = 1;
for r = 1:patchSize(1):size(groundTruth, 1)-patchSize(1)+1
    for c = 1:patchSize(2):size(groundTruth, 2)-patchSize(2)+1
        % Calculate end row and column indices for the patch
        endRow = r + patchSize(1) - 1;
        endCol = c + patchSize(2) - 1;
        
        % Extract label patch
        AN_labelPatch = uint8(groundTruth(r:endRow, c:endCol));
        
        % Save label patch as grayscale PNG image
labelFileName = sprintf('AN_Dilated_label_%05d.png', patchIndex);
imwrite(AN_labelPatch, fullfile(labelsDir, labelFileName));


        % Increment patch index
        patchIndex = patchIndex + 1;
    end
end
%%
%Top half for training, bottom half for testing.

% Create folder to save training images
trainingImagesDir = fullfile(imageDir, 'AN_trainingImages');
mkdir(trainingImagesDir);

% Extract and save patches from the top half of the image for training
patchIndex = 1;
for r = 1:patchSize(1):512-patchSize(1)+1 % Adjusted loop indices for top half
    for c = 1:patchSize(2):1024-patchSize(2)+1 % Adjusted loop indices for top half
        % Calculate end row and column indices for the patch
        endRow = r + patchSize(1) - 1;
        endCol = c + patchSize(2) - 1;
        
        % Extract patch
        patch = image(r:endRow, c:endCol, :);
        
        % Save patch
        patchFileName = sprintf('AN_patch_%05d.png', patchIndex);
       
            imwrite(patch, fullfile(trainingImagesDir, patchFileName));
        
        
        % Increment patch index
        patchIndex = patchIndex + 1;
    end
end


%%
% Create folder to save training labels
trainingLabelsDir = fullfile(imageDir, 'AN_Dilated_trainingLabels');
mkdir(trainingLabelsDir);

% Extract and save label patches for training
patchIndex = 1;
for r = 1:patchSize(1):512-patchSize(1)+1 % Adjusted loop indices for top half
    for c = 1:patchSize(2):1024-patchSize(2)+1 % Adjusted loop indices for top half
        % Calculate end row and column indices for the patch
        endRow = r + patchSize(1) - 1;
        endCol = c + patchSize(2) - 1;
        
        % Extract label patch
        AN_labelPatch = uint8(groundTruth(r:endRow, c:endCol));
        
        labelFileName = sprintf('AN_Dilated_label_%05d.png', patchIndex);
        
        % Save label patch for training (using top half)
        imwrite(AN_labelPatch, fullfile(trainingLabelsDir, labelFileName));
        
        % Increment patch index
        patchIndex = patchIndex + 1;
    end
end


%%
%Augmented Images and Labels

% Define folder paths
trainingImagesDir = fullfile(imageDir, 'AN_trainingImages');
trainingLabelsDir = fullfile(imageDir, 'AN_trainingLabels');
augmentedImagesDir = fullfile(imageDir, 'AN_augmentedImages');
augmentedLabelsDir = fullfile(imageDir, 'AN_augmentedLabels');

% Load all images from the trainingImages directory
imageFiles = dir(fullfile(trainingImagesDir, '*.png'));

% Create augmented images and labels directories if they don't exist
if ~exist(augmentedImagesDir, 'dir')
    mkdir(augmentedImagesDir);
end
if ~exist(augmentedLabelsDir, 'dir')
    mkdir(augmentedLabelsDir);
end

% Loop through each image
for i = 1:length(imageFiles)
    % Read the image
    image = imread(fullfile(trainingImagesDir, imageFiles(i).name));
    
    % Apply data augmentation
    augmentedImages = cell(1, 5); % Store augmented images (including original)
    augmentedImages{1} = image; % Store original image
    
    % Data augmentation: horizontal flip
    augmentedImages{2} = flip(image, 2);
    
    % Data augmentation: vertical flip
    augmentedImages{3} = flip(image, 1);
    
    % Data augmentation: rotation (90 degrees)
    augmentedImages{4} = imrotate(image, 90);
    
    % Data augmentation: Gaussian noise
    noise = 0.05 * randn(size(image)); % Adjust the noise level as needed
    noise_var = var(noise(:)); % Compute variance of the noise vector
    augmentedImages{5} = imnoise(image, 'gaussian', 0, noise_var);
    
    % Save augmented images
    for j = 1:numel(augmentedImages)
        [~, name, ext] = fileparts(imageFiles(i).name);
        augmentedImageName = sprintf('%s_augmented_%d%s', name, j, ext);
        imwrite(augmentedImages{j}, fullfile(augmentedImagesDir, augmentedImageName));
    end
    
%     % Load corresponding label
%     labelFileName = fullfile(trainingLabelsDir, ['label_' num2str(i) '.mat']);
%    disp(labelFileName);
%     load(labelFileName, 'labelPatch');
%     
%     % Apply the same data augmentation to the label
%     augmentedLabels = cell(1, 5); % Store augmented labels (including original)
%     augmentedLabels{1} = labelPatch; % Store original label
%     
%     % Apply the same data augmentation to the label
%     for j = 2:numel(augmentedImages)
%         augmentedLabels{j} = flip(labelPatch, 2);
%     end
%     
%     % Save augmented labels
% for j = 1:numel(augmentedLabels)
%     augmentedLabelName = sprintf('label_%d_augmented_%d.mat', i, j);
%     labelToSave = augmentedLabels{j};
%     save(fullfile(augmentedLabelsDir, augmentedLabelName), 'labelToSave');
% end

end

%%
% Define folder paths
trainingLabelsDir = fullfile(imageDir, 'AN_Dilated_trainingLabels');
augmentedLabelsDir = fullfile(imageDir, 'AN_Dilated_augmentedLabels');

% Load all label files from the trainingLabels directory
labelFiles = dir(fullfile(trainingLabelsDir, '*.png'));

% Create augmented labels directory if it doesn't exist
if ~exist(augmentedLabelsDir, 'dir')
    mkdir(augmentedLabelsDir);
end

% Loop through each label file
for i = 1:length(labelFiles)
    % Extract label number from filename
    [~, fileName, ~] = fileparts(labelFiles(i).name);
    labelNumber = str2double(extractAfter(fileName, "label_"));
    
    % Display the label file path
    labelFileName = fullfile(trainingLabelsDir, labelFiles(i).name);
    disp(labelFileName);
    
    % Load corresponding label
    AN_labelPatch = imread(labelFileName);
    
    % Convert labelPatch to double
    AN_labelPatch = double(AN_labelPatch);
    
    % Apply data augmentation to the label
    augmentedLabels = cell(1, 5); % Store augmented labels (including original)
    augmentedLabels{1} = AN_labelPatch; % Store original label
    
    % Data augmentation: horizontal flip
    augmentedLabels{2} = flip(AN_labelPatch, 2);
    
    % Data augmentation: vertical flip
    augmentedLabels{3} = flip(AN_labelPatch, 1);
    
    % Data augmentation: rotation (90 degrees)
    augmentedLabels{4} = imrotate(AN_labelPatch, 90);
    
    % Data augmentation: Gaussian noise
    %noise = 0.05 * randn(size(labelPatch));
    %noise_var = var(noise(:)); 
    augmentedLabels{5} = AN_labelPatch;
    
    % Save 
    for j = 1:numel(augmentedLabels)
        augmentedLabelName = sprintf('AN_Dilated_label_%05d_augmented_%d.png', labelNumber, j);
        augmentedLabel = augmentedLabels{j}; % Temporary variable to hold the augmented label
        imwrite(uint8(augmentedLabel), fullfile(augmentedLabelsDir, augmentedLabelName));
    end
end
%%
%Remember to reload variable image
% Create folder to save test images
testImagesDir = fullfile(imageDir, 'AN_Dilated_testImages');
mkdir(testImagesDir);

% Create folder to save test labels
testLabelsDir = fullfile(imageDir, 'AN_Dilated_testLabels');
mkdir(testLabelsDir);

% Extract and save patches for testing
patchIndex = 513;
for r = 513:patchSize(1):1024-patchSize(1)+1
    for c = 1:patchSize(2):1024-patchSize(2)+1
        % Calculate end row and column indices for the patch
        endRow = r + patchSize(1) - 1; % Ensure endRow does not exceed image dimensions
        endCol = c + patchSize(2) - 1; % Ensure endCol does not exceed image dimensions
        
        % Extract patch
        patch = image(r:endRow, c:endCol, :);
        
        % Save patch for testing (using bottom half)
        patchFileName = sprintf('AN_Dilated_patch_%05d.png', patchIndex);
        imwrite(patch, fullfile(testImagesDir, patchFileName));
        
        % Extract label patch
        AN_labelPatch = uint8(groundTruth(r:endRow, c:endCol));
        
        % Save label patch for testing (using bottom half)
        labelFileName = sprintf('AN_Dilated_label_%05d.png', patchIndex);
        imwrite(AN_labelPatch, fullfile(testLabelsDir, labelFileName))
        
        % Increment patch index
        patchIndex = patchIndex + 1;
    end
end
