% Define image folder
imageFolder = '/MATLAB/TMA/Patches_ViT'; % Update this path

% Create an image datastore
imds = imageDatastore(imageFolder, 'FileExtensions', '.jpg'); % Adjust extension if needed

% Assign labels based on filename
numImages = numel(imds.Files);
labels = strings(numImages, 1);

for i = 1:numImages
    disp(i)
    fileName = imds.Files{i};
    if contains(fileName, 'BrNormal')
        labels(i) = "Normal";
    elseif contains(fileName, 'BrTUM')
        labels(i) = "TUM";
    else
        warning('Unknown label for image: %s', fileName);
    end
end

% Convert the string labels to categorical type
imds.Labels = categorical(labels);
%%
% Split the dataset into training and testing sets
splitRatio = 0.8; 
[trainImgs, testImgs] = splitEachLabel(imds, splitRatio, 'randomize');

% Define image input size
inputSize = [384 384 3]; % Image size
%%
% Load the pre-trained Vision Transformer model
net = visionTransformer("tiny-16-imagenet-384"); % Load the model

% Convert to layer graph to inspect layers
lgraph = layerGraph(net); % Convert the network to a layer graph

% Display layer names for debugging
layerNames = {lgraph.Layers.Name}; % Get the names of all layers
disp(layerNames); % Display all layer names
%%
% Define the new fully connected layer and classification layer
numClasses = 2; % Number of classes

% Create the new fully connected layer
newFC = fullyConnectedLayer(numClasses, ...
    'WeightLearnRateFactor', 10, ...
    'BiasLearnRateFactor', 10, ...
    'Name', 'new_fc'); % Name for new fully connected layer
%%
% Create the new softmax layer
newSoftmax = softmaxLayer('Name', 'new_softmax'); % New softmax layer
newClassLayer = classificationLayer('Name', 'new_classification'); % New classification layer
%%
% Remove the 'cls_index', 'head', and 'softmax' layers
lgraph = removeLayers(lgraph, {'cls_index', 'head', 'softmax'}); % Remove layers by name
%%
% Add the new layers to the graph
lgraph = addLayers(lgraph, newFC);
lgraph = addLayers(lgraph, newSoftmax);
lgraph = addLayers(lgraph, newClassLayer);
%%
% Connect the new layers
lastLayerName = 'encoder_norm'; % The layer before removed layers
lgraph = connectLayers(lgraph, lastLayerName, newFC.Name); % Connect to the new FC layer
lgraph = connectLayers(lgraph, newFC.Name, newSoftmax.Name); % Connect the new FC layer to the new softmax layer
lgraph = connectLayers(lgraph, newSoftmax.Name, newClassLayer.Name); % Connect softmax to classification layer

%%
% Define data augmentation options
imageAugmenter = imageDataAugmenter( ...
    'RandXReflection', true, ...         % Randomly flip images horizontally
    'RandRotation', [-20, 20], ...       % Randomly rotate images between -20 and +20 degrees
    'RandXTranslation', [-10, 10], ...    % Randomly translate images horizontally
    'RandYTranslation', [-10, 10], ...    % Randomly translate images vertically
    'RandScale', [0.8, 1.2]);            % Randomly scale images

% Create an augmented image datastore for training data
augmentedTrainImgs = augmentedImageDatastore(inputSize(1:2), trainImgs, ...
    'DataAugmentation', imageAugmenter);
%%
options = trainingOptions('adam', ... % Use Adam optimizer
    'MaxEpochs', 15, ...             % Set maximum number of epochs
    'InitialLearnRate', 1e-4, ...    % Set initial learning rate
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 3, ...
    'LearnRateDropFactor', 0.5, ...
    'MiniBatchSize', 32, ...          % Set mini-batch size
    'Plots', 'training-progress', ...  % Plot training progress
    'Verbose', false);                 % Suppress output 
%%
trainedNet = trainNetwork(augmentedTrainImgs, lgraph, options);
%%
% After training, evaluate the model on test images
%testPredictions = classify(trainedNet, testImgs); % Classify test images
testLabels = testImgs.Labels; % Get the actual labels for test images

% Calculate accuracy
%accuracy = sum(testPredictions == testLabels) / numel(testLabels);
%disp(['Test Accuracy: ', num2str(accuracy * 100), '%']);

% Assuming testImgs is an ImageDatastore
expectedSize = [384, 384, 3];

% Use augmentedImageDatastore to resize images to 384x384
testImgsResized = augmentedImageDatastore(expectedSize(1:2), testImgs);

% Now classify the resized images
predictedLabels = classify(trainedNet, testImgsResized);

% Step 1: Classify the test data
%predictedLabels = classify(trainedNet, testImgs);

% Step 2: Calculate accuracy
accuracy = sum(predictedLabels == testImgs.Labels) / numel(testImgs.Labels);
disp(['Test Accuracy: ', num2str(accuracy)]);

% Step 3: Generate confusion matrix
trueLabels = testImgs.Labels;
confMat = confusionmat(trueLabels, predictedLabels);
disp('Confusion Matrix:');
disp(confMat);

% Step 4: Plot confusion matrix for better visualization
figure;
confusionchart(trueLabels, predictedLabels);
title('Confusion Matrix');

% Step 5: Compute precision, recall, and F1-score
% Extract the values from the confusion matrix
TP = diag(confMat); % True Positives (for each class)
FP = sum(confMat, 1)' - TP; % False Positives (for each class)
FN = sum(confMat, 2) - TP; % False Negatives (for each class)

% Precision: TP / (TP + FP)
precision = TP ./ (TP + FP);

% Recall (Sensitivity): TP / (TP + FN)
recall = TP ./ (TP + FN);

% F1 Score: 2 * (Precision * Recall) / (Precision + Recall)
F1 = 2 * (precision .* recall) ./ (precision + recall);

% Display results
disp('Precision for each class:');
disp(precision);
disp('Recall for each class:');
disp(recall);
disp('F1 Score for each class:');
disp(F1);

