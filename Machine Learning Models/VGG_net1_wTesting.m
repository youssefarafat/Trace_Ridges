imageFolder = 'MATLAB/TMA/Fibre_Patches_noBack'; % Replace with the path to your image folder
imds = imageDatastore(imageFolder, 'FileExtensions', '.jpg'); % Adjust extension if needed

% Assign labels based on filename
numImages = numel(imds.Files);
labels = strings(numImages, 1);

for i = 1:numImages
    fileName = imds.Files{i};
    disp(i)
    
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
splitRatio = 0.8; 
[trainImgs, testImgs] = splitEachLabel(imds, splitRatio, 'randomize');

%%
% Use imageDatastore directly for training and testing images
augmentedTrainImgs = trainImgs;
augmentedTestImgs = testImgs;
%%
net = vgg16;
%%
layers = net.Layers;
numClasses = numel(categories(trainImgs.Labels));

% Modify the last fully connected layer
layers(end-2) = fullyConnectedLayer(numClasses, 'WeightLearnRateFactor', 10, 'BiasLearnRateFactor', 10);
layers(end) = classificationLayer;
%%
options = trainingOptions('sgdm', ...
    'InitialLearnRate', 1e-4, ...
    'MaxEpochs', 10, ...
    'MiniBatchSize', 16, ...
    'Plots', 'training-progress');
%%
trainedNet = trainNetwork(augmentedTrainImgs, layers, options);
%%
% Step 1: Classify the test data
predictedLabels = classify(trainedNet, augmentedTestImgs);

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

