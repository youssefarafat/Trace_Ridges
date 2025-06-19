classNames = ["BACK", "FIB"];
labelIDs = (1:numel(classNames));

% Create image and label datastores for testing
imdsTest = imageDatastore("testingImages_Dilated");

pxdsTest = pixelLabelDatastore("testingLabels_Dilated", classNames, labelIDs);
%%
load('UNET_Dilated_top_Fibre_4.mat');
%%
predictedLabels = semanticseg(imdsTest, Unet, ...
    'MiniBatchSize', 16, ...
    'Verbose', false);
%%
metrics = evaluateSemanticSegmentation(predictedLabels, pxdsTest);

% Display overall dataset-level metrics (accuracy, mean IoU, etc.)
disp(metrics.DataSetMetrics)