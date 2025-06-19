imageSize = [32, 32];
dataDir = '\MATLAB\UNET\trainImages\';
labelDir = '\MATLAB\UNET\trainLabels\';

%%
class_Names = ["BACK","FIB"];

labelIDs = (1:numel(class_Names));
imds = imageDatastore(dataDir);
pxds = pixelLabelDatastore(labelDir,class_Names,labelIDs);

pximds = pixelLabelImageDatastore(imds,pxds)
%%
encoderDepth = 3;
numClasses = numel(class_Names);
[lgraph, OutputSize] = unetLayers(imageSize, numClasses, 'EncoderDepth', encoderDepth);
%%
options = trainingOptions('adam', ...
    'InitialLearnRate', 0.0001, ...
    'MaxEpochs', 15, ...
    'VerboseFrequency', 5, ...
    'MiniBatchSize', 16, ...
    'Shuffle', 'every-epoch', ...
    'Verbose', true);
%%
Unet = trainNetwork(pximds, lgraph, options)
%%
save('UNET_top_Fibre_1.mat', 'Unet');