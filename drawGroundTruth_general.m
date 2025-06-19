%% Draw ground truth line by line
clf
%% Display image to be delineated
dataIn = imread('/flattened_imageV3.tif');
dataOut = dataIn;
groundTruth = zeros (size(dataIn));
imagesc(dataOut);
k=0;
%% draw lines, keep in a stack 
k=k+1;
lineDrawn                   = drawpolyline;
lineAsMask                  = createMask(lineDrawn);
groundTruth(:,:,k)          = k*lineAsMask;
imagesc(dataIn.*uint8(((1-imdilate(max(groundTruth,[],3),ones(5))))))

%% Use this to measure distance between GT and output

imagesc(((fibronectinOut2.edges>0)).*(bwdist(max(groundTruth,[],3)>0)))
imagesc((bwdist(fibronectinOut2.edges>0)).*(max(groundTruth,[],3)>0))


qqq=regionprops(fibronectinOut2.edges,bwdist(max(groundTruth,[],3)>0),'Area','MeanIntensity');
qqq=regionprops(max(groundTruth,[],3)>0,bwdist(fibronectinOut2.edges>0),'Area','MaxIntensity');
mean([qqq.MeanIntensity])



%%

save('groundTruth_flattened_imageV3.mat','groundTruth');