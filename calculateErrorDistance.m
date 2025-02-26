function [DistanceErrors] = calculateErrorDistance(ActualResults,groundTruth)

% flatten both parameters in case they are a 3D matrix
% Also take values >0 in case they are labelled
ActualResults   = max(ActualResults,[],3)>0;
groundTruth   = max(groundTruth,[],3)>0;


distMap_GT  = bwdist(groundTruth);
distMap_Ac  = bwdist(ActualResults);
DistanceErrors.comparison       = (groundTruth)-(ActualResults);
DistanceErrors.comparison2(:,:,1) = (groundTruth);
DistanceErrors.comparison2(:,:,2) =((groundTruth)+(ActualResults))==2; 
DistanceErrors.comparison2(:,:,3) =(ActualResults);

DistanceErrors.DAG_map          = (ActualResults).*(distMap_GT);
DistanceErrors.DGA_map          = (distMap_Ac).*(groundTruth);

DAG_all          = regionprops(ActualResults,distMap_GT,'Area','MeanIntensity','MaxIntensity');
DGA_all          = regionprops(groundTruth,distMap_Ac,'Area','MeanIntensity','MaxIntensity');

DistanceErrors.DAG_av           = mean([DAG_all.MeanIntensity]);
DistanceErrors.DGA_av           = mean([DGA_all.MeanIntensity]);

DistanceErrors.DAG_avmax        = mean([DAG_all.MaxIntensity]);
DistanceErrors.DGA_avmax        = mean([DGA_all.MaxIntensity]);

DistanceErrors.DAG_max          = max([DAG_all.MaxIntensity]);
DistanceErrors.DGA_max          = max([DGA_all.MaxIntensity]);

DistanceErrors.D_av             = DistanceErrors.DAG_av+DistanceErrors.DGA_av;
DistanceErrors.D_max            = DistanceErrors.DAG_max+DistanceErrors.DGA_max;

DistanceErrors.Tot2              = sum(sum(DistanceErrors.DAG_map)) + sum(sum(DistanceErrors.DGA_map));



