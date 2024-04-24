function [fibronectinOut,fibronectinOut2,dataOut] =  Trace_Ridges(dataIn,cannySize,distGap)


if isa(dataIn,'char')
    % input parameter is a name, read
    try
        dataIn                  = imread(dataIn);
    catch
        disp('Could not read input file')
        return;
    end
end
[rows,cols,levs]                = size(dataIn);
[~,maxChan]                     = max(mean(mean(dataIn)));


% input is the 3D data or the maximum intensity projection OR an RGB
if levs==1  
    dataIn_MIP                  = dataIn;
elseif levs ==3
    % This is an RGB image
    dataIn_MIP                  = rgb2gray(dataIn);
else
    % this is a multiple level image, take Maximum Intensity Projection
    [dataIn_MIP]                = max(dataIn(:,:,:),[],3);
end

%detect the background dark/bright and leave with intensity always as high
[y,x]                           = hist(double(dataIn_MIP(:)),10);

if (sum(y(6:10))/sum(y(1:5)))>2
    dataIn_MIP                  = max(max(dataIn_MIP))-dataIn_MIP;
end

% Detect background once inverted
[y,x]                           = hist(double(dataIn_MIP(:)),50);
[y2,x2,w,p]                     = findpeaks(y);
backgroundIntensityLev          = ceil(x(x2(1))+2);
background                      = bwmorph(dataIn_MIP<=backgroundIntensityLev,'majority');




%distGap                         = 1.5;
if ~exist('distGap','var')
    distGap                     = 1.5;
elseif isempty(distGap)
    distGap                     = 1.5;
end

sizeFilt                    = 3;
% if ~exist('sizeFilt','var')
%     sizeFilt                       = 3;
% elseif isempty(sizeFilt)
%     sizeFilt                       = 3;
% end

if ~exist('cannySize','var')
    cannySize                       = 2;
elseif isempty(cannySize)
    cannySize                       = 2;
end
% if ~exist('displayData','var')
%     displayData                     = 0;
% elseif isempty(displayData)
%         displayData                 = 0;
% end
if ~exist('distGap','var')
    distGap                         = 13;
elseif isempty(distGap)
        distGap                     = 13;
end


%dataIn_MIP                 = imfilter(dataIn_MIP,gaussF(sizeFilt,sizeFilt,1),'replicate');
% lowpass filter
%dataIn_MIP                 = medfilt2(dataIn_MIP,[sizeFilt sizeFilt]);
dataIn_MIP                  = imfilter(dataIn_MIP,fspecial('Gaussian',sizeFilt,1),'replicate');

% find the watershed
allWatersheds               = (watershed(dataIn_MIP))==0;

% find edges that will break the watersheds
allEdges                  = edge(dataIn_MIP,'canny',[],cannySize);



% remove edges from the watersheds
allRidgesNoEdges            = (allWatersheds).*(allEdges==0);
% remove all watersheds that are far from edges (i.e. basins)

distMapEdges                = bwdist(allEdges);


allRidgesNoEdges2           = allRidgesNoEdges.*(distMapEdges<5);
% remove darker regions of the watersheds
lowLevelRidges              = 0.9*255*graythresh(dataIn_MIP);
allRidgesNoEdges3           = allRidgesNoEdges2.*double(dataIn_MIP>(lowLevelRidges));


% SPUR removes part of the actual ridges so it should only be removed in
% the middle bit not in the 
% clean small spurious ridges that are coming from the watershed
%allRidgesNoEdges4           = bwmorph(allRidgesNoEdges3,'spur',5);
allRidgesNoEdges4           = bwskel(allRidgesNoEdges3>0,'MinBranchLength',5);


% remove very small ridges
allRidgesNoEdges4_L         = bwlabel(allRidgesNoEdges3);
allRidgesNoEdges4_P         = regionprops(allRidgesNoEdges4_L,'MajoraxisLength');

% fill holes and thin


% discard anything that is smaller than 5 in major axis
largeRidges1                = ismember(allRidgesNoEdges4_L,find([allRidgesNoEdges4_P.MajorAxisLength]>5));
largeRidges1_L              = bwlabel(largeRidges1);
largeRidges1_P              = regionprops(largeRidges1_L,'EulerNumber');


% fill and thin This is particularly useful with very busy images with
% small holes that appear between parallel fibres. For large holes, this
% can be formed by valid ridges that close a polygon.

% this are ridges with no holes
largeRidges2_NoHoles        = (ismember(largeRidges1_L,find([largeRidges1_P.EulerNumber]==1)));

% find the holes and keep only the smaller ones
largeRidges2_WithHoles      = imfill(ismember(largeRidges1_L,find([largeRidges1_P.EulerNumber]<1)),'holes');
largeRidges2_JustHoles      = largeRidges2_WithHoles.*(1-largeRidges1);
largeRidges2_JustHoles2     = imerode(imclose(largeRidges2_WithHoles,ones(3)),ones(5));

largeRidges3                = largeRidges1.*(1-largeRidges2_JustHoles2);

largeRidges4                = bwmorph(largeRidges3+largeRidges2_NoHoles,'spur',2);


largeRidges4_L              = bwlabel(largeRidges4);
largeRidges4_P              = regionprops(largeRidges4_L,'MajoraxisLength','MinoraxisLength');
%largeRidges3_P              = regionprops(largeRidges3_L,'EulerNumber','Area','MajoraxisLength','MinoraxisLength','Eccentricity');

% keep anything that is straight and long, break otherwise in junctions
largeRidges4A               = ismember(largeRidges4_L,find( ([largeRidges4_P.MajorAxisLength]>10) & ([largeRidges4_P.MinorAxisLength]<3)  )    );
largeRidges4B               = ismember(largeRidges4_L,find( ([largeRidges4_P.MajorAxisLength]>20) & ([largeRidges4_P.MinorAxisLength]<10)  )    );

largeRidges4C               = ismember(largeRidges4_L,find( ([largeRidges4_P.MajorAxisLength]>20) & ([largeRidges4_P.MinorAxisLength]>=10)  )    );
largeRidges4D               = bwmorph(bwmorph(largeRidges4C,'thin','inf'),'branch');
largeRidges4E               = bwmorph(largeRidges4C-largeRidges4D,'spur',1);
largeRidges4E_L             = bwlabel(largeRidges4E);
largeRidges4E_P             = regionprops(largeRidges4E_L,'MajoraxisLength');
largeRidges4F               = ismember(largeRidges4E_L,find( ([largeRidges4E_P.MajorAxisLength]>4)));


% all fibres

allFibres1                   = (largeRidges4F+largeRidges4A+largeRidges4B)>0;
[allFibres1_L,numEdges]      = bwlabel(allFibres1);
allFibres1_P                 = regionprops(allFibres1_L,'Orientation','MajorAxisLength','MinorAxisLength','Centroid','EquivDiameter');

%% Check fibres that would be missed 
% as they are not watershed, single line that finishes in the basin, these
% can be identified from the edges, but only if they are significant.

allEdgesNotFibres           = allEdges.*(1-imdilate(allFibres1,ones(7)));
allEdgesNotFibres2          = imopen(imclose(allEdgesNotFibres,ones(5)),ones(5));
allEdgesNotFibres3          = bwmorph(allEdgesNotFibres2,'thin','inf');
allEdgesNotFibres4          = bwmorph(allEdgesNotFibres3,'branch');
allEdgesNotFibres5          = bwmorph(allEdgesNotFibres3-allEdgesNotFibres4,'spur',1);
[allEdgesNotFibres5_L,numNew]        = bwlabel(allEdgesNotFibres5);
allEdgesNotFibres5_P        = regionprops(allEdgesNotFibres5_L,'MajoraxisLength','MinorAxisLength');


% keep anything that is straight and long, break otherwise in junctions
if numNew<100
    allEdgesNotFibres6          = ismember(allEdgesNotFibres5_L,find( ([allEdgesNotFibres5_P.MajorAxisLength]>15) & ([allEdgesNotFibres5_P.MinorAxisLength]<10)  )    );
else
    allEdgesNotFibres6          = ismember(allEdgesNotFibres5_L,find( ([allEdgesNotFibres5_P.MajorAxisLength]>25) & ([allEdgesNotFibres5_P.MinorAxisLength]<10)  )    );
end

%% Final combination of the fibres

allFibres2                  = bwmorph(bwmorph(allFibres1,'thin','inf'),'branch');
allFibres3                  = bwmorph(allFibres1-allFibres2,'spur',1);

% break corners
points                      = detectHarrisFeatures(allFibres3,'MinQuality',0.45);

for counterP = 1:points.Count
    currentPos              = round(points(counterP).Location);
    rr                      = max(1,-1+currentPos(2)):min(rows,1+currentPos(2));
    cc                      = max(1,-1+currentPos(1)):min(rows,1+currentPos(1));
    try
    allFibres3(rr,cc)       = 0;
    catch
        qqq=1;
    end
end

allFibres3_L                = bwlabel(allFibres3);
allFibres3_P                = regionprops(allFibres3,'MajoraxisLength');
allFibres4                 = ismember(allFibres3_L,find( ([allFibres3_P.MajorAxisLength]>10)));


allFibres                   = allEdgesNotFibres6+allFibres4;
[allFibres_L,numEdges]      = bwlabel(allFibres);
allFibres_P                 = regionprops(allFibres_L,'Orientation','MajorAxisLength','MinorAxisLength','Centroid','EquivDiameter','MaxFeretProperties','Circularity','Area','Eccentricity','Perimeter');

%%
for k = 1:numEdges
    curvature(k) = (allFibres_P(k).MajorAxisLength)/(allFibres_P(k).Area);
end
%%
for k = 1:numEdges
    curvature_P(k) = (allFibres_P(k).MajorAxisLength)/(allFibres_P(k).Perimeter);
end
%%
for k = 1:numEdges
    curvature_P2(k) = (allFibres_P(k).MajorAxisLength)/(0.5*(allFibres_P(k).Perimeter));
end
%%
for k = 1:numEdges
    AspectRatio(k) = (allFibres_P(k).MinorAxisLength)/(allFibres_P(k).MajorAxisLength);
end


%%


dataOut                     = dataIn;

if levs==1
    dataOut(:,:,3)          = dataIn;
end
dataOut(:,:,3)              = dataOut(:,:,3).*uint8(imerode(1-allFibres,ones(3)));
dataOut(:,:,2)              = dataOut(:,:,3).*uint8(imerode(1-allFibres,ones(3)));




% different metrics extracted from the edges, more structured cells will
% have longer edges, more oriented, closer together
fibronectinOut.avOrientation        = mean([allFibres_P.Orientation]);
fibronectinOut.stdOrientation       = std([allFibres_P.Orientation]);
fibronectinOut.avMajorAxisLength    = mean([allFibres_P.MajorAxisLength]);
fibronectinOut.avMinorAxisLength    = mean([allFibres_P.MinorAxisLength]);
fibronectinOut.stdMajorAxisLength   = std([allFibres_P.MajorAxisLength]);
fibronectinOut.stdMinorAxisLength   = std([allFibres_P.MinorAxisLength]);
fibronectinOut.gapArea              = sum(distMapEdges(:)>distGap);

fibronectinOut.gapAreaRel           = fibronectinOut.gapArea/rows/cols ;
fibronectinOut.avMajorAxisRel       = fibronectinOut.avMajorAxisLength/rows;
fibronectinOut.avMinorAxisRel       = fibronectinOut.avMinorAxisLength/rows;
fibronectinOut.MinAxMajAx           = fibronectinOut.avMinorAxisLength /fibronectinOut.avMajorAxisLength ;



fibronectinOut.distgaps             = double(mean(distMapEdges(distMapEdges>0)));

fibronectinOut.numEdges             = numEdges;
% fibronectinOut.avOrientation2        = mean([Fibronect_5.Orientation]);
% fibronectinOut.stdOrientation2       = std([Fibronect_5.Orientation]);

%fibronectinOut.stdOrientation3       = mean(abs([allFibres_P.Orientation]-mean([allFibres_P.Orientation])));
% fibronectinOut.stdOrientation4       = mean(abs([Fibronect_5.Orientation]-mean([Fibronect_5.Orientation])));

% Since orientations of 89 and -89 are very close and produce incorrect
% orientations and standard deviations, calculate via histograms

[yOrient,xOrient]                   = hist([allFibres_P.Orientation],[-90:5:90]);
[maxO,maxO_l]                       = max(yOrient);
yOrient2                            = circshift(yOrient,-maxO_l+18);

fibronectinOut.avOrientation5        = xOrient(maxO_l);
fibronectinOut.stdOrientation5       = 5*sum(yOrient2>(0.7071*maxO));


fibronectinOut.largestGap            = double(max(distMapEdges(:)));

% add intensity metrics 
fibronectinOut.meanIntensity         = mean2(dataIn(:,:,maxChan));
fibronectinOut.stdIntensity          = std2(dataIn(:,:,maxChan));

% find how many regions are above an intensity threshold

brightRegions                       = (dataIn(:,:,maxChan)>(fibronectinOut.meanIntensity+1*fibronectinOut.stdIntensity));
[brightRegions2, numRegions]        = bwlabel(bwmorph(brightRegions,'majority'));


fibronectinOut.numRegions           = numRegions;
fibronectinOut.areaFibres           = sum(allFibres_L(:)>0)/rows/cols; 


fibronectinOut2.edges                = allFibres_L;
fibronectinOut2.dist                 = distMapEdges;
fibronectinOut2.regions              = brightRegions2;

%% Calculate the density of traces and the orientation

mapOrientation1                     = zeros(rows,cols);
mapDensity                          = zeros(rows,cols);
mapOrientation2                     = zeros(rows,cols);
for k=1: numEdges
    sizeSq                          = round(15*allFibres_P(k).EquivDiameter);
    rr                              = max(1,-sizeSq+round(allFibres_P(k).Centroid(2))):min(rows,+sizeSq+round(allFibres_P(k).Centroid(2)) );
    cc                              = max(1,-sizeSq+round(allFibres_P(k).Centroid(1))):min(rows,+sizeSq+round(allFibres_P(k).Centroid(1)) );
    mapOrientation1(rr,cc)          = mapOrientation1(rr,cc) + allFibres_P(k).Orientation;
    mapDensity(rr,cc)               = mapDensity(rr,cc) + 1;    
end

mapOrientation                      = mapOrientation1./(1+mapDensity);

for k=1: numEdges
    sizeSq                          = round(15*allFibres_P(k).EquivDiameter);
    rr                              = max(1,-sizeSq+round(allFibres_P(k).Centroid(2))):min(rows,+sizeSq+round(allFibres_P(k).Centroid(2)) );
    cc                              = max(1,-sizeSq+round(allFibres_P(k).Centroid(1))):min(rows,+sizeSq+round(allFibres_P(k).Centroid(1)) );
    mapOrientation2(rr,cc)          = mapOrientation2(rr,cc) + (allFibres_P(k).Orientation- mapOrientation(rr,cc)).^2;
end


fibronectinOut2.mapDensity          = mapDensity;
fibronectinOut2.mapOrientation      = mapOrientation;
fibronectinOut2.mapVariability      = sqrt(mapOrientation2./(1+mapDensity));
%%
fibronectinOut.avgCircularity       = mean([allFibres_P.Circularity]);
fibronectinOut.stdCircularity       = std([allFibres_P.Circularity]);
fibronectinOut.avgCurvature         = mean(curvature);
fibronectinOut.stdCurvature         = std(curvature);
fibronectinOut.avgEccentricity      = mean([allFibres_P.Eccentricity]);
fibronectinOut.stdEccentricity      = std([allFibres_P.Eccentricity]);
fibronectinOut.avgCurvature_P         = mean(curvature_P);
fibronectinOut.stdCurvature_P         = std(curvature_P);
fibronectinOut.avgCurvature_P2         = mean(curvature_P2);
fibronectinOut.stdCurvature_P2         = std(curvature_P2);
fibronectinOut.avgAspectRatio       = mean(AspectRatio);
fibronectinOut.stdAspectRatio       = std(AspectRatio);