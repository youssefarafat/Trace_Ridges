function [fibronectinOut,fibronectinOut2,dataOut] =  Trace_Ridges(dataIn,cannySize,distGap)


%function [allFibres,allFibres_P,dataOut] = traceRidges(dataIn,sizeFilt,cannySize)



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
    dataIn_MIP                  = rgb2gray(dataIn);
else
    [dataIn_MIP]                = max(dataIn(:,:,:),[],3);
end

%detect the background dark/bright and leave with intensity always as high
[y,x]                           = hist(double(dataIn_MIP(:)));

if (sum(y(6:10))/sum(y(1:5)))>2
    dataIn_MIP                  = max(max(dataIn_MIP))-dataIn_MIP;
end

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
[allEdges]                  = edge(dataIn_MIP,'canny',[],cannySize);
% remove edges from the watersheds
allRidgesNoEdges            = (allWatersheds).*(allEdges==0);
% remove all watersheds that are far from edges (i.e. basins)

distMapEdges                = bwdist(allEdges);


allRidgesNoEdges2           = allRidgesNoEdges.*(distMapEdges<5);
% remove darker regions of the watersheds
lowLevelRidges              = 0.9*255*graythresh(dataIn_MIP);
allRidgesNoEdges3           = allRidgesNoEdges2.*double(dataIn_MIP>(lowLevelRidges));

% clean small spurious ridges that are coming from the watershed
allRidgesNoEdges4           = bwmorph(allRidgesNoEdges3,'spur',5);

% remove very small ridges
allRidgesNoEdges4_L         = bwlabel(allRidgesNoEdges4);
allRidgesNoEdges4_P         = regionprops(allRidgesNoEdges4_L,'MajoraxisLength');
%allRidgesNoEdges4_P         = regionprops(allRidgesNoEdges4_L,'EulerNumber','Area','MajoraxisLength','MinoraxisLength','Eccentricity');

% fill holes and thin


% discard anything that is smaller than 10 in major axis
largeRidges1                = ismember(allRidgesNoEdges4_L,find([allRidgesNoEdges4_P.MajorAxisLength]>10));
largeRidges1_L              = bwlabel(largeRidges1);
largeRidges1_P              = regionprops(largeRidges1_L,'EulerNumber');
%largeRidges1_P              = regionprops(largeRidges1_L,'EulerNumber','Area','MajoraxisLength','MinoraxisLength','Eccentricity');


% fill and thin This is particularly useful with very busy images with
% small holes that appear between parallel fibres. For large holes, this
% can be formed by valid ridges that close a polygon.

% this are ridges with no holes
largeRidges2_NoHoles        = (ismember(largeRidges1_L,find([largeRidges1_P.EulerNumber]==1)));

% find the holes and keep only the smaller ones
largeRidges2_WithHoles      = imfill(ismember(largeRidges1_L,find([largeRidges1_P.EulerNumber]<1)),'holes');
largeRidges2_JustHoles      = largeRidges2_WithHoles.*(1-largeRidges1);
largeRidges2_JustHoles_L    = bwlabel(largeRidges2_JustHoles);
largeRidges2_JustHoles_P    = regionprops(largeRidges2_JustHoles_L,'Area');
% This calculates a threshold beyond which the holes will not be filled and
% thinned. But it depends on sufficient number to get mean + 3 * std
if (max(largeRidges2_JustHoles_L(:))>100)
    largeHoleArea               = mean([largeRidges2_JustHoles_P.Area])+3*std([largeRidges2_JustHoles_P.Area]);
else
    largeHoleArea               = mean([largeRidges2_JustHoles_P.Area])+1*std([largeRidges2_JustHoles_P.Area]);
end
largeRidges2_LargeHoles     = ismember(largeRidges2_JustHoles_L,find([largeRidges2_JustHoles_P.Area]>largeHoleArea));



largeRidges2_thin           = bwmorph(bwmorph(largeRidges2_WithHoles-largeRidges2_LargeHoles,'thin','inf'),'spur',7);

largeRidges3                = largeRidges2_NoHoles + largeRidges2_thin;
largeRidges3_L              = bwlabel(largeRidges3);
largeRidges3_P              = regionprops(largeRidges3_L,'MajoraxisLength','MinoraxisLength');
%largeRidges3_P              = regionprops(largeRidges3_L,'EulerNumber','Area','MajoraxisLength','MinoraxisLength','Eccentricity');

% keep anything that is straight and long, break otherwise in junctions
largeRidges4A               = ismember(largeRidges3_L,find( ([largeRidges3_P.MajorAxisLength]>10) & ([largeRidges3_P.MinorAxisLength]<3)  )    );
largeRidges4B               = ismember(largeRidges3_L,find( ([largeRidges3_P.MajorAxisLength]>20) & ([largeRidges3_P.MinorAxisLength]<10)  )    );

largeRidges4C               = ismember(largeRidges3_L,find( ([largeRidges3_P.MajorAxisLength]>20) & ([largeRidges3_P.MinorAxisLength]>=10)  )    );
largeRidges4D               = bwmorph(bwmorph(largeRidges4C,'thin','inf'),'branch');
largeRidges4E               = bwmorph(largeRidges4C-largeRidges4D,'spur',1);
largeRidges4E_L             = bwlabel(largeRidges4E);
largeRidges4E_P             = regionprops(largeRidges4E_L,'MajoraxisLength');
largeRidges4F               = ismember(largeRidges4E_L,find( ([largeRidges4E_P.MajorAxisLength]>10)));


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

%imagesc(allEdgesNotFibres6+2*allFibres1+0*allEdges)
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
allFibres_P                 = regionprops(allFibres_L,'Orientation','MajorAxisLength','MinorAxisLength','Centroid','EquivDiameter','Area','MaxFeretProperties');
allFibres_dist              = bwdist(allFibres_L);

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
fibronectinOut.gapArea_allFibres    = sum(allFibres_dist(:)>distGap);

fibronectinOut.gapAreaRel           = fibronectinOut.gapArea/rows/cols ;
fibronectinOut.avMajorAxisRel       = fibronectinOut.avMajorAxisLength/rows;
fibronectinOut.avMinorAxisRel       = fibronectinOut.avMinorAxisLength/rows;
fibronectinOut.MinAxMajAx           = fibronectinOut.avMinorAxisLength /fibronectinOut.avMajorAxisLength ;
fibronectinOut.avArea               = mean([allFibres_P.Area]);
fibronectinOut.stdArea              = std([allFibres_P.Area]);
fibronectinOut.avMaxFeretD          = mean([allFibres_P.MaxFeretDiameter]);
fibronectinOut.stdMaxFeretD         = std([allFibres_P.MaxFeretDiameter]);



fibronectinOut.distgaps             = double(mean(distMapEdges(distMapEdges>0)));
fibronectinOut.distgaps_allFibres   = double(mean(allFibres_dist(allFibres_dist>0)));

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
fibronectinOut.largestGap_allFibres  = double(max(allFibres_dist(:)));

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
fibronectinOut2.allFibres_dist       = allFibres_dist; 

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


% imagesc(mapOrientation1./(1+mapDensity))    ;colorbar
% caxis([-90 90])
% cHot = hot;
% cCool = cHot(end:-1:1,end:-1:1);
% colormap([cHot(end-100:end,:);cCool(1:100,:)])





%% 
% 
% % discard small regions, divide in 3 to use different scales below
% largeRidges             = ismember(allRidgesNoEdges2_L,find([allRidgesNoEdges2_P.Area]>150));
% mediumRidges            = ismember(allRidgesNoEdges2_L,find([allRidgesNoEdges2_P.Area]>50))-largeRidges;
% smallRidges             = ismember(allRidgesNoEdges2_L,find([allRidgesNoEdges2_P.Area]>25))-largeRidges-mediumRidges;
% 
% % Clean ridges from spurious lines due to the watershed
% largeRidges2_Holes            = bwmorph(largeRidges,'spur',15);
% mediumRidges2           = bwmorph(mediumRidges,'spur',10);
% smallRidges2            = bwmorph(smallRidges,'spur',5);
% 
% % find ridges with holes
% 
% allRidges_1             = largeRidges2_Holes+mediumRidges2+smallRidges2;
% allRidges_1_L           = bwlabel(allRidges_1);
% allRidges_1_P           = regionprops(allRidges_1_L,'EulerNumber','MajoraxisLength','MinoraxisLength');
% 
% % fill and thin
% ridgesHoles_Filled      = imfill(ismember(allRidges_1_L,find([allRidges_1_P.EulerNumber]<1)),'holes');
% ridgesHoles_thin        = bwmorph(ridgesHoles_Filled,'thin','inf');
% ridgesNoHoles           = (ismember(allRidges_1_L,find([allRidges_1_P.EulerNumber]==1)));
% 
% ridgesHoles_thin_L      = bwlabel(ridgesHoles_thin);
% ridgesHoles_thin_P      = regionprops(allRidges_1_L,'area','Eccentricity','MajoraxisLength');
% ridgesNoHoles_L         = bwlabel(ridgesNoHoles);
% ridgesNoHoles_P         = regionprops(ridgesNoHoles_L,'area','Eccentricity','MajoraxisLength');
% 
% % Analyse the structure, break if there are branching points, keep straight segments,
% % spur again, remove sections that overlap
% 
% 
% allRidges_2             = ridgesHoles_thin.*(1-ridgesNoHoles)+ridgesNoHoles;
% allRidges_2_L           = bwlabel(allRidges_2);
% 
% 
% imagesc(ridgesHoles_thin*2+ridgesNoHoles)
% 
% imagesc((allRidges_1>0) +ridgesNoHoles)
% 
% 
% %bwlabel((1-allEdges).*(allWatersheds==0).*(bwdist(allEdges)<5))