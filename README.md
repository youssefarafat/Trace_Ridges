**Introduction**
 A fully automatic and fast algorithm that traces fibre-like structures to investigate the nature of such structures.
 Trace Ridges is a method that combines methods of Watershed for ridge detection, Edge detection which breaks any ridges that run from a main ridge towards the sides of the basins, separating them from the minor ones and morphological properties to delineate fibre-like structures.

 Trace Ridges was compared against six other algoritms: Edge detection, CT Fire, Scale Space, Twombli, U-Net and Graph based, with three filtering options: no filtering, Gaussian filtering and DnCnn which is a denoising convolutional neural network. Four images of varying characterisitcs were used for the comparison, Second Harmonic Generation (SHG) images of tumour bearing mouse mammary glands, fluorescently labelled fibronectin images, Breast Cancer Biopsy images of collagen and disease mimicking extracellular matrix (Fig 1).

Trace Ridges outperformed other techniques in total and average distance error metrics and fastest time after Edge detection.

The results of the comparison can be found in this pre-print  https://biorxiv.org/cgi/content/short/2024.04.15.589548v1


 

<img width="572" alt="Screenshot 2024-04-24 at 14 25 52" src="https://github.com/youssefarafat/Trace_Ridges/assets/90700104/24a92e7d-1a9e-40ce-b481-32ebb231c451">

**Using Trace Ridges**
 
The following will act as a manual on how to use the algorithm.

Firstly ensure that the images you would like to use are in the same folder as the algorithm, alternatively you can call on the directory of images in matlab as such:

foldername1 = ('/Users/youssefarafat/Documents/GitHub/Trace_Ridges/ImagesforAnalysis')
addpath(foldername1);

or open the specific image as such:

uiopen('/Users/youssefarafat/Documents/GitHub/Trace_Ridges/ImagesforAnalysis/image1.png',1)


you can load an image as such:

dataIn = imread('/Users/youssefarafat/Documents/GitHub/Trace_Ridges/ImagesforAnalysis/image1.png')

This will open and read any .png, .tif and .jpg image.

you can view dataIn as such:

fig(1)

imagesc(dataIn)

<img width="546" alt="Screenshot 2024-04-24 at 18 43 56" src="https://github.com/youssefarafat/Trace_Ridges/assets/90700104/fb4b7206-eeda-4862-9ef0-39053c15aad3">


Trace ridges can be run in the following way by calling the function:

[fibronectinOut,fibronectinOut2,dataOut] =  Trace_Ridges(dataIn,cannySize,distGap)

dataIn is the image being fed to the algorithm. The default cannySize is 2 but that can be replaced by choosing your own value. cannySize is the size of the structural element used for Edge detection. distGap is the threshold above which an area is considered a gap. The default distGap is 13. If cannySize and distGap are left empty, the default values will stand.


A number of variables are outputed


fibronectinOut2.edges                = allFibres_L;

fibronectinOut2.dist                 = distMapEdges;

fibronectinOut2.regions              = brightRegions2;

The final trace of fibres in an image can be viewed in fibronectinOut2.edges

fig(2)

imagesc(fibronectinOut2.edges)

colormap([0 0 0 ; rand(256,3)]);

<img width="544" alt="Screenshot 2024-04-24 at 18 44 28" src="https://github.com/youssefarafat/Trace_Ridges/assets/90700104/6125913e-b435-4ad7-ae9e-8f27842563cb">


Colours are assigned randomly to discern between different fibres.

The gaps can be seen in fibronectinOut2.dist which is the distance map of the final trace.

fig(3)

imagesc(fibronectinOut2.dist)

colormap hot

<img width="544" alt="Screenshot 2024-04-24 at 18 45 33" src="https://github.com/youssefarafat/Trace_Ridges/assets/90700104/8f480109-44af-420f-8e58-80d810441b23">


Yellow regions means that the pixel is further away from a fibre.

**Metrics**
A number of metrics are extracted:

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

[yOrient,xOrient]                   = hist([allFibres_P.Orientation],[-90:5:90]);

[maxO,maxO_l]                       = max(yOrient);

yOrient2                            = circshift(yOrient,-maxO_l+18);

fibronectinOut.avOrientation5        = xOrient(maxO_l);

fibronectinOut.stdOrientation5       = 5*sum(yOrient2>(0.7071*maxO));

fibronectinOut.largestGap            = double(max(distMapEdges(:)));

fibronectinOut.meanIntensity         = mean2(dataIn(:,:,maxChan));

fibronectinOut.stdIntensity          = std2(dataIn(:,:,maxChan));

fibronectinOut.numRegions           = numRegions;

fibronectinOut.areaFibres           = sum(allFibres_L(:)>0)/rows/cols; 

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

**Multiple image processing**

If you would like to run more than one image at the same time you can specify the directory off images as follows 

baseDir1              = '/Users/amrarafat/Documents/GitHub/Trace_Ridges/ImagesforAnalysis';

dir0                = dir(strcat(baseDir1,filesep,'.png'));

numFiles            = numel(dir0);

then run a for loop 

figure

for k =1:numFiles

    h(k) = subplot(2,5,k);
    
    dataIn = imread(strcat(baseDir1,filesep,dir0(k).name))
    
    [fibronectinOut,fibronectinOut2,dataOut]=  Trace_Ridges(dataIn);
    
 results(k) = fibronectinOut.gapArea  ;
 
 imagesc(fibronectinOut2.dist)
 
 ylabel(dir0(k).name,'interpreter','none')   
 
 colormap hot  
 
 end

This will store the values (results) which can be used to run statistical comparisons if comparing two populations of images.
