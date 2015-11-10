clc; close all hidden; fclose('all'); clear all; clear java;
pause(2);

%% =====DESCRIPTION=====

% Plot spherical scatter plot and spherical histogram of populations

% ==Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output file: "BatchchromtfmplotArea*.txt'"
% Row 1) FileName
% Row 2) DataPtCt
% Row 3) OrgPHITHETAClonalArea (occupied area (px) in original THETA-PHI grid)
% Row 4) CorrPHITHETAClonalArea (occupied area (px) in transformed THETA-PHI grid)
% Row 5) OutofMagCorrTHETAPlotRangeFlag: 0 if clone area within partial transformed THETA-PHI grid, 1 otherwise 

% ==Subfunction
% chromtfmplot.m


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


 %% USER INPUT
 
% Color data folder
ColorDataInputFolder='ClnColorDEMO/RCB output/';

% Output folder
chromtfmplotOutputFolder='ClnColorDEMO/chromtfmplot output/';

% # cells to analyze (-1=all)
DataPtCtChoice=-1;

% Color data columns 
RchannelChoice=1;
GchannelChoice=2;
BchannelChoice=3;

% Grid step size: Delta_THETA, Delta_PHI (deg)
% Should match "Transformed THETA-PHI Grid Info" below
THETAchoice=0.01;
PHIchoice=0.01;

% 'y' to calculate relative cell brightness of each cell
RelCellBrightQuery='n';

% RelCellBrightLoQuery = lowest relative cell brightness (xAF) to include in analysis
% RelCellBrightHiQuery = highest relative cell brightness (xAF) (-1=max)
RelCellBrightLoQuery=20;
RelCellBrightHiQuery=-1;

% Partial Tranformed THETA-PHI Grid range
TargetImgCorrTHETASpan=32400*2;     % THETA range
TargetImgCorrPHISpan=8100*2;           % PHI range

% Clone mask img output size (pixels)
TargetImgCorrTHETAPxWidth=3240;     % Width (THETA)         
TargetImgCorrPHIPxHeight=810;           % Height (PHI)    
  

 %% For DEMO: Copy color data to "chromtfmplotOutputFolder/ColorData"

if ~exist(strcat(chromtfmplotOutputFolder,'ColorData'),'dir')
    mkdir(strcat(chromtfmplotOutputFolder,'ColorData'))
end;

ColorDataInputFileList=dir(fullfile(char(ColorDataInputFolder),'*cln*RelCellBright.txt'));
NumColorDataInputFiles=length(ColorDataInputFileList);

BatchInputFolder=strcat(chromtfmplotOutputFolder,'ColorData/');

for i=1:NumColorDataInputFiles
    ColorDataInputFileNameString=ColorDataInputFileList(i).name;
    copyfile(strcat(ColorDataInputFolder,ColorDataInputFileNameString),strcat(BatchInputFolder));
end;  

 
%% Transformed THETA-PHI Grid Info

CorrTHETAbinMatrix_FileNameString='ClnColorDEMO/T0.01P0.01Step20tfmgrid/THETA0.01PHI0.01 IntensityRangeCorrTHETAbinMatrix.txt';
CorrTHETAbinEdge_FileNameString='ClnColorDEMO/T0.01P0.01Step20tfmgrid/THETA0.01PHI0.01 IntensityRangeCorrTHETAbinEdge.txt';
CorrPHIbinMatrix_FileNameString='ClnColorDEMO/T0.01P0.01Step20tfmgrid/THETA0.01PHI0.01 IntensityRangeCorrPHIbinMatrix.txt';
CorrPHIbinEdge_FileNameString='ClnColorDEMO/T0.01P0.01Step20tfmgrid/THETA0.01PHI0.01 IntensityRangeCorrPHIbinEdge.txt';
THETACTRxcorrgrid_FileNameString='ClnColorDEMO/T0.01P0.01Step20tfmgrid/THETA0.01PHI0.01 THETACTRxcorrgrid.txt';
PHICTRycorrgrid_FileNameString='ClnColorDEMO/T0.01P0.01Step20tfmgrid/THETA0.01PHI0.01 PHICTRycorrgrid.txt';

IntensityRangeCorrOutputSize_FileNameString='ClnColorDEMO/T0.01P0.01Step20tfmgrid/THETA0.01PHI0.01 IntensityRangeCorrOutputSize.txt';
IntensityRangeCorrOutputSize_FileID=fopen( IntensityRangeCorrOutputSize_FileNameString);
IntensityRangeCorrOutputSize=dlmread(IntensityRangeCorrOutputSize_FileNameString,'\t',1,0);
fclose('all');


%% Function 

InputFiles=dir(fullfile(BatchInputFolder,'*.txt'));
NumInputFiles=length(InputFiles);

UnfmClnAreaFuncOutputMatrix=zeros(4,NumInputFiles);

for k =1:NumInputFiles
    
    FileNameString=fullfile(BatchInputFolder,InputFiles(k).name)
    
    [FileNameSpec,DataPtCt,OrgPHITHETAClonalArea,CorrPHITHETAClonalArea,OutofMagCorrTHETAPlotRangeFlag]=chromtfmplot(FileNameString,chromtfmplotOutputFolder,DataPtCtChoice,RchannelChoice,GchannelChoice,BchannelChoice,THETAchoice,PHIchoice,CorrTHETAbinMatrix_FileNameString,CorrTHETAbinEdge_FileNameString,CorrPHIbinMatrix_FileNameString,CorrPHIbinEdge_FileNameString,THETACTRxcorrgrid_FileNameString,PHICTRycorrgrid_FileNameString,IntensityRangeCorrOutputSize,RelCellBrightQuery,RelCellBrightLoQuery,RelCellBrightHiQuery,TargetImgCorrTHETASpan,TargetImgCorrPHISpan,TargetImgCorrTHETAPxWidth,TargetImgCorrPHIPxHeight);
    
     UnfmClnAreaFuncOutput_FileHeaderMtx{k}=FileNameSpec;      
     UnfmClnAreaFuncOutputMatrix(1,k)=DataPtCt;
     UnfmClnAreaFuncOutputMatrix(2,k)=OrgPHITHETAClonalArea;
     UnfmClnAreaFuncOutputMatrix(3,k)=CorrPHITHETAClonalArea;
     UnfmClnAreaFuncOutputMatrix(4,k)=OutofMagCorrTHETAPlotRangeFlag;
     
end;

time=clock;
 
UnfmClnAreaFuncOutput_FileNameString=strcat(chromtfmplotOutputFolder,'BatchchromtfmplotArea',sprintf('%02d',time(2)),sprintf('%02d',time(3)),sprintf('%02d',time(4)),sprintf('%02d',time(5)),'.txt');
UnfmClnAreaFuncOutput_FileID=fopen(UnfmClnAreaFuncOutput_FileNameString,'w');

for i=1:numel(UnfmClnAreaFuncOutput_FileHeaderMtx);
    fprintf(UnfmClnAreaFuncOutput_FileID,'%s\t',UnfmClnAreaFuncOutput_FileHeaderMtx{i});
end

fprintf(UnfmClnAreaFuncOutput_FileID,'\n');

dlmwrite(UnfmClnAreaFuncOutput_FileNameString,UnfmClnAreaFuncOutputMatrix,'delimiter','\t','-append');

fclose(UnfmClnAreaFuncOutput_FileID);

fclose('all');



