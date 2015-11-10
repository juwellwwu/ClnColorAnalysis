 clc; close all hidden; fclose('all'); clear all; clear java;

 %% =====DESCRIPTION=====

% Plot spherical scatter plot and spherical histogram of populations

% ==Usage: 
% User specifies variables in "USER INPUT" section.

% ==Subfunction
% chromsphplot.m


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUTS

% Parallel computing setup
if matlabpool('size')==0
    matlabpool open 4;
end;
 
% Color data folder
ColorDataInputFolder='ClnColorDEMO/cloneassign output/';

% Output folder
chromsphplotOutputFolder='ClnColorDEMO/chromsphplot output/';

% # cells to analyze (-1=all)
DataPtCtChoice=-1;

% Color data columns 
RchannelChoice=1;
GchannelChoice=2;
BchannelChoice=3;

% Clone ID column
% For pre-clonal assigned population, input CloneID channel
% For non-assigned population, or want to plot cells by their chromaticity, input 9999
CloneColorchannelChoice=4;

% Valid only when CloneColorchannelChoice =/= 9999
% 'y' to load CloneColor file, Specify "CloneColor_FileNameString"
% 'n' to color using colormap (MATLAB or user-defined). Specify # clones = ClonalColorCt
CloneColorFileQuery='y';
CloneColor_FileNameString='ClnColorAnalysis/mcln color code.txt';
ClonalColorCt=2;

% Spherical grid step size: Delta_THETA, Delta_PHI (deg)
sphTHETAchoice=0.2;
sphPHIchoice=0.2;

% 'y' to calculate relative cell brightness of each cell
RelCellBrightQuery='n';

% RelCellBrightLoQuery = lowest relative cell brightness (xAF) to include in analysis
% RelCellBrightHiQuery = highest relative cell brightness (xAF) (-1=max)
RelCellBrightLoQuery=20;
RelCellBrightHiQuery=-1;

% 'y' to plot adjusted cell count histogram. 
% Specify correction factor file.
AdjCellCtHistQuery='y';
UfmAreaCorrFactorMatrix_FileNameString='ClnColorDEMO/T0.2P0.2Step1tfmgrid/THETA0.2PHI0.2 IntensityRangeCorrFactorMatrix.txt';


%% For DEMO: Copy assigned color data to "chromsphplotOutputFolder/ColorData"

if ~exist(strcat(chromsphplotOutputFolder,'ColorData'),'dir')
    mkdir(strcat(chromsphplotOutputFolder,'ColorData'))
end;

BatchInputFolder=strcat(chromsphplotOutputFolder,'ColorData/');

% 1st level: "CloneAssignOutputFolder/Fsk*"
BatchFlaskSampleFolderDirectory=dir(fullfile(ColorDataInputFolder,'Fsk*'));
BatchFlaskSampleFolder_SubFolderIndex=[BatchFlaskSampleFolderDirectory(:).isdir];
BatchFlaskSampleFolder_SubFolderName={BatchFlaskSampleFolderDirectory(BatchFlaskSampleFolder_SubFolderIndex).name}';
BatchFlaskSampleFolder_SubFolderName(ismember(BatchFlaskSampleFolder_SubFolderName,{'.','..'})) = [];
BatchFlaskSampleFolder_NumSubFolder=size(BatchFlaskSampleFolder_SubFolderName,1);

% 2nd level folder: "CloneAssignOutputFolder/Fsk*/Result"
for jj=1:BatchFlaskSampleFolder_NumSubFolder
    
        ClnIDFileList=dir(fullfile(strcat(ColorDataInputFolder,BatchFlaskSampleFolder_SubFolderName{jj},'/Result/*CloneID.txt')));
        NumClnIDFiles=length(ClnIDFileList);
        for kk=1:NumClnIDFiles
            copyfile(strcat(ColorDataInputFolder,BatchFlaskSampleFolder_SubFolderName{jj},'/Result/',ClnIDFileList(kk).name),strcat(BatchInputFolder));
        end;

end;


%% Clone Coloring

if CloneColorchannelChoice<9999 && CloneColorFileQuery=='y';
    
    [ClonalMaskColorSpec_FilePath,ClonalMaskColorSpec_FileName,ClonalMaskColorSpec_FileExt]=fileparts(CloneColor_FileNameString);
    ClonalMaskColorSpec_FileID=fopen(CloneColor_FileNameString);
    ClonalMaskColorSpec_HeaderRow=fgets(ClonalMaskColorSpec_FileID);
    
    ClonalColorCt=0;
    
    while (fgets(ClonalMaskColorSpec_FileID) ~= -1),
       ClonalColorCt=ClonalColorCt+1;
    end;

    XFP_CloneColorBGR=zeros(ClonalColorCt,3);
    XFP_CloneColorBGR=dlmread(CloneColor_FileNameString,'\t',[1,1,ClonalColorCt,3]);
    
    XFP_CloneColorBGR=bsxfun(@rdivide,XFP_CloneColorBGR,sqrt(XFP_CloneColorBGR(:,1).^2+XFP_CloneColorBGR(:,2).^2+XFP_CloneColorBGR(:,3).^2));
    
    fclose(ClonalMaskColorSpec_FileID);
    
    % Unassigned clone color
    XFP_CloneColorBGR=vertcat(XFP_CloneColorBGR,[0.1,0.1,0.1]);
    
elseif CloneColorchannelChoice<9999 && CloneColorFileQuery=='n';
    
    XFP_CloneColorBGR=[];    
    
    % Option 1: MATLAB colormap
    XFP_CloneColorBGR=fliplr(colormap(jet(ClonalColorCt)));
    
    % Option 2: User-defined colormap
    XFP_CloneColorBGR_TargetColor=[1,0.75,0.8]; 
    XFP_CloneColorBGR_OriginColor=[0.25,0.0,0.05];    
    XFP_CloneColorBGR=zeros(ClonalColorCt,3);
    for k=1:ClonalColorCt
        XFP_CloneColorBGR(k,1)=XFP_CloneColorBGR_OriginColor(1,1)+((XFP_CloneColorBGR_TargetColor(1,1)-XFP_CloneColorBGR_OriginColor(1,1))/ClonalColorCt)*(k-1);
        XFP_CloneColorBGR(k,2)=XFP_CloneColorBGR_OriginColor(1,2)+((XFP_CloneColorBGR_TargetColor(1,2)-XFP_CloneColorBGR_OriginColor(1,2))/ClonalColorCt)*(k-1);
        XFP_CloneColorBGR(k,3)=XFP_CloneColorBGR_OriginColor(1,3)+((XFP_CloneColorBGR_TargetColor(1,3)-XFP_CloneColorBGR_OriginColor(1,3))/ClonalColorCt)*(k-1);
    end; 
    colormap(fliplr(XFP_CloneColorBGR));
    colorbar('EastOutside');
        
elseif CloneColorchannelChoice==9999
    
    XFP_CloneColorBGR=[];    

end;

clear ClonalMaskColorSpec_FileNameString ClonalMaskColorSpec_FileID ClonalMaskColorSpec_FilePath ClonalMaskColorSpec_FileExt ClonalMaskColorSpec_HeaderRow;


%% Load correction factors for adjusted cell ct histogram

if AdjCellCtHistQuery=='y';        
            UfmAreaCorrFactorMatrix_FileID=fopen(UfmAreaCorrFactorMatrix_FileNameString);
            UfmAreaCorrFactorMatrix=dlmread(UfmAreaCorrFactorMatrix_FileNameString,'\t');
            UfmAreaCorrFactorMatrix=padarray(UfmAreaCorrFactorMatrix,[1 1],'post');
else   
             UfmAreaCorrFactorMatrix=[]; 
end;


%% Function

InputFiles=dir(fullfile(BatchInputFolder,'*.txt')); 
NumInputFiles=length(InputFiles); 

SphMeshGrid=zeros(90/sphTHETAchoice+1,90/sphTHETAchoice+1,3);
[SphMeshGrid(:,:,1),SphMeshGrid(:,:,2),SphMeshGrid(:,:,3)]=sphere_1stoct(90/sphTHETAchoice);

for k =1:NumInputFiles
    
    FileNameString=fullfile(BatchInputFolder,InputFiles(k).name)   
    chromsphplot(FileNameString,chromsphplotOutputFolder,XFP_CloneColorBGR,UfmAreaCorrFactorMatrix,DataPtCtChoice,RchannelChoice,GchannelChoice,BchannelChoice,CloneColorchannelChoice,sphTHETAchoice,sphPHIchoice,SphMeshGrid,RelCellBrightQuery,RelCellBrightLoQuery,RelCellBrightHiQuery)
       
end;

fclose('all');

matlabpool close;
