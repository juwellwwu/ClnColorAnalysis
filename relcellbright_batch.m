clc; close all hidden; fclose('all'); clear all; clear java;

%% =====DESCRIPTION=====

% Determine relative cell brightness of datapts in all txt files in working directory

% ==Usage: 
% User specifies variables in "USER INPUT" section.

% ==Subfunction
% relcellbright.m


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUTS

% Parallel pool setup
if matlabpool('size')==0
    matlabpool open 4;
end;

% Folder with color data files 
BatchInputFolder='ClnColorAnalysis/ColorData/Polyclone';  

% 1xAF property file (*1xAF.txt)
AFIsosurfConvexHull_FileNameString='ClnColorDEMO/xAF output/NoFP R1G2B3 1xAF.txt';

% Folder with color data files 
OutputFolderString='ClnColorDEMO/RCB output';

% Color data columns 
RchannelChoice=1;
GchannelChoice=2;
BchannelChoice=3;


%%  Load Autofluorescence 1xAF Data

fprintf('Loading 1xAF properties for "relcellbright_batch.m" ...\n');

[AFIsosurfConvexHull_FilePath,AFIsosurfConvexHull_FileName,AFIsosurfConvexHull_FileExt]=fileparts(AFIsosurfConvexHull_FileNameString);
AFIsosurfConvexHull_FileNameTrunc= regexprep(AFIsosurfConvexHull_FileName,' 1xAF', '');

AFIsosurfConvexHull_FileID=fopen(AFIsosurfConvexHull_FileNameString);

AFIsosurfConvexHull_HeaderRow=fgets(AFIsosurfConvexHull_FileID);

AFIsosurfConvexHull_RowCt=0;
while (fgets(AFIsosurfConvexHull_FileID) ~= -1),
  AFIsosurfConvexHull_RowCt=AFIsosurfConvexHull_RowCt+1;
end


% Structure AF_IsoSurf_VerticeFace holds AF 1xAF properties

AF_IsoSurf_VerticeFaceArray=struct('FileName',{''},'vertices',zeros(AFIsosurfConvexHull_RowCt,3,'double'),'faces',zeros(AFIsosurfConvexHull_RowCt,3,'double'),'index',zeros(AFIsosurfConvexHull_RowCt,3,'double'),'LinBGRMode',zeros(1,3,'double'));

AF_IsoSurf_VerticeFaceArray.FileName=AFIsosurfConvexHull_FileNameTrunc;

AF_IsoSurf_VerticeFaceArray.LinBGRMode=10.^(dlmread(AFIsosurfConvexHull_FileNameString,'\t',[1,9,1,11]));
AF_IsoSurf_VerticeFaceArray.LinBGRMode(isnan(AF_IsoSurf_VerticeFaceArray.LinBGRMode(:,1)),:)=[];

AF_IsoSurf_VerticeFaceArray.vertices=dlmread(AFIsosurfConvexHull_FileNameString,'\t',[1,0,AFIsosurfConvexHull_RowCt,2]);
AF_IsoSurf_VerticeFaceArray.vertices(isnan(AF_IsoSurf_VerticeFaceArray.vertices(:,1)),:)=[];
AF_IsoSurf_VerticeFaceArray.vertices=10.^(AF_IsoSurf_VerticeFaceArray.vertices);
AF_IsoSurf_VerticeFaceArray.vertices=bsxfun(@minus,AF_IsoSurf_VerticeFaceArray.vertices,AF_IsoSurf_VerticeFaceArray.LinBGRMode);  

AF_IsoSurf_VerticeFaceArray.faces=dlmread(AFIsosurfConvexHull_FileNameString,'\t',[1,3,AFIsosurfConvexHull_RowCt,5]);
AF_IsoSurf_VerticeFaceArray.faces(isnan(AF_IsoSurf_VerticeFaceArray.faces(:,1)),:)=[];

AF_IsoSurf_VerticeFaceArray.index=dlmread(AFIsosurfConvexHull_FileNameString,'\t',[1,6,AFIsosurfConvexHull_RowCt,8]);
AF_IsoSurf_VerticeFaceArray.index(isnan(AF_IsoSurf_VerticeFaceArray.index(:,1)),:)=[];

fclose(AFIsosurfConvexHull_FileID);
fclose('all');

fprintf('Done.\n\n');


%% Run Function

InputFiles=dir(fullfile(BatchInputFolder,'*.txt')); 
NumInputFiles=length(InputFiles);

BatchInputFileNameStringArray=cell(NumInputFiles,1);
for k =1:NumInputFiles
    BatchInputFileNameStringArray{k}=fullfile(BatchInputFolder,InputFiles(k).name);
end;

% If 'y' for Fig2Query in relcellbright.m, use FOR. If 'n', can use parfor
for k =1:NumInputFiles   
    
relcellbright(AF_IsoSurf_VerticeFaceArray,BatchInputFileNameStringArray{k},OutputFolderString,-1,RchannelChoice,GchannelChoice,BchannelChoice);
          
end;

fclose('all');

matlabpool close;
