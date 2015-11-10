clc; close all hidden; fclose('all'); clear all; clear java;
pause(2);
 
%% =====DESCRIPTION=====

% Create clonal masks of multiple clones for clonal assignment on partial THETA-PHI grid
% Calculate chromatic stability for single clones

% ==Usage: 
% User specifies variables in "USER INPUT" section.
% User specifies unique clone list's identifier
% ([ClnName_StartIndex,ClnName_EndIndex]=regexp...)
% User specifies location of tfmTHETA-PHI grid property files
% (in section "%% Transformed THETA-PHI Grid Info"

% ==Output files: *ChromSpdPxArea*.txt'"
% Save list of chromatic spread area (of specific %nmax) analyzed single clones 
% Stored in sub-directory "Working_directory_name CloneMasks" 

% ==Output files: "ChromStability*.txt'"
% Save chromatic stability (of specific %nmax) of unique clones in working directory
% Stored in sub-directory "Working_directory_name CloneMasks 

% ==Subfunction
% chromspdanalysis.m


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

% Directory, chromatic mode and spread properties of single %nmax value
% Do not include '/' at end
BatchCvxHullInputFolder='ClnColorDEMO/cmodspd output/Isofrac0.50';

% Grid step size: Delta_THETA, Delta_PHI (deg)
% Should match "Transformed THETA-PHI Grid Info" below
THETAchoice=0.01;
PHIchoice=0.01;

% Partial Tranformed THETA-PHI Grid range
TargetImgCorrTHETASpan=32400*2;     % THETA range
TargetImgCorrPHISpan=8100*2;           % PHI range

% Clone mask img output size (pixels)
TargetImgCorrTHETAPxWidth=3240;     % Width (THETA)         
TargetImgCorrPHIPxHeight=810;           % Height (PHI)    

% 'y' to save chromatic spread pixel area (on transformed THETA-PHI grid) for each clone, each measurement
ChromSpdPxArea_SaveQuery='y';

% 'y' to save chromatic stability data for each clone
ChromStability_SaveQuery='y';


%% Load chromatic mode and spread properties

% Clonal masks to be stored in "BatchCvxHullInputFolder Name/CloneMasks"
BatchCvxHullClonalMaskOutputFolder=strcat(BatchCvxHullInputFolder,'/CloneMasks');
mkdir(BatchCvxHullClonalMaskOutputFolder);

BatchCvxHullInputFiles=dir(fullfile(BatchCvxHullInputFolder,'*ChromModeSpdProperty.txt'));
NumCvxHullInputFiles=length(BatchCvxHullInputFiles);


%% List of unique clones

ClnNameList_Struct=struct('ClnName',{''});

for k =1:NumCvxHullInputFiles
    [ClnName_StartIndex,ClnName_EndIndex]=regexp(BatchCvxHullInputFiles(k).name,'mcln[0-9][1-9]*');
    ClnNameList_Struct(k).ClnName=BatchCvxHullInputFiles(k).name(ClnName_StartIndex:ClnName_EndIndex);
end;

ClnNameList=rot90(reshape(struct2cell(ClnNameList_Struct),1,[],1),3);
[UniqueClnNameList,~,UniqueClnID]=unique(ClnNameList);


%% Load Single Clone Convex Hull Information

fprintf('Loading Chromatic Mode & Spread Info...\n');

for k =1:NumCvxHullInputFiles
     
    CvxHull_FileNameString=fullfile(BatchCvxHullInputFolder,BatchCvxHullInputFiles(k).name)
    
    [CvxHull_FilePath,CvxHull_FileName,CvxHull_FileExt]=fileparts(CvxHull_FileNameString);
    CvxHull_FileName=regexprep(CvxHull_FileName,' ChromModeSpdProperty', '');

    CvxHull_FileID=fopen(CvxHull_FileNameString);

    CvxHull_HeaderRow=fgets(CvxHull_FileID);
    CvxHull_HeaderMtx=textscan(CvxHull_HeaderRow,'%s','delimiter','\t');

    CvxHull_RowCt=0;
    while (fgets(CvxHull_FileID) ~= -1),
      CvxHull_RowCt=CvxHull_RowCt+1;
    end

    CvxHull_VerticeFace=struct('FileName',{''},'ClnName',{''},'UniqueClnID',zeros(NumCvxHullInputFiles,1,'double'),'vertices',zeros(CvxHull_RowCt,3,'double'),'faces',zeros(CvxHull_RowCt,3,'double'),'index',zeros(CvxHull_RowCt,3,'double'),'THETAPHIvertices',zeros(CvxHull_RowCt,2,'double'),'LinBGRMode',zeros(1,3,'double'));
    
    CvxHull_VerticeFace.FileName=CvxHull_FileName;
    CvxHull_VerticeFace.ClnName=char(ClnNameList(k));
    CvxHull_VerticeFace.UniqueClnID=UniqueClnID(k);
    
    CvxHull_VerticeFace.vertices=dlmread(CvxHull_FileNameString,'\t',[1,0,CvxHull_RowCt,2]);
    CvxHull_VerticeFace.vertices(isnan(CvxHull_VerticeFace.vertices(:,1)),:)=[];

    CvxHull_VerticeFace.faces=dlmread(CvxHull_FileNameString,'\t',[1,3,CvxHull_RowCt,5]);
    CvxHull_VerticeFace.faces(isnan(CvxHull_VerticeFace.faces(:,1)),:)=[];

    CvxHull_VerticeFace.index=dlmread(CvxHull_FileNameString,'\t',[1,6,CvxHull_RowCt,8]);
    CvxHull_VerticeFace.index(isnan(CvxHull_VerticeFace.index(:,1)),:)=[];
    
    CvxHull_VerticeFace.THETAPHIvertices=dlmread(CvxHull_FileNameString,'\t',[1,9,CvxHull_RowCt,10]);
    CvxHull_VerticeFace.THETAPHIvertices(isnan(CvxHull_VerticeFace.THETAPHIvertices(:,1)),:)=[];

    CvxHull_VerticeFace.LinBGRMode=dlmread(CvxHull_FileNameString,'\t',[1,11,1,13]);
    CvxHull_VerticeFace.LinBGRMode(isnan(CvxHull_VerticeFace.LinBGRMode(:,1)),:)=[];
    
    fclose(CvxHull_FileID);
    fclose('all');

    Batch_CvxHull_VerticeFace_CellArray{k}=CvxHull_VerticeFace;
         
end;

fprintf('Done.\n\n\n');


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

ClnMaskPixelArea=zeros(NumCvxHullInputFiles,1);

[ClnMaskPixelArea,UniqueClnMaskPixelArea_OrAndMin_MaskCt]=chromspdanalysis(UniqueClnNameList,Batch_CvxHull_VerticeFace_CellArray,BatchCvxHullClonalMaskOutputFolder,THETAchoice,PHIchoice,CorrTHETAbinMatrix_FileNameString,CorrTHETAbinEdge_FileNameString,CorrPHIbinMatrix_FileNameString,CorrPHIbinEdge_FileNameString,THETACTRxcorrgrid_FileNameString,PHICTRycorrgrid_FileNameString,IntensityRangeCorrOutputSize,TargetImgCorrTHETASpan,TargetImgCorrPHISpan,TargetImgCorrTHETAPxWidth,TargetImgCorrPHIPxHeight);

% Chromatic stability definition : [(Or-Min)/#Masks]/[Min]
UniqueClnChromStability=((UniqueClnMaskPixelArea_OrAndMin_MaskCt(:,1)-UniqueClnMaskPixelArea_OrAndMin_MaskCt(:,3))./(UniqueClnMaskPixelArea_OrAndMin_MaskCt(:,4)))./UniqueClnMaskPixelArea_OrAndMin_MaskCt(:,3);


%% Save chromatic spread pixel area, chromatic stability per clone

if ChromSpdPxArea_SaveQuery=='y';

            time=clock;
            ClnMaskPixelArea_FileNameString=strcat(BatchCvxHullClonalMaskOutputFolder,'/ChromSpdPxArea',sprintf('%02d',time(1)),sprintf('%02d',time(2)),sprintf('%02d',time(3)),sprintf('%02d',time(4)),sprintf('%02d',time(5)),'.txt');

            for kk=1:NumCvxHullInputFiles
                ClnMaskPixelArea_FileHeaderRow{kk}=BatchCvxHullInputFiles(kk).name;
            end;

            fid=fopen(ClnMaskPixelArea_FileNameString,'w');

            for i=1:numel(ClnMaskPixelArea_FileHeaderRow);
                fprintf(fid,'%s\t',ClnMaskPixelArea_FileHeaderRow{i});
            end;

            fprintf(fid,'\n');

            dlmwrite(ClnMaskPixelArea_FileNameString,ClnMaskPixelArea','delimiter','\t','-append','precision',16);

            fclose(fid);

end;


if ChromStability_SaveQuery=='y';
    
            time=clock;
            UniqueClnConsistency_FileNameString=strcat(BatchCvxHullClonalMaskOutputFolder,'/ChromStability',sprintf('%02d',time(1)),sprintf('%02d',time(2)),sprintf('%02d',time(3)),sprintf('%02d',time(4)),sprintf('%02d',time(5)),'.txt');

            for kk=1:length(UniqueClnNameList)
                UniqueClnConsistency_FileHeaderRow{kk}=char(UniqueClnNameList(kk));
            end;

            fid=fopen(UniqueClnConsistency_FileNameString,'w');

            for i=1:numel(UniqueClnConsistency_FileHeaderRow);
                fprintf(fid,'%s\t',UniqueClnConsistency_FileHeaderRow{i});
            end;

            fprintf(fid,'\n');

            dlmwrite(UniqueClnConsistency_FileNameString,UniqueClnMaskPixelArea_OrAndMin_MaskCt','delimiter','\t','-append');
            dlmwrite(UniqueClnConsistency_FileNameString,UniqueClnChromStability','delimiter','\t','-append');
            
            fclose(fid);

end;    
    

clc; close all hidden; fclose('all'); clear all; clear java;
matlabpool close;