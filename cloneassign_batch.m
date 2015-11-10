clc; close all hidden; fclose('all'); clear all; clear java;
pause(2);
 
%% =====DESCRIPTION=====

% Clonal assignment of multiple polyclonal populations

% ==Usage: 
% User specifies variables in "USER INPUT" section.
% User specifies location of tfmTHETA-PHI grid property files
% (in section "%% Transformed THETA-PHI Grid Info"

% ==Output files: "cloneassign output/Fsk*/Result/*CloneCellCount.tif'"
% Save cell ct assigned to each clone for each polyclonal population

%==Subfunctions:
% cloneassign.m
% clonemaskreg.m


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
RelCellBrightQuery='y';

% RelCellBrightLoQuery = lowest relative cell brightness (xAF) to include in analysis
% RelCellBrightHiQuery = highest relative cell brightness (xAF) (-1=max)
RelCellBrightLoQuery=20;
RelCellBrightHiQuery=-1;

% Directory, color data to be clonally assigned populations
AssignpopFolder='ClnColorDEMO/RCB output/';

% Directory, chromatic spread data (with .../Isofrac*/CloneMasks)
cmodspdFolder='ClnColorDEMO/cmodspd output/';

% Directory, output 
BatchRbwInCloneProcessFolder='ClnColorDEMO/cloneassign output/';

% Clone mask img output size (pixels)
TargetImgCorrTHETAPxWidth=3240;     % Width (THETA)         
TargetImgCorrPHIPxHeight=810;           % Height (PHI)   

% 'y' to use create registered clonal masks (to be stored in "Working_directory/Fsk*/Isofrac*/REG)
RegShiftMaskQuery='y';

% Directory, bUnwarpJ data (with .../Isofrac*/CloneMasks)
RegShiftSpecFolder='ClnColorDEMO/bUnwarpJ output/';

% Registered clonal mask img size (in pixels)
RegImg_Width=TargetImgCorrTHETAPxWidth;      % Width (THETA)      
RegImg_Height=TargetImgCorrPHIPxHeight;       % Height (PHI)  

% # Header rows, columns in bUnwarpJ "*inverse_transf RAW.txt" file
% file ("Raw" format) from bUnwarpJ. Values used to properly load the matrix.
NumXHeaderRow=4;
NumYHeaderRow=2;
NumHeaderCol=1;

% Directory, registered clone mask 
RegShiftMaskFolderName='/REG/';

% Partial Tranformed THETA-PHI Grid range
TargetImgCorrTHETASpan=32400*2;     % THETA range
TargetImgCorrPHISpan=8100*2;           % PHI range
 
% 'y' to remove overlapped regions between all clonal masks per isoval
RemoveClnMaskOverlapQuery='y';

% 'y' if want to use Clonal mask color file to color inclonemask cells.
% Otherwise will color using MATLAB colormap
ClonalMaskColorSpecQuery='y';

% Clone color code file
ClonalMaskColorSpec_FileNameString='ClnColorAnalysis/mcln color code.txt';


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


%% Create putput folder structure for DEMO

if ~exist(BatchRbwInCloneProcessFolder, 'dir')
    mkdir(BatchRbwInCloneProcessFolder);
end;

AssignpopFileList=dir(fullfile(char(AssignpopFolder),'pcln0*RelCellBright.txt'));
NumAssignpopFiles=length(AssignpopFileList);

if RegShiftMaskQuery=='y'
    RegFiles=dir(fullfile(char(RegShiftSpecFolder),'*inverse_transf RAW.txt'));
    RegFileNameString=RegFiles(1).name;
end;

for i=1:NumAssignpopFiles
    mkdir(strcat(BatchRbwInCloneProcessFolder,'Fsk',sprintf('%02d',i)));
    AssignpopFileNameString=AssignpopFileList(i).name;
    copyfile(strcat(AssignpopFolder,AssignpopFileNameString),strcat(BatchRbwInCloneProcessFolder,'Fsk',sprintf('%02d',i)));
    if RegShiftMaskQuery=='y'
        RegFiles=dir(fullfile(char(RegShiftSpecFolder),'*inverse_transf RAW.txt'));
        copyfile(strcat(RegShiftSpecFolder,RegFileNameString),strcat(BatchRbwInCloneProcessFolder,'Fsk',sprintf('%02d',i)));
    end;
end;    

cspdIsofracDirectory=dir(fullfile(cmodspdFolder,'Isofrac*'));
cspdIsofrac_SubFolderIndex=[cspdIsofracDirectory(:).isdir];
cspdIsofrac_SubFolderName={cspdIsofracDirectory(cspdIsofrac_SubFolderIndex).name}';
cspdIsofrac_SubFolderName(ismember(cspdIsofrac_SubFolderName,{'.','..'}))=[];
NumcspdIsofracDir=length(cspdIsofrac_SubFolderName);

for i=1:NumAssignpopFiles
    for j=1:NumcspdIsofracDir
        mkdir(strcat(BatchRbwInCloneProcessFolder,'Fsk',sprintf('%02d',i),'/',cspdIsofrac_SubFolderName{j}));
        copyfile(strcat(cmodspdFolder,cspdIsofrac_SubFolderName{j},'/CloneMasks/mcln0*CloneMask.tif'),strcat(BatchRbwInCloneProcessFolder,'Fsk',sprintf('%02d',i),'/',cspdIsofrac_SubFolderName{j}));
    end;
end;



%% Load clonal masks

% 1st level sub-directory: one for each population to be assigned. Registry file is stored @ this level.

BatchFlaskSampleFolderDirectory=dir(fullfile(BatchRbwInCloneProcessFolder,'Fsk*'));
BatchFlaskSampleFolder_SubFolderIndex=[BatchFlaskSampleFolderDirectory(:).isdir];
BatchFlaskSampleFolder_SubFolderName={BatchFlaskSampleFolderDirectory(BatchFlaskSampleFolder_SubFolderIndex).name}';
BatchFlaskSampleFolder_SubFolderName(ismember(BatchFlaskSampleFolder_SubFolderName,{'.','..'})) = [];

BatchFlaskSampleFolder_NumSubFolder=size(BatchFlaskSampleFolder_SubFolderName,1);

XRegGridCoordinate=[];
YRegGridCoordinate=[];

for jj=1:BatchFlaskSampleFolder_NumSubFolder
    
        % Color data files  
        RbwFileList=dir(fullfile(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj))),'*RelCellBright.txt'));
        NumRbwFiles=length(RbwFileList);
        
        % Registry specs
        if RegShiftMaskQuery=='y' && isempty(XRegGridCoordinate)            
                XRegGridCoordinate=dlmread(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName{jj},'/',RegFileNameString),'',[NumXHeaderRow NumHeaderCol-1  RegImg_Height+NumXHeaderRow-1 RegImg_Width+NumHeaderCol-2]);
                YRegGridCoordinate=dlmread(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName{jj},'/',RegFileNameString),'',[RegImg_Height+NumXHeaderRow+NumYHeaderRow NumHeaderCol-1  RegImg_Height*2+NumXHeaderRow+NumYHeaderRow-1 RegImg_Width+NumHeaderCol-2]);
        end;
        
        % 2nd level folder: Isofrac
        BatchClonalMaskInputFolderDirectory=dir(fullfile(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj))),'Isofrac*'));
        BatchClonalMaskInputFolder_SubFolderIndex=[BatchClonalMaskInputFolderDirectory(:).isdir];
        BatchClonalMaskInputFolder_SubFolderName={BatchClonalMaskInputFolderDirectory(BatchClonalMaskInputFolder_SubFolderIndex).name}';
        BatchClonalMaskInputFolder_SubFolderName(ismember(BatchClonalMaskInputFolder_SubFolderName,{'.','..'})) = [];
        
        BatchClonalMaskInputFolder_NumSubFolder=size(BatchClonalMaskInputFolder_SubFolderName,1);
        BatchClonalMaskInputFiles=cell(BatchClonalMaskInputFolder_NumSubFolder,1);
        NumClonalMaskInputFiles=zeros(BatchClonalMaskInputFolder_NumSubFolder,1);
        
        % Search clone masks in 2nd (unregistered) level sub-directory
        for i=1:BatchClonalMaskInputFolder_NumSubFolder
            BatchClonalMaskInputFiles{i}=dir(fullfile(char(strcat(BatchRbwInCloneProcessFolder,'/',BatchFlaskSampleFolder_SubFolderName(jj),'/',BatchClonalMaskInputFolder_SubFolderName(i))),'*.tif'));
            NumClonalMaskInputFiles(i)=length(BatchClonalMaskInputFiles{i});
        end;
        
        % Create registered clone masks      
        if RegShiftMaskQuery=='y' && ~isempty(XRegGridCoordinate)   
                for ii=1:BatchClonalMaskInputFolder_NumSubFolder    
                    fprintf('\nRegistering clone masks...\n');
                    ClonalMask_FileNameString=cell(NumClonalMaskInputFiles(ii),1);
                    for k=1:NumClonalMaskInputFiles(ii)
                        ClonalMask_FileNameString{k}=fullfile(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj),'/',BatchClonalMaskInputFolder_SubFolderName(ii))),BatchClonalMaskInputFiles{ii}(k).name);
                    end;
                    parfor k =1:NumClonalMaskInputFiles(ii)   
                            clonemaskreg(ClonalMask_FileNameString{k},XRegGridCoordinate,YRegGridCoordinate,RegImg_Width,RegImg_Height);
                    end;
                end;
                
        end;
              
        % Read clone masks in 2nd (unregistered) or 3rd (registered) level sub-directory
        if RegShiftMaskQuery=='y';
            BatchClonalMaskInputFolder_SubFolderName=strcat(BatchClonalMaskInputFolder_SubFolderName,repmat({RegShiftMaskFolderName},BatchClonalMaskInputFolder_NumSubFolder,1));
        else
            BatchClonalMaskInputFolder_SubFolderName=strcat(BatchClonalMaskInputFolder_SubFolderName,repmat({'/'},BatchClonalMaskInputFolder_NumSubFolder,1));
        end;
        
        for i=1:BatchClonalMaskInputFolder_NumSubFolder
            BatchClonalMaskInputFiles{i}=dir(fullfile(char(strcat(BatchRbwInCloneProcessFolder,'/',BatchFlaskSampleFolder_SubFolderName(jj),'/',BatchClonalMaskInputFolder_SubFolderName(i))),'*.tif'));
            NumClonalMaskInputFiles(i)=length(BatchClonalMaskInputFiles{i});
        end;
        
        
        fprintf('Loading Clonal Masks Information...\n');
        
        RbwInCloneMaskCellCtFile_HeaderRow={'RbwFileName'};
        
        for ii=1:BatchClonalMaskInputFolder_NumSubFolder
            
            ClonalMask_Img=zeros(TargetImgCorrPHIPxHeight,TargetImgCorrTHETAPxWidth);          

                    if RemoveClnMaskOverlapQuery=='y'
                        ClonalMask_ImgSum=zeros(TargetImgCorrPHIPxHeight,TargetImgCorrTHETAPxWidth);
                        nOverlapCloneMask=zeros(TargetImgCorrPHIPxHeight,TargetImgCorrTHETAPxWidth,BatchClonalMaskInputFolder_NumSubFolder);
                        AllOverlapCloneMask=zeros(TargetImgCorrPHIPxHeight,TargetImgCorrTHETAPxWidth,BatchClonalMaskInputFolder_NumSubFolder);                       
                            
                            for k =1:NumClonalMaskInputFiles(ii)
                                ClonalMask_FileNameString=fullfile(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj),'/',BatchClonalMaskInputFolder_SubFolderName(ii))),BatchClonalMaskInputFiles{ii}(k).name);
                                [ClonalMask_FilePath,ClonalMask_FileName,ClonalMask_FileExt]=fileparts(ClonalMask_FileNameString);
                                ClonalMask_FileName=regexprep(ClonalMask_FileName,'CloneMask','');
                                ClonalMask_Img=imread(ClonalMask_FileNameString);
                                ClonalMask_ImgSum=ClonalMask_ImgSum+ClonalMask_Img/NumClonalMaskInputFiles(ii);
                            end;
                                           
                            MaxNumCloneOverlap=max(max(ClonalMask_ImgSum))/(1/NumClonalMaskInputFiles(ii));
                                                        
                            for kk=2:MaxNumCloneOverlap
                                nOverlapCloneMask(:,:,ii)=min(ClonalMask_ImgSum>(1/NumClonalMaskInputFiles(ii)*(kk-1)),ClonalMask_ImgSum<=(1/NumClonalMaskInputFiles(ii)*(kk)));
                                AllOverlapCloneMask(:,:,ii)=AllOverlapCloneMask(:,:,ii)+ nOverlapCloneMask(:,:,ii);
                            end;

                            clearvars  ClonalMask_ImgSum;

                    end;


                    for k =1:NumClonalMaskInputFiles(ii)
                        ClonalMask_FileNameString=fullfile(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj),'/',BatchClonalMaskInputFolder_SubFolderName(ii))),BatchClonalMaskInputFiles{ii}(k).name);
                        [ClonalMask_FilePath,ClonalMask_FileName,ClonalMask_FileExt]=fileparts(ClonalMask_FileNameString);
                        ClonalMask_FileName=regexprep(ClonalMask_FileName,'CloneMask','');
                        ClonalMask_Img=imread(ClonalMask_FileNameString);
                        
                        if RemoveClnMaskOverlapQuery=='y'
                                ClonalMask_Img(find(AllOverlapCloneMask(:,:,ii)))=0;
                        end;
                        
                        ClonalMask_boundary=bwboundaries(ClonalMask_Img,'noholes');
                        ClonalMask_YX=cell2mat(ClonalMask_boundary);
                        
                        if size(ClonalMask_YX,1)>0  
                                ClonalMask_HullCorrPHI=(TargetImgCorrPHISpan-(ClonalMask_YX(:,1)*TargetImgCorrPHISpan/TargetImgCorrPHIPxHeight))-TargetImgCorrPHISpan/2;
                                ClonalMask_HullCorrTHETA=ClonalMask_YX(:,2)*TargetImgCorrTHETASpan/TargetImgCorrTHETAPxWidth-TargetImgCorrTHETASpan/2;
                        else
                                ClonalMask_HullCorrPHI=[-TargetImgCorrPHISpan/2-1000;-TargetImgCorrPHISpan/2-10000;-TargetImgCorrPHISpan/2-1000;-TargetImgCorrPHISpan/2-10000];
                                ClonalMask_HullCorrTHETA=[-TargetImgCorrTHETASpan/2-1000;-TargetImgCorrTHETASpan/2-1000;-TargetImgCorrTHETASpan/2-10000;-TargetImgCorrTHETASpan/2-10000];
                        end;

                        if ii==1
                            HeaderRow_ClonalMask_FileName=char(regexprep(ClonalMask_FileName,'Isoval(\d*).(\d*)',cell2mat(BatchClonalMaskInputFolder_SubFolderName(:)')));
                            RbwInCloneMaskCellCtFile_HeaderRow=horzcat(RbwInCloneMaskCellCtFile_HeaderRow,char(strcat('CellCt #, ',HeaderRow_ClonalMask_FileName)),char(strcat('CellCt %, ',HeaderRow_ClonalMask_FileName)));
                        end;

                        ClonalMask_HullCorrTHETAPHI=struct('FileName',{''},'CorrTHETAvertices',zeros(size(ClonalMask_HullCorrTHETA),'double'),'CorrPHIvertices',zeros(size(ClonalMask_HullCorrPHI),'double'));

                        ClonalMask_HullCorrTHETAPHI.FileName=ClonalMask_FileName;
                        ClonalMask_HullCorrTHETAPHI.CorrTHETAvertices=ClonalMask_HullCorrTHETA;
                        ClonalMask_HullCorrTHETAPHI.CorrPHIvertices=ClonalMask_HullCorrPHI;

                        Batch_ClonalMask_HullCorrTHETAPHI_CellArray{BatchClonalMaskInputFolder_NumSubFolder-ii+1,k}=ClonalMask_HullCorrTHETAPHI;

                        clearvars ClonalMask_boundary ClonalMask_YX ClonalMask_HullCorrTHETA ClonalMask_HullCorrPHI

                    end;

      end;

      fprintf('Done.\n\n\n');

      clearvars ClonalMask_Img AllOverlapCloneMask ClonalMask_HullCorrTHETAPHI;
      
      
        %% Clone Coloring

        XFP_CloneColorBGR=zeros(size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2),3);

        if ClonalMaskColorSpecQuery=='y';
            [ClonalMaskColorSpec_FilePath,ClonalMaskColorSpec_FileName,ClonalMaskColorSpec_FileExt]=fileparts(ClonalMaskColorSpec_FileNameString);
            ClonalMaskColorSpec_FileID=fopen(ClonalMaskColorSpec_FileNameString);
            XFP_CloneColorBGR=dlmread(ClonalMaskColorSpec_FileNameString,'\t',[1,1,size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2),3]);

            XFP_CloneColorBGR=bsxfun(@rdivide,XFP_CloneColorBGR,sqrt(XFP_CloneColorBGR(:,1).^2+XFP_CloneColorBGR(:,2).^2+XFP_CloneColorBGR(:,3).^2));
 
        else
            XFP_CloneColorBGR=colormap(jet(size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2)));

        end;

      
      %% Load files to clonal assign

        RbwFileList=dir(fullfile(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj))),'pcln0*RelCellBright.txt'));
        NumRbwFiles=length(RbwFileList);

        RbwInCloneMaskOutputMatrixArray=cell(NumRbwFiles,1);
        RbwInCloneMaskCellCtFile_HeaderColArray=cell(NumRbwFiles,1);
        RbwFileNameString=cell(NumRbwFiles,1);

        for k=1:NumRbwFiles
            RbwInCloneMaskOutputMatrixArray{k}=zeros(1,(NumClonalMaskInputFiles(ii)+1)*2);
            RbwInCloneMaskCellCtFile_HeaderColArray{k}='';
            RbwFileNameString{k}=fullfile(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj))),RbwFileList(k).name);
        end;

        RbwInCloneMaskCellCtFile_HeaderRow={RbwInCloneMaskCellCtFile_HeaderRow{1},RbwInCloneMaskCellCtFile_HeaderRow{2:2:end-1},char('CellCt# OutofMasks'),RbwInCloneMaskCellCtFile_HeaderRow{3:2:end},char('CellCt% OutofMasks')};

        for k =1:NumRbwFiles
            
            % Function
            [FileNameSpec,DataPtCt,XFP_InCloneMaskCellCt]=cloneassign(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,XFP_CloneColorBGR,RbwFileNameString{k},DataPtCtChoice,RchannelChoice,GchannelChoice,BchannelChoice,RegShiftMaskQuery,THETAchoice,PHIchoice,CorrTHETAbinMatrix_FileNameString,CorrTHETAbinEdge_FileNameString,CorrPHIbinMatrix_FileNameString,CorrPHIbinEdge_FileNameString,THETACTRxcorrgrid_FileNameString,PHICTRycorrgrid_FileNameString,IntensityRangeCorrOutputSize,RelCellBrightQuery,RelCellBrightLoQuery,RelCellBrightHiQuery,TargetImgCorrTHETASpan,TargetImgCorrPHISpan,TargetImgCorrTHETAPxWidth,TargetImgCorrPHIPxHeight);

            RbwInCloneMaskOutputMatrixArray{k}=[XFP_InCloneMaskCellCt',DataPtCt-sum(XFP_InCloneMaskCellCt),XFP_InCloneMaskCellCt'/DataPtCt*100,(DataPtCt-sum(XFP_InCloneMaskCellCt))/DataPtCt*100];

            RbwInCloneMaskCellCtFile_HeaderColArray{k}={FileNameSpec};

        end;

        RbwInCloneMaskCellCtFile_HeaderCol='';
        RbwInCloneMaskCellCtFile_HeaderCol=horzcat(RbwInCloneMaskCellCtFile_HeaderColArray{:});

        RbwInCloneMaskOutputMatrix=zeros(NumRbwFiles,(NumClonalMaskInputFiles(ii)+1)*2);
        RbwInCloneMaskOutputMatrix=cell2mat(RbwInCloneMaskOutputMatrixArray);


        % Save clonal assignment
        
        time=clock;

        mkdir(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj))),'Result');

        if NumRbwFiles>1
            RbwInCloneMaskCellCtFile_FileNameString=strcat(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj))),'/Result/',sprintf('%02d',time(2)),sprintf('%02d',time(3)),sprintf('%02d',time(4)),sprintf('%02d',time(5)),' BatchCloneCellCt.txt');
        elseif NumRbwFiles==1
            RbwInCloneMaskCellCtFile_FileNameString=strcat(char(strcat(BatchRbwInCloneProcessFolder,BatchFlaskSampleFolder_SubFolderName(jj))),'/Result/',char(RbwInCloneMaskCellCtFile_HeaderColArray{1}),' CloneCellCt.txt');    
        end;

        RbwInCloneMaskCellCtFile_FileID=fopen(RbwInCloneMaskCellCtFile_FileNameString,'w');

        for i=1:numel(RbwInCloneMaskCellCtFile_HeaderRow)
            fprintf(RbwInCloneMaskCellCtFile_FileID,'%s\t',RbwInCloneMaskCellCtFile_HeaderRow{i});
        end;

        fprintf(RbwInCloneMaskCellCtFile_FileID,'\n');

        for k =1:NumRbwFiles 
            fprintf(RbwInCloneMaskCellCtFile_FileID,'%s\t',RbwInCloneMaskCellCtFile_HeaderCol{k});
            fprintf(RbwInCloneMaskCellCtFile_FileID,repmat('%f\t',1,size(RbwInCloneMaskOutputMatrix(k,:),2)),RbwInCloneMaskOutputMatrix(k,:));
            fprintf(RbwInCloneMaskCellCtFile_FileID,'\n');
        end;

        fclose(RbwInCloneMaskCellCtFile_FileID);

        fclose('all');
        
end;

matlabpool close;
