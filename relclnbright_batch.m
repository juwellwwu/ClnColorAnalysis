clc; close all hidden; fclose('all'); clear all; clear java;
 
 %% =====DESCRIPTION=====

% Calculate relative clone brightness: count # cells within a relative cell brightness range

% ==Usage: 
% User specifies variables in "USER INPUT" section
% User specifies output file location for "BatchRelClnBright*.txt",

% ==Output file: "*BatchRelClnBright*.txt'"
% Save list of relative clonal brightness
% Row1: Filename
% Row2: Total # data points in file
% Row3: Total # data points with relative cell brightness between RelCellBrightLoQuery and RelCellBrightHiQuery
% Row4: Relative clonal brightness = Row2/Row1

% ==Subfunction:
% relclnbright.m


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUT

% Folder with "*RelCellBright.txt" files
BatchInputFolder='ClnColorDEMO/RCB output'

% Color data columns 
RchannelChoice=1;
GchannelChoice=2;
BchannelChoice=3;

% 'y' to calculate relative cell brightness of each cell
RelCellBrightQuery='y';

% RelCellBrightLoQuery = b* = lowest relative cell brightness (xAF) to include in analysis
% RelCellBrightHiQuery = highest relative cell brightness (xAF) (-1=max)
RelCellBrightLoQuery=20;
RelCellBrightHiQuery=-1;


%% Function

InputFiles=dir(fullfile(BatchInputFolder,'mcln*RelCellBright.txt'));
NumInputFiles=length(InputFiles);

DataPtCtMatrix=zeros(3,NumInputFiles);

for k =1:NumInputFiles
    
    FileNameString=fullfile(BatchInputFolder,InputFiles(k).name)    
    [FileNameSpec,DataPtCtAll,DataPtCtRCB,RelClnBright]=relclnbright(FileNameString,RchannelChoice,GchannelChoice,BchannelChoice,RelCellBrightQuery,RelCellBrightLoQuery,RelCellBrightHiQuery);     
    DataPtCtFileHeaderMtx{k}=FileNameSpec;     
    DataPtCtMatrix(:,k)=[DataPtCtAll;DataPtCtRCB;RelClnBright];
     
end;

time=clock;

mkdir(strcat('ClnColorDEMO/relclnbright output/'));

DataPtCtFileNameString=strcat('ClnColorDEMO/relclnbright output/BatchRelClnBright',sprintf('%02d',time(1)),sprintf('%02d',time(2)),sprintf('%02d',time(3)),sprintf('%02d',time(4)),sprintf('%02d',time(5)),'.txt');
DataPtCtFileID=fopen(DataPtCtFileNameString,'w');

for i=1:numel(DataPtCtFileHeaderMtx);
    fprintf(DataPtCtFileID,'%s\t',DataPtCtFileHeaderMtx{i});
end

fprintf(DataPtCtFileID,'\n');

dlmwrite(DataPtCtFileNameString,DataPtCtMatrix,'delimiter','\t','-append');

fclose(DataPtCtFileID);
fclose('all');

