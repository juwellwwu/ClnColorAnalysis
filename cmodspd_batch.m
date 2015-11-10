clc; close all hidden; fclose('all'); clear all; clear java;
pause(2);

%% =====DESCRIPTION=====

% Calculate, save clonal chromatic mode, spread of multiple single clones

% ==Usage: 
% User specifies variables in "USER INPUT" section.
% User specifies output file location for "BatchChromMode*.txt",

% ==Output file: "*BatchChromMode*.txt'"
% Save list of chromatic modes of analyzed single clones in RGB Cartesian coordinates and spherical coordinates

% ==Subfunction
% cmodspd.m


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
BatchInputFolder='ClnColorDEMO/RCB output'

% Color data columns 
RchannelChoice=1;
GchannelChoice=2;
BchannelChoice=3;

% Chromatic spread's isovalue, in fraction of max cell count @ chromatic mode 
IsovalFrac=0.5;

% 'y' to calculate relative cell brightness of each cell
RelCellBrightQuery='y';

% RelCellBrightLoQuery = b* = lowest relative cell brightness (xAF) to include in analysis
% RelCellBrightHiQuery = highest relative cell brightness (xAF) (-1=max)
RelCellBrightLoQuery=20;
RelCellBrightHiQuery=-1;


%% Prepare File for Chromatic Mode data Output

time=clock;

mkdir(strcat('ClnColorDEMO/cmodspd output/Isofrac',sprintf('%.2f',IsovalFrac)));

ModeData_FileNameString=strcat('ClnColorDEMO/cmodspd output/Isofrac',sprintf('%.2f',IsovalFrac),'/BatchChromMode',sprintf('%02d',time(1)),sprintf('%02d',time(2)),sprintf('%02d',time(3)),sprintf('%02d',time(4)),sprintf('%02d',time(5)),'.txt');

ModeDataFileID=fopen(ModeData_FileNameString,'w');

ModeData_HeaderRow={strcat('FileName');strcat('Chromatic Mode R');strcat('Chromatic Mode G');strcat('Chromatic Mode B');strcat('Chromatic Mode THETA');strcat('Chromatic Mode PHI');strcat('Chromatic Mode RADIUS')};

for i=1:numel(ModeData_HeaderRow);
        fprintf(ModeDataFileID,'%s\t',ModeData_HeaderRow{i});
end;


%% Run Function

InputFiles=dir(fullfile(BatchInputFolder,'mcln*RelCellBright.txt'));
NumInputFiles=length(InputFiles);

for k =1:NumInputFiles
     
    FileNameString=fullfile(BatchInputFolder,InputFiles(k).name)
    
    [Rmode,Gmode,Bmode,Mode_THETA,Mode_PHI,Mode_RADIUS]=cmodspd(FileNameString,RchannelChoice,GchannelChoice,BchannelChoice,IsovalFrac,RelCellBrightQuery,RelCellBrightLoQuery,RelCellBrightHiQuery);
    
    fprintf(ModeDataFileID,'\n\r');
    fprintf(ModeDataFileID,'%s\t',ModeDataFileID,InputFiles(k).name);
    fprintf(ModeDataFileID,'%f\t',Rmode);
    fprintf(ModeDataFileID,'%f\t',Gmode);
    fprintf(ModeDataFileID,'%f\t',Bmode);
    fprintf(ModeDataFileID,'%f\t',Mode_THETA);
    fprintf(ModeDataFileID,'%f\t',Mode_PHI);
    fprintf(ModeDataFileID,'%f\t',Mode_RADIUS);
       
end;

fclose(ModeDataFileID);
fclose('all');

matlabpool close;