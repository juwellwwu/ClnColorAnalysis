function [FileNameSpec,DataPtCtAll,DataPtCtRCB,RelClnBright] = relclnbright(FileNameString,RchannelChoice,GchannelChoice,BchannelChoice,AFTieringQuery,AFMultiplierLoQuery,AFMultiplierHiQuery)

%% =====DESCRIPTION=====

% Calculate relative clonal brightness for one clone

% ==Usage: 
% Function call from relclnbright_batch


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% Load color data

[FilePath,FileName,FileExt]=fileparts(FileNameString);

FileID=fopen(FileNameString);

HeaderRow=fgets(FileID);
HeaderMtx=textscan(HeaderRow,'%s','delimiter','\t');

DataPtCt=0;
while (fgets(FileID) ~= -1),
  DataPtCt=DataPtCt+1;
end

DataPtCtAll=DataPtCt;

for i=1:numel(HeaderMtx{1})
    i_str=num2str(i);
    dispstring=[i_str,': ',HeaderMtx{1}{i}];
    disp(dispstring)
end

fprintf('9999: *NONE* \n\n')

ColorChoice(1)=RchannelChoice;
ColorChoice(2)=GchannelChoice;
ColorChoice(3)=BchannelChoice;

 fprintf(strcat('R set to:\t',num2str(ColorChoice(1)),'\n'));
 fprintf(strcat('G set to:\t',num2str(ColorChoice(2)),'\n'));
 fprintf(strcat('B set to:\t',num2str(ColorChoice(3)),'\n\n\n'));

XFP_RGB=zeros(DataPtCt,3,'double');

for i=1:3
    if ColorChoice(i)<9999 
        XFP_RGB(:,i)=dlmread(FileNameString,'\t',[1,ColorChoice(i)-1,DataPtCt,ColorChoice(i)-1]);
    else
        XFP_RGB(:,i)=zeros(DataPtCt,1,'double');
        HeaderMtx{1}{ColorChoice(i)}=' ';
    end
end


%% Relative cell brightness, <1xAF exclusion

FlagLessThanAFHullQuery='';

if isempty(AFTieringQuery)
        AFTieringQuery='n';
end;
    
if AFTieringQuery=='y'
    
    Mode0_XFP_RGB_AFMultiplier=dlmread(FileNameString,'\t',[1,(numel(HeaderMtx{1})-3)-1,size(XFP_RGB,1),(numel(HeaderMtx{1})-3)-1]);
       
        FlagLessThanAFHullQuery='y';
        Flag_LessThanAFConvexHull=dlmread(FileNameString,'\t',[1,(numel(HeaderMtx{1})-2)-1,size(XFP_RGB,1),(numel(HeaderMtx{1})-2)-1]);
        XFP_RGB(find(Flag_LessThanAFConvexHull),:)=[];
        Mode0_XFP_RGB_AFMultiplier(find(Flag_LessThanAFConvexHull),:)=[];
        
        fprintf(strcat('\nRel Cell Brightness range of data set: \t',num2str(min(Mode0_XFP_RGB_AFMultiplier)),'-',num2str(max(Mode0_XFP_RGB_AFMultiplier)),' (xAF).\n'));    
    
        if AFMultiplierLoQuery<0
            AFMultiplierLoQuery=min(Mode0_XFP_RGB_AFMultiplier);
        end;
 
        if AFMultiplierHiQuery<0
            AFMultiplierHiQuery=max(Mode0_XFP_RGB_AFMultiplier);
        end;
    
        fprintf(strcat('Lowest relative cell brightness (xAF) included:\t',num2str(AFMultiplierLoQuery),'\n'));
        fprintf(strcat('Highest relative cell brighness (xAF) included:\t',num2str(AFMultiplierHiQuery),'\n'));

        AFMultiplier_Lo_LessThanIndex=bsxfun(@lt,Mode0_XFP_RGB_AFMultiplier,AFMultiplierLoQuery);
        AFMultiplier_Hi_GreaterThanIndex=bsxfun(@gt,Mode0_XFP_RGB_AFMultiplier,AFMultiplierHiQuery);

        AFMultiplier_RemoveIndex=find(AFMultiplier_Lo_LessThanIndex+AFMultiplier_Hi_GreaterThanIndex);
        XFP_RGB(AFMultiplier_RemoveIndex,:)=[];
        Mode0_XFP_RGB_AFMultiplier(AFMultiplier_RemoveIndex,:)=[];
        
        fprintf(strcat('\n# Data points after relative cell brightness exclusion: \t',num2str(size(XFP_RGB,1)))); 
    
end;

fclose(FileID);

DataPtCt=size(XFP_RGB,1);
DataPtCtRCB=DataPtCt;

clearvars AFMultiplier_RemoveIndex Mode0_XFP_RGB_AFMultiplier AFMultiplier_Hi_GreaterThanIndex AFMultiplier_Lo_LessThanIndex  


%% Remove data points near 0 or channel saturation (262144)

RGB_SatRemoveIndex=zeros(size(XFP_RGB,1),1);

RGB_SatRemoveIndex=bsxfun(@gt,XFP_RGB(:,1),2E5)+bsxfun(@lt,XFP_RGB(:,1),10)+bsxfun(@gt,XFP_RGB(:,2),2E5)+bsxfun(@lt,XFP_RGB(:,2),10)+bsxfun(@gt,XFP_RGB(:,3),2E5)+bsxfun(@lt,XFP_RGB(:,3),10);
XFP_RGB(find(RGB_SatRemoveIndex),:)=[];

fprintf(strcat('\n\n# Data points after near-0 or saturated R/G/B intensities removed: \t',num2str(size(XFP_RGB,1)),'\n\n\n')); 

DataPtCt=size(XFP_RGB,1);
DataPtCtRCBSatRmv=DataPtCt;

clearvars RGB_SatRemoveIndex;


%% Relative Clone Brightness

RelClnBright=DataPtCtRCB/DataPtCtAll;


%% FILENAME Designate 

if AFTieringQuery=='y'
        FileNameSpec=strcat(FileName,' RCBLo',num2str(AFMultiplierLoQuery),'Hi',num2str(AFMultiplierHiQuery));
else
        FileNameSpec=strcat(FileName,' R',num2str(ColorChoice(1)),'G',num2str(ColorChoice(2)),'B',num2str(ColorChoice(3)));
end;

close all hidden;  fclose('all'); clearvars -except FileNameSpec DataPtCtAll DataPtCtRCB RelClnBright;
clc;


end

