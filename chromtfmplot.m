function [FileNameSpec,DataPtCt,OrgPHITHETAClonalArea,CorrPHITHETAClonalArea,OutofMagCorrTHETAPlotRangeFlag] = chromtfmplot(FileNameString,chromtfmplotOutputFolder,DataPtCtChoice,RchannelChoice,GchannelChoice,BchannelChoice,THETAchoice,PHIchoice,CorrTHETAbinMatrix_FileNameString,CorrTHETAbinEdge_FileNameString,CorrPHIbinMatrix_FileNameString,CorrPHIbinEdge_FileNameString,THETACTRxcorrgrid_FileNameString,PHICTRycorrgrid_FileNameString,IntensityRangeCorrOutputSize,AFTieringQuery,AFMultiplierLoQuery,AFMultiplierHiQuery,TargetImgCorrTHETASpan,TargetImgCorrPHISpan,TargetImgCorrTHETAPxWidth,TargetImgCorrPHIPxHeight)

%% =====DESCRIPTION=====

% Plot histogram of one population on partial transformed THETA-PHI grid

% ==Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output files: "*pTfmHistogram.tif'"
% Plot and save histogram on partial transformed THETA-PHI grid


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUT

% 1 to plot Fig5 cts in linear scale, 2 in log10 scale
Fig5Lin1Log2Query=2;


%% Load Color Data

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
fprintf(strcat('B set to:\t',num2str(ColorChoice(3)),'\n'));

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
    
            AFMode0_XFP_LogRGB_AFMultiplier=dlmread(FileNameString,'\t',[1,(numel(HeaderMtx{1})-3)-1,size(XFP_RGB,1),(numel(HeaderMtx{1})-3)-1]);
    
            FlagLessThanAFHullQuery='y';    
            Flag_LessThanAFConvexHull=dlmread(FileNameString,'\t',[1,(numel(HeaderMtx{1})-2)-1,size(XFP_RGB,1),(numel(HeaderMtx{1})-2)-1]);
            XFP_RGB(find(Flag_LessThanAFConvexHull),:)=[];
            AFMode0_XFP_LogRGB_AFMultiplier(find(Flag_LessThanAFConvexHull),:)=[];
        
            fprintf(strcat('\nRel Cell Brightness range of data set: \t',num2str(min(AFMode0_XFP_LogRGB_AFMultiplier)),'-',num2str(max(AFMode0_XFP_LogRGB_AFMultiplier)),' (x Autofluorescence).\n'));

            if AFMultiplierLoQuery<0
                AFMultiplierLoQuery=min(AFMode0_XFP_LogRGB_AFMultiplier);
            end;
            
            if AFMultiplierHiQuery<0
                AFMultiplierHiQuery=max(AFMode0_XFP_LogRGB_AFMultiplier);
            end;
    
            fprintf(strcat('Lowest relative cell brightness (xAF) included (=b*):\t',num2str(AFMultiplierLoQuery),'\n'));
            fprintf(strcat('Highest relative cell brightness (xAF) included:\t',num2str(AFMultiplierHiQuery),'\n'));
    
            AFMultiplier_Lo_LessThanIndex=bsxfun(@lt,AFMode0_XFP_LogRGB_AFMultiplier,AFMultiplierLoQuery);
            AFMultiplier_Hi_GreaterThanIndex=bsxfun(@gt,AFMode0_XFP_LogRGB_AFMultiplier,AFMultiplierHiQuery);
    
            AFMultiplier_RemoveIndex=find(AFMultiplier_Lo_LessThanIndex+AFMultiplier_Hi_GreaterThanIndex);
            XFP_RGB(AFMultiplier_RemoveIndex,:)=[];
            AFMode0_XFP_LogRGB_AFMultiplier(AFMultiplier_RemoveIndex,:)=[];
    
            fprintf(strcat('\n# Data points after relative cell brightness exclusion: \t',num2str(size(XFP_RGB,1)))); 

            fprintf('\n');
    
end;

fclose(FileID);
fclose('all');

DataPtCt=size(XFP_RGB,1);

clearvars AFMultiplier_Lo_LessThanIndex AFMultiplier_Hi_GreaterThanIndex AFMode0_XFP_LogRGB_AFMultiplier Flag_LessThanAFConvexHull AFMultiplier_RemoveIndex;


%% Remove data points near 0 or channel saturation (262144)

RGB_SatRemoveIndex=zeros(size(XFP_RGB,1),1);

RGB_SatRemoveIndex=bsxfun(@gt,XFP_RGB(:,1),2E5)+bsxfun(@lt,XFP_RGB(:,1),10)+bsxfun(@gt,XFP_RGB(:,2),2E5)+bsxfun(@lt,XFP_RGB(:,2),10)+bsxfun(@gt,XFP_RGB(:,3),2E5)+bsxfun(@lt,XFP_RGB(:,3),10);
XFP_RGB(find(RGB_SatRemoveIndex),:)=[];

fprintf(strcat('\n# Data points after near-0 or saturated R/G/B intensities removed: \t',num2str(size(XFP_RGB,1)),'\n\n\n')); 

DataPtCt=size(XFP_RGB,1);

clearvars RGB_SatRemoveIndex;


%% DataPtCtChoice

if DataPtCtChoice>0
    if DataPtCt>DataPtCtChoice
        XFP_RGB=XFP_RGB(randperm(DataPtCt,DataPtCtChoice),:);
        DataPtCt=DataPtCtChoice;
    end;
end;


%% Filename 

if AFTieringQuery=='y'
        % FileNameSpec=strcat(FileName,' R',num2str(ColorChoice(1)),'G',num2str(ColorChoice(2)),'B',num2str(ColorChoice(3)),' T',num2str(THETAchoice),' P',num2str(PHIchoice),' RCBLo',num2str(AFMultiplierLoQuery),'Hi',num2str(AFMultiplierHiQuery));
        FileNameSpec=strcat(FileName,' T',num2str(THETAchoice),'P',num2str(PHIchoice),' RCBLo',num2str(AFMultiplierLoQuery),'Hi',num2str(AFMultiplierHiQuery));
else
        % FileNameSpec=strcat(FileName,' R',num2str(ColorChoice(1)),'G',num2str(ColorChoice(2)),'B',num2str(ColorChoice(3)),' T',num2str(THETAchoice),'P',num2str(PHIchoice));
        FileNameSpec=strcat(FileName,' T',num2str(THETAchoice),'P',num2str(PHIchoice));
end;


%% Convert to spherical coordinates

XFP_THETArad=zeros([size(DataPtCt),1],'double');
XFP_PHIrad=zeros([size(DataPtCt),1],'double');
XFP_RADIUS=zeros([size(DataPtCt),1],'double');
XFP_THETAPHI=zeros([size(DataPtCt),2],'double');

[XFP_THETArad,XFP_PHIrad,XFP_RADIUS]=cart2sph(XFP_RGB(:,3),XFP_RGB(:,2),XFP_RGB(:,1));
XFP_THETAPHI=[180/pi().*(XFP_THETArad),180/pi().*(XFP_PHIrad)];

clear XFP_RGB XFP_THETArad XFP_PHIrad XFP_RADIUS;


%% Histogram, Original->transformed THETA-PHI grid

sphTHETAchoice=THETAchoice;
sphPHIchoice=PHIchoice;
sphTHETAgrid=0:sphTHETAchoice:90;
sphPHIgrid=0:sphPHIchoice:90;

sphctrs{1}=sphPHIgrid;
sphctrs{2}=sphTHETAgrid;

Hist3Data=zeros(size(sphPHIgrid,2),size(sphTHETAgrid,2));
Hist3Pos{1}=zeros(1, ceil(90/sphPHIchoice+1)); 
Hist3Pos{2}=zeros(1, ceil(90/sphTHETAchoice+1));

[Hist3Data,Hist3Pos]=hist3([XFP_THETAPHI(:,2),XFP_THETAPHI(:,1)],sphctrs);

OrgPHITHETAClonalArea=size(find(Hist3Data),1);

clearvars Hist3Data Hist3Pos;

OrgPHIDataPtHistInd=zeros(size(DataPtCt,1),'uint32');
OrgTHETADataPtHistInd=zeros(size(DataPtCt,1),'uint32');
CorrPHIDataPtHistInd=zeros(size(DataPtCt,1),'uint32');
CorrTHETADataPtHistInd=zeros(size(DataPtCt,1),'uint32');

[OrgPHIDataPtCt,OrgPHIDataPtHistInd]=histc(XFP_THETAPHI(:,2),sphPHIgrid-PHIchoice/2);
[OrgTHETADataPtCt,OrgTHETADataPtHistInd]=histc(XFP_THETAPHI(:,1),sphTHETAgrid-THETAchoice/2);

minOrgPHIDataPtHistInd=min(OrgPHIDataPtHistInd(:));
maxOrgPHIDataPtHistInd=max(OrgPHIDataPtHistInd(:));
minOrgTHETADataPtHistInd=min(OrgTHETADataPtHistInd(:));
maxOrgTHETADataPtHistInd=max(OrgTHETADataPtHistInd(:));

CorrPHIbinMatrix=zeros([maxOrgPHIDataPtHistInd-minOrgPHIDataPtHistInd+1,maxOrgTHETADataPtHistInd-minOrgTHETADataPtHistInd+1],'uint32');
CorrPHIbinMatrix_FileID=fopen(CorrPHIbinMatrix_FileNameString);
CorrPHIbinMatrix=dlmread(CorrPHIbinMatrix_FileNameString,'\t',[minOrgPHIDataPtHistInd-1 minOrgTHETADataPtHistInd-1 maxOrgPHIDataPtHistInd-1 maxOrgTHETADataPtHistInd-1]);
fclose('all');

parfor i=1:DataPtCt;    
    CorrPHIDataPtHistInd(i,1)=CorrPHIbinMatrix(OrgPHIDataPtHistInd(i,1)-(minOrgPHIDataPtHistInd-1),OrgTHETADataPtHistInd(i,1)-(minOrgTHETADataPtHistInd-1));
end;

clearvars CorrPHIbinMatrix

CorrTHETAbinMatrix=zeros([maxOrgPHIDataPtHistInd-minOrgPHIDataPtHistInd+1,maxOrgTHETADataPtHistInd-minOrgTHETADataPtHistInd+1],'uint32');
CorrTHETAbinMatrix_FileID=fopen(CorrTHETAbinMatrix_FileNameString);
CorrTHETAbinMatrix=dlmread(CorrTHETAbinMatrix_FileNameString,'\t',[minOrgPHIDataPtHistInd-1 minOrgTHETADataPtHistInd-1 maxOrgPHIDataPtHistInd-1 maxOrgTHETADataPtHistInd-1]);
fclose('all');

parfor i=1:DataPtCt;    
    CorrTHETADataPtHistInd(i,1)=CorrTHETAbinMatrix(OrgPHIDataPtHistInd(i,1)-(minOrgPHIDataPtHistInd-1),OrgTHETADataPtHistInd(i,1)-(minOrgTHETADataPtHistInd-1));
end;

clearvars CorrTHETAbinMatrix;
clearvars OrgPHIDataPtCt OrgPHIDataPtHistInd OrgTHETADataPtCt OrgTHETADataPtHistInd
clearvars minOrgPHIDataPtHistInd maxOrgPHIDataPtHistInd minOrgTHETADataPtHistInd maxOrgTHETADataPtHistInd

minCorrPHIDataPtHistInd=min(CorrPHIDataPtHistInd(:));
maxCorrPHIDataPtHistInd=max(CorrPHIDataPtHistInd(:));
minCorrTHETADataPtHistInd=min(CorrTHETADataPtHistInd(:));
maxCorrTHETADataPtHistInd=max(CorrTHETADataPtHistInd(:));

CorrPHITHETAsubscriptEdge{1}=((minCorrPHIDataPtHistInd-1):1:(maxCorrPHIDataPtHistInd+1));
CorrPHITHETAsubscriptEdge{2}=((minCorrTHETADataPtHistInd-1):1:(maxCorrTHETADataPtHistInd+1));

CorrPHIbinEdge=zeros([maxCorrPHIDataPtHistInd-minCorrPHIDataPtHistInd+2,1]);
CorrPHIbinEdge_FileID=fopen(CorrPHIbinEdge_FileNameString);
CorrPHIbinEdge=dlmread(CorrPHIbinEdge_FileNameString,'\t',[((minCorrPHIDataPtHistInd-1)-1) 0 ((maxCorrPHIDataPtHistInd+1)-1) 0]);
fclose('all');

CorrPHITHETAEdge{1}=CorrPHIbinEdge(1:end); 

ClnCorrPHImin=CorrPHIbinEdge(2,1);
ClnCorrPHImax=CorrPHIbinEdge(end-1,1);

clearvars CorrPHIbinEdge;

CorrTHETAbinEdge=zeros([maxCorrTHETADataPtHistInd-minCorrTHETADataPtHistInd+2,1]);
CorrTHETAbinEdge_FileID=fopen(CorrTHETAbinEdge_FileNameString);
CorrTHETAbinEdge=dlmread(CorrTHETAbinEdge_FileNameString,'\t',[((minCorrTHETADataPtHistInd-1)-1) 0 ((maxCorrTHETADataPtHistInd+1)-1) 0]);
fclose('all');

CorrPHITHETAEdge{2}=CorrTHETAbinEdge(1:end);

ClnCorrTHETAmin=CorrTHETAbinEdge(2,1);
ClnCorrTHETAmax=CorrTHETAbinEdge(end-1,1);

clearvars CorrTHETAbinEdge;

CorrPHITHETAsubscriptHist3Data=zeros([IntensityRangeCorrOutputSize(1,7),IntensityRangeCorrOutputSize(1,8)],'uint32');
CorrPHITHETAsubscriptHist3Pos{1}=zeros(1,IntensityRangeCorrOutputSize(1,8),'uint32'); 
CorrPHITHETAsubscriptHist3Pos{2}=zeros(1,IntensityRangeCorrOutputSize(1,7),'uint32');

[CorrPHITHETAsubscriptHist3Data,CorrPHITHETAsubscriptHist3Pos]=hist3([CorrPHIDataPtHistInd,CorrTHETADataPtHistInd],CorrPHITHETAsubscriptEdge);

CorrPHITHETAClonalArea=size(find(CorrPHITHETAsubscriptHist3Data),1); 

clearvars CorrPHIDataPtHistInd CorrTHETADataPtHistInd CorrPHITHETAsubscriptHist3Pos CorrPHITHETAsubscriptEdge;
clearvars minCorrPHIDataPtHistInd maxCorrPHIDataPtHistInd minCorrTHETADataPtHistInd maxCorrTHETADataPtHistInd;

CorrTHETAaxismagRangeMultiplier=(TargetImgCorrTHETASpan)/IntensityRangeCorrOutputSize(1,17);
CorrTHETAaxismagRange=IntensityRangeCorrOutputSize(1,17)*CorrTHETAaxismagRangeMultiplier;

CorrPHIaxismagRangeMultiplier=(TargetImgCorrPHISpan)/IntensityRangeCorrOutputSize(1,18);
CorrPHIaxismagRange=IntensityRangeCorrOutputSize(1,18)*CorrPHIaxismagRangeMultiplier;

OutofMagCorrTHETAPlotRangeFlag=0;  

OutofMagCorrTHETAPlotRangeFlag=gt((lt(min(CorrPHITHETAEdge{1}(:)),-CorrPHIaxismagRange)+gt(min(CorrPHITHETAEdge{1}(:)),CorrPHIaxismagRange)+lt(min(CorrPHITHETAEdge{2}(:)),-CorrTHETAaxismagRange)+gt(min(CorrPHITHETAEdge{2}(:)),CorrTHETAaxismagRange)),0);

THETAxPHIycorrgridReadSubscriptMin=round(IntensityRangeCorrOutputSize(1,1)/2-IntensityRangeCorrOutputSize(1,15)*(CorrTHETAaxismagRangeMultiplier+1));
THETAxPHIycorrgridReadSubscriptMax=round(IntensityRangeCorrOutputSize(1,1)/2+IntensityRangeCorrOutputSize(1,15)*(CorrTHETAaxismagRangeMultiplier+1));

THETACTRxcorrgrid=zeros(THETAxPHIycorrgridReadSubscriptMax-THETAxPHIycorrgridReadSubscriptMin+1,IntensityRangeCorrOutputSize(1,2));
PHICTRycorrgrid=zeros(THETAxPHIycorrgridReadSubscriptMax-THETAxPHIycorrgridReadSubscriptMin+1,IntensityRangeCorrOutputSize(1,2));

THETACTRxcorrgrid_FileID=fopen(THETACTRxcorrgrid_FileNameString);
fclose('all');

PHICTRycorrgrid_FileID=fopen(PHICTRycorrgrid_FileNameString);
fclose('all');

clearvars THETACTRxcorrgrid PHICTRycorrgrid


%% FIGURE5:  Histogram, partial transformed THETA-PHI grid 
 
scrsz = get(0,'ScreenSize');

figHandle5=figure('Position',[10 420 TargetImgCorrTHETAPxWidth TargetImgCorrPHIPxHeight],'Color','black','InvertHardCopy','off','PaperType','B0','PaperPositionMode','auto','Visible','off');

hold all
grid off

axis off
axis([-TargetImgCorrTHETASpan/2 TargetImgCorrTHETASpan/2 -TargetImgCorrPHISpan/2 TargetImgCorrPHISpan/2]);
axis fill; 

set(gca,'Position',[0 0 1 1]);

if Fig5Lin1Log2Query==1
        set(gca,'Clim',[0,max(CorrPHITHETAsubscriptHist3Data(:))]);
elseif Fig5Lin1Log2Query==2
        set(gca,'Clim',[0,log10(max(CorrPHITHETAsubscriptHist3Data(:)))]);
end;

Gray8btCMcol=0:1:255;
Gray8btCM=repmat(Gray8btCMcol',[1,3])/255;
colormap(Gray8btCM);
clearvars Gray8btCMcol Gray8btCM;

shading flat;

if Fig5Lin1Log2Query==1
        Fig5CellDensityhandle=pcolor(CorrPHITHETAEdge{2},CorrPHITHETAEdge{1},CorrPHITHETAsubscriptHist3Data);
elseif Fig5Lin1Log2Query==2
        Fig5CellDensityhandle=pcolor(CorrPHITHETAEdge{2},CorrPHITHETAEdge{1},log10(CorrPHITHETAsubscriptHist3Data));
end;

set(Fig5CellDensityhandle,'EdgeColor','none');

hold off

if Fig5Lin1Log2Query==1
        export_fig(figHandle5,strcat(chromtfmplotOutputFolder,FileNameSpec,' pTfmHistogram.tif'),'-nocrop','-a1');
elseif Fig5Lin1Log2Query==2
        export_fig(figHandle5,strcat(chromtfmplotOutputFolder,FileNameSpec,' LOG10 pTfmHistogram.tif'),'-nocrop','-a1');
end;    

fprintf('Computation Complete.\n\n\n');

fprintf('\n\n');

close all hidden;  fclose('all'); clearvars -except FileNameSpec chromtfmplotOutputFolder DataPtCt OrgPHITHETAClonalArea CorrPHITHETAClonalArea OutofMagCorrTHETAPlotRangeFlag TargetImgCorrTHETASpan TargetImgCorrPHISpan TargetImgCorrTHETAPxWidth TargetImgCorrPHIPxHeight;

