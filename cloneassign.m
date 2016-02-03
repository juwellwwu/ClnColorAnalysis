function [FileNameSpec,DataPtCt,XFP_InCloneMaskCellCt] = cloneassign(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,XFP_CloneColorBGR,FileNameString,DataPtCtChoice,RchannelChoice,GchannelChoice,BchannelChoice,RegShiftMaskQuery,THETAchoice,PHIchoice,CorrTHETAbinMatrix_FileNameString,CorrTHETAbinEdge_FileNameString,CorrPHIbinMatrix_FileNameString,CorrPHIbinEdge_FileNameString,THETACTRxcorrgrid_FileNameString,PHICTRycorrgrid_FileNameString,IntensityRangeCorrOutputSize,AFTieringQuery,AFMultiplierLoQuery,AFMultiplierHiQuery,TargetImgCorrTHETASpan,TargetImgCorrPHISpan,TargetImgCorrTHETAPxWidth,TargetImgCorrPHIPxHeight)

 %% =====DESCRIPTION=====

% Clonal assignment of one population

% ==Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output files: *Clone ID.txt*
% Save list of [R,G,B] data point and their assigned clonal ID

% ==Output files: *CloneCellCt.txt*
% Save list # cells assigned to each clonal ID

% ==Output files: "*ClonalAssignPlot.tif'"
% Plot data pts, colored according their assigned clonal ID
% Clone color specified by clone color code file 


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUT

% 'y' to plot clonal assignment 
CloneAssignPlotQuery='y';
DataPtPlotFig=2.5E4; % Approx # data pts to plot, -1 to plot all


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
    
        AFMode0_XFP_LogRGB_AFMultiplier=dlmread(FileNameString,'\t',[1,(numel(HeaderMtx{1})-3)-1,size(XFP_RGB,1),(numel(HeaderMtx{1})-3)-1]);
    
        FlagLessThanAFHullQuery='y';
        Flag_LessThanAFConvexHull=dlmread(FileNameString,'\t',[1,(numel(HeaderMtx{1})-2)-1,size(XFP_RGB,1),(numel(HeaderMtx{1})-2)-1]);
        XFP_RGB(find(Flag_LessThanAFConvexHull),:)=[];
        AFMode0_XFP_LogRGB_AFMultiplier(find(Flag_LessThanAFConvexHull),:)=[];
    
        fprintf(strcat('\nRel Cell Brightness range of data set: \t',num2str(min(AFMode0_XFP_LogRGB_AFMultiplier)),'-',num2str(max(AFMode0_XFP_LogRGB_AFMultiplier)),' (xAF).\n'));
    
        if AFMultiplierLoQuery<0
            AFMultiplierLoQuery=min(AFMode0_XFP_LogRGB_AFMultiplier);
        end;

        if AFMultiplierHiQuery<0
            AFMultiplierHiQuery=max(AFMode0_XFP_LogRGB_AFMultiplier);
        end;

        fprintf(strcat('Lowest relative cell brightness (xAF) included:\t',num2str(AFMultiplierLoQuery),'\n'));
        fprintf(strcat('Highest relative cell brighness (xAF) included:\t',num2str(AFMultiplierHiQuery),'\n'));

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


%% XFP Data: Cartesian to spherical coordinates. 
 
XFP_THETArad=zeros([size(DataPtCt),1],'double');
XFP_PHIrad=zeros([size(DataPtCt),1],'double');
XFP_RADIUS=zeros([size(DataPtCt),1],'double');
XFP_THETAPHI=zeros([size(DataPtCt),2],'double');

[XFP_THETArad,XFP_PHIrad,XFP_RADIUS]=cart2sph(XFP_RGB(:,3),XFP_RGB(:,2),XFP_RGB(:,1));
XFP_THETAPHI=[180/pi().*(XFP_THETArad),180/pi().*(XFP_PHIrad)];

clearvars XFP_THETArad XFP_PHIrad XFP_RADIUS;


%% FileName 

if AFTieringQuery=='y'
    FileNameSpec=strcat(FileName,' RCBLo',num2str(AFMultiplierLoQuery),'Hi',num2str(AFMultiplierHiQuery));
else
    FileNameSpec=strcat(FileName);
end;

if RegShiftMaskQuery=='y'
    FileNameSpec=strcat(FileNameSpec,' REG');
end;


%% Convert to transformed THETA-PHI coordinates

sphTHETAchoice=THETAchoice;
sphPHIchoice=PHIchoice;
sphTHETAgrid=0:sphTHETAchoice:90;
sphPHIgrid=0:sphPHIchoice:90;

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
CorrPHIbinMatrix=dlmread(CorrPHIbinMatrix_FileNameString,'\t',[minOrgPHIDataPtHistInd minOrgTHETADataPtHistInd maxOrgPHIDataPtHistInd maxOrgTHETADataPtHistInd]);
fclose('all');

CorrPHIDataPtHistInd=CorrPHIbinMatrix(sub2ind(size(CorrPHIbinMatrix),OrgPHIDataPtHistInd(:,1)-(minOrgPHIDataPtHistInd-1),OrgTHETADataPtHistInd(:,1)-(minOrgTHETADataPtHistInd-1)));

clearvars CorrPHIbinMatrix

CorrTHETAbinMatrix=zeros([maxOrgPHIDataPtHistInd-minOrgPHIDataPtHistInd+1,maxOrgTHETADataPtHistInd-minOrgTHETADataPtHistInd+1],'uint32');
CorrTHETAbinMatrix_FileID=fopen(CorrTHETAbinMatrix_FileNameString);
CorrTHETAbinMatrix=dlmread(CorrTHETAbinMatrix_FileNameString,'\t',[minOrgPHIDataPtHistInd minOrgTHETADataPtHistInd maxOrgPHIDataPtHistInd maxOrgTHETADataPtHistInd]);
fclose('all');

CorrTHETADataPtHistInd=CorrTHETAbinMatrix(sub2ind(size(CorrTHETAbinMatrix),OrgPHIDataPtHistInd(:,1)-(minOrgPHIDataPtHistInd-1),OrgTHETADataPtHistInd(:,1)-(minOrgTHETADataPtHistInd-1)));

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

CorrPHIDataPtbinEdge=CorrPHIbinEdge(CorrPHIDataPtHistInd(:,1)-(minCorrPHIDataPtHistInd-1),1);

fclose('all');

CorrPHITHETAEdge{1}=CorrPHIbinEdge(1:end); 

clearvars CorrPHIbinEdge;

CorrTHETAbinEdge=zeros([maxCorrTHETADataPtHistInd-minCorrTHETADataPtHistInd+2,1]);
CorrTHETAbinEdge_FileID=fopen(CorrTHETAbinEdge_FileNameString);
CorrTHETAbinEdge=dlmread(CorrTHETAbinEdge_FileNameString,'\t',[((minCorrTHETADataPtHistInd-1)-1) 0 ((maxCorrTHETADataPtHistInd+1)-1) 0]);

CorrTHETADataPtbinEdge=CorrTHETAbinEdge(CorrTHETADataPtHistInd(:,1)-(minCorrTHETADataPtHistInd-1),1);

fclose('all');

CorrPHITHETAEdge{2}=CorrTHETAbinEdge(1:end);

clearvars CorrTHETAbinEdge;


%% Clonal assignment

XFP_InClonalMaskQuery=zeros(DataPtCt,1);
XFP_InCloneMaskCellCt=zeros(size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2),1);

if matlabpool('size')>0
    CorrDataPtbinEdgeBlockNum=matlabpool('size');
else
    CorrDataPtbinEdgeBlockNum=1;
end;


if CorrDataPtbinEdgeBlockNum>1 && DataPtCt>5000

    CorrDataPtbinEdgeBlockRowCt=floor(size(CorrTHETADataPtbinEdge,1)/CorrDataPtbinEdgeBlockNum);

    for j=1:CorrDataPtbinEdgeBlockNum
        if j<CorrDataPtbinEdgeBlockNum
            CorrDataPtbinEdgeBlockRowStartEnd(j,:)=[CorrDataPtbinEdgeBlockRowCt*(j-1)+1,CorrDataPtbinEdgeBlockRowCt*j];
        elseif j==CorrDataPtbinEdgeBlockNum
            CorrDataPtbinEdgeBlockRowStartEnd(j,:)=[CorrDataPtbinEdgeBlockRowCt*(j-1)+1,size(CorrTHETADataPtbinEdge,1)];
        end;
    end;

    XFP_InClonalMaskQuery_temp_Array=cell(CorrDataPtbinEdgeBlockNum,1);
    XFP_InClonalMaskQuery_Array=cell(CorrDataPtbinEdgeBlockNum,1);
    XFP_InCloneMaskCellCt_Array=cell(CorrDataPtbinEdgeBlockNum,1);

    CloneMaskForLoop_Counter_Array=cell(CorrDataPtbinEdgeBlockNum,1);
    IsofracForLoop_Counter_Array=cell(CorrDataPtbinEdgeBlockNum,1);
    ForLoop_RemainCellCt=cell(CorrDataPtbinEdgeBlockNum,1);

    for j=1:CorrDataPtbinEdgeBlockNum
        XFP_InClonalMaskQuery_temp_Array{j}=horzcat([1:1:CorrDataPtbinEdgeBlockRowStartEnd(j,2)-CorrDataPtbinEdgeBlockRowStartEnd(j,1)+1]',zeros(CorrDataPtbinEdgeBlockRowStartEnd(j,2)-CorrDataPtbinEdgeBlockRowStartEnd(j,1)+1,1));
        XFP_InClonalMaskQuery_Array{j}=zeros(CorrDataPtbinEdgeBlockRowStartEnd(j,2)-CorrDataPtbinEdgeBlockRowStartEnd(j,1)+1,1);
        XFP_InCloneMaskCellCt_Array{j}=zeros(size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,1),size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2));
        CloneMaskForLoop_Counter_Array{j}=ones(1,1);
        IsofracForLoop_Counter_Array{j}=ones(1,1);
        ForLoop_RemainCellCt{j}=CorrDataPtbinEdgeBlockRowStartEnd(j,2)-CorrDataPtbinEdgeBlockRowStartEnd(j,1)+1;
    end;

    
    jj=0;
    parfor jj=1:CorrDataPtbinEdgeBlockNum

        while IsofracForLoop_Counter_Array{jj}(1,1)<size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,1)+1
 
            while CloneMaskForLoop_Counter_Array{jj}(1,1)<size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2)+1

                XFP_InClonalMaskQuery_temp_Array{jj}(1:ForLoop_RemainCellCt{jj}(1,1),2)=inpolygon(CorrTHETADataPtbinEdge(XFP_InClonalMaskQuery_temp_Array{jj}(1:ForLoop_RemainCellCt{jj}(1,1),1)+CorrDataPtbinEdgeBlockRowStartEnd(jj,1)-1,1),CorrPHIDataPtbinEdge(XFP_InClonalMaskQuery_temp_Array{jj}(1:ForLoop_RemainCellCt{jj}(1,1),1)+CorrDataPtbinEdgeBlockRowStartEnd(jj,1)-1,1),Batch_ClonalMask_HullCorrTHETAPHI_CellArray{IsofracForLoop_Counter_Array{jj}(1,1),CloneMaskForLoop_Counter_Array{jj}(1,1)}.CorrTHETAvertices,Batch_ClonalMask_HullCorrTHETAPHI_CellArray{IsofracForLoop_Counter_Array{jj}(1,1),CloneMaskForLoop_Counter_Array{jj}(1,1)}.CorrPHIvertices);
                XFP_InClonalMaskQuery_temp_Array{jj}(ForLoop_RemainCellCt{jj}(1,1)+1:end,2)=0;

                XFP_InClonalMaskQuery_Array{jj}(XFP_InClonalMaskQuery_temp_Array{jj}(find(XFP_InClonalMaskQuery_temp_Array{jj}(:,2)),1),1)=CloneMaskForLoop_Counter_Array{jj}(1,1);

                XFP_InCloneMaskCellCt_Array{jj}(IsofracForLoop_Counter_Array{jj}(1,1),CloneMaskForLoop_Counter_Array{jj}(1,1))=size(find(XFP_InClonalMaskQuery_temp_Array{jj}(:,2)),1);

                ForLoop_RemainCellCt{jj}=ForLoop_RemainCellCt{jj}-XFP_InCloneMaskCellCt_Array{jj}(IsofracForLoop_Counter_Array{jj}(1,1),CloneMaskForLoop_Counter_Array{jj}(1,1));

                XFP_InClonalMaskQuery_temp_Array{jj}(find(XFP_InClonalMaskQuery_temp_Array{jj}(:,2)),1)=-9999;
                XFP_InClonalMaskQuery_temp_Array{jj}=sortrows(XFP_InClonalMaskQuery_temp_Array{jj},[-1]);

                CloneMaskForLoop_Counter_Array{jj}(1,1)=CloneMaskForLoop_Counter_Array{jj}(1,1)+1;

            end;
            
           IsofracForLoop_Counter_Array{jj}(1,1)=IsofracForLoop_Counter_Array{jj}(1,1)+1;
           CloneMaskForLoop_Counter_Array{jj}(1,1)=1;
           
        end;
    
    end;

    clearvars XFP_InClonalMaskQuery_temp_Array ForLoop_Counter_Array ForLoop_RemainCellCt
    
    XFP_InClonalMaskQuery=cell2mat(XFP_InClonalMaskQuery_Array);
    XFP_InCloneMaskCellCt=sum(cell2mat(XFP_InCloneMaskCellCt_Array),1)';
    
    clearvars XFP_InClonalMaskQuery_Array XFP_InCloneMaskCellCt_Array XFP_InClonalMaskQuery_temp_Array ForLoop_Counter_Array ForLoop_RemainCellCt;
    pause(1);

else

    XFP_InClonalMaskQuery_temp=horzcat([1:1:DataPtCt]',zeros(DataPtCt,1));
       
    for ii=1:size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,1)
        for kk=1:size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2)

            XFP_InClonalMaskQuery_temp(:,2)=inpolygon(CorrTHETADataPtbinEdge(XFP_InClonalMaskQuery_temp(:,1),1),CorrPHIDataPtbinEdge(XFP_InClonalMaskQuery_temp(:,1),1),Batch_ClonalMask_HullCorrTHETAPHI_CellArray{ii,kk}.CorrTHETAvertices,Batch_ClonalMask_HullCorrTHETAPHI_CellArray{ii,kk}.CorrPHIvertices);
            XFP_InClonalMaskQuery(XFP_InClonalMaskQuery_temp(find(XFP_InClonalMaskQuery_temp(:,2)),1),1)=kk;

            XFP_InCloneMaskCellCt(kk,1)=XFP_InCloneMaskCellCt(kk,1)+size(find(XFP_InClonalMaskQuery_temp(:,2)),1);

            XFP_InClonalMaskQuery_temp(find(XFP_InClonalMaskQuery_temp(:,2)),:)=[];

        end;
    end;

end

XFP_InClonalMaskQuery(XFP_InClonalMaskQuery==0)=size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2)+1;
XFP_CloneColorBGR(size(Batch_ClonalMask_HullCorrTHETAPHI_CellArray,2)+1,1:3)=[0.25,0.25,0.25];

clearvars XFP_InClonalMaskQuery_temp;


%% Save clonal assignment

XFP_RGB_InCloneMask_HeaderRow={strcat('Cell R: ',HeaderMtx{1}{ColorChoice(1)});strcat('Cell G: ',HeaderMtx{1}{ColorChoice(2)});strcat('Cell B: ',HeaderMtx{1}{ColorChoice(3)});strcat('Clone ID')};

XFP_RGB_InCloneMask_FileNameString=strcat(FilePath,'/Result/',FileNameSpec,' CloneID.txt');

mkdir(FilePath,'Result');
fid=fopen(XFP_RGB_InCloneMask_FileNameString,'w');

for i=1:numel(XFP_RGB_InCloneMask_HeaderRow);
       fprintf(fid,'%s\t',XFP_RGB_InCloneMask_HeaderRow{i});
end

fprintf(fid,'\n');

dlmwrite(XFP_RGB_InCloneMask_FileNameString,horzcat(XFP_RGB,XFP_InClonalMaskQuery),'delimiter','\t','-append');

fclose(fid);
fclose('all');
    
clearvars XFP_RGB;
    

%% FIGURE1: Plot clonal assignment 

 if CloneAssignPlotQuery=='y'
 
    scrsz = get(0,'ScreenSize');

    figHandle1=figure('Position',[10 100 TargetImgCorrTHETAPxWidth TargetImgCorrPHIPxHeight],'Color','black','InvertHardCopy','off','PaperType','B0','PaperPositionMode','auto','Visible','off');
    
    hold all
    grid off
    axis off
    
    axis([-TargetImgCorrTHETASpan/2 TargetImgCorrTHETASpan/2 -TargetImgCorrPHISpan/2 TargetImgCorrPHISpan/2]);
    axis fill;

    set(gca,'Position',[0 0 1 1]);

    PrintDataPtSkip=max(1,round(DataPtCt/DataPtPlotFig));
    DataPtSize=1;
    DataPtColor=fliplr(XFP_CloneColorBGR(XFP_InClonalMaskQuery(:,1),1:3));

    Fig1ScatterPlothandle=scatter(CorrTHETADataPtbinEdge(1:PrintDataPtSkip:end,1),CorrPHIDataPtbinEdge(1:PrintDataPtSkip:end,1),DataPtSize,DataPtColor(1:PrintDataPtSkip:end,:));

    export_fig(figHandle1,strcat(FilePath,'/Result/',FileNameSpec,' ClonalAssignPlot.tif'),'-nocrop','-a1');
    
 end;

 
clearvars -except FileNameSpec DataPtCt XFP_InCloneMaskCellCt Batch_ClonalMask_HullCorrTHETAPHI_CellArray XFP_CloneColorBGR FileNameString DataPtCtChoice RchannelChoice GchannelChoice BchannelChoice THETAchoice PHIchoice CorrTHETAbinMatrix_FileNameString CorrTHETAbinEdge_FileNameString CorrPHIbinMatrix_FileNameString CorrPHIbinEdge_FileNameString THETACTRxcorrgrid_FileNameString PHICTRycorrgrid_FileNameString IntensityRangeCorrOutputSize AFTieringQuery AFMultiplierLoQuery AFMultiplierHiQuery TargetImgCorrTHETASpan TargetImgCorrPHISpan TargetImgCorrTHETAPxWidth TargetImgPHIPxWidth


