 clc; close all hidden; clear all; clear java; fclose('all');

%% =====DESCRIPTION=====

% Define autofluorescence (1xAF). Plot and save properties.

% ==Usage: 
% User specifies variables in "USER INPUT" section.
% User follows prompts for further input
% User specifies output file location for "*1xAFPlot.tif", "*1xAF.txt"

% ==Output file: "*1xAFplot.tif"
% Plots 1xAF, NoFP data pts in Cartesian space

% ==Output file: "*1xAF.txt"
% Saves properties of 1xAF


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUT

AFFileNameString='ClnColorAnalysis/ColorData/NoFP.txt';

OutputFolderString='ClnColorDEMO/xAF output';


%% Load Color data

[FilePath,FileName,FileExt]=fileparts(AFFileNameString);

FileID=fopen(AFFileNameString);

HeaderRow=fgets(FileID);
HeaderMtx=textscan(HeaderRow,'%s','delimiter','\t');

DataPtCt=0;
while (fgets(FileID) ~= -1),
  DataPtCt=DataPtCt+1;
end

DataPtCtChoice=DataPtCt;

fprintf('Autofluorescence analysis. Choose color data columns: \n\n');

for i=1:numel(HeaderMtx{1})
    i_str=num2str(i);
    dispstring=[i_str,': ',HeaderMtx{1}{i}];
    disp(dispstring)
end

fprintf('9999: *NONE* \n\n')
fprintf('1:Red, 2:Green, 3:Blue.\n\n')

for i=1:3
   
    ColorChoice_i=input(horzcat(num2str(i),': '));
    
    if isempty(ColorChoice_i)
        if i==1 
            ColorChoice(i)=1; 
            fprintf(strcat('R set to default:\t',num2str(ColorChoice(i)),'\n'));
        elseif i==2 
            ColorChoice(i)=2; 
            fprintf(strcat('G set to default:\t',num2str(ColorChoice(i)),'\n'));
        elseif i==3 
            ColorChoice(i)=3; 
            fprintf(strcat('B set to default:\t',num2str(ColorChoice(i)),'\n'));
        end;
    else
        ColorChoice(i)=ColorChoice_i;
    end;
    
end;

fprintf('\n\n');

XFP_RGB=zeros(DataPtCt,3,'double');

for i=1:3
    if ColorChoice(i)<9999 
        XFP_RGB(:,i)=dlmread(AFFileNameString,'\t',[1,ColorChoice(i)-1,DataPtCt,ColorChoice(i)-1]);
    else
        XFP_RGB(:,i)=zeros(DataPtCt,1,'double');
        HeaderMtx{1}{ColorChoice(i)}=' ';
    end
end

fclose(FileID);


%% Delete data pts near 0 and channel saturation (262144)

RGB_SatRemoveIndex=zeros(size(XFP_RGB,1),1);

RGB_SatRemoveIndex=bsxfun(@gt,XFP_RGB(:,1),2E5)+bsxfun(@lt,XFP_RGB(:,1),0.01)+bsxfun(@gt,XFP_RGB(:,2),2E5)+bsxfun(@lt,XFP_RGB(:,2),0.01)+bsxfun(@gt,XFP_RGB(:,3),2E5)+bsxfun(@lt,XFP_RGB(:,3),0.01);
XFP_RGB(find(RGB_SatRemoveIndex),:)=[];

DataPtCt=size(XFP_RGB,1);

clearvars RGB_SatRemoveIndex;


%% FIGURE 1-3: Autofluorescence log plot

%== Fig1: Blue (x-axis) vs Green (y-axis)

MinB=1;
MinG=1;
MinR=1;

MaxB=1E6;
MaxG=1E6;
MaxR=1E6;

Bgrid=linspace(log10(MinB),log10(MaxB),(log10(MaxB)-log10(MinB))/0.02+1);
Ggrid=linspace(log10(MinG),log10(MaxG),(log10(MaxG)-log10(MinG))/0.02+1);
Rgrid=linspace(log10(MinR),log10(MaxR),(log10(MaxR)-log10(MinR))/0.02+1);

scrsz = get(0,'ScreenSize');

scrszWscalefactor=500/scrsz(3);
scrszHscalefactor=500/scrsz(4);

figHandle1=figure('Position',[5 scrsz(4)-(scrsz(4)*scrszHscalefactor)-100 scrsz(3)*scrszWscalefactor scrsz(4)*scrszHscalefactor],'Color','w','PaperPositionMode', 'auto');

hold on
axis on
grid on

set(gca,'XColor',[0,0,0.5],'YColor',[0,0.5,0],'Position',[0.1 0.1 0.8 0.8]);
axis([log10(MinB) log10(MaxB) log10(MinG) log10(MaxG)]);
axis square;

xlabel(strcat('B: log10 ( ',HeaderMtx{1}{ColorChoice(3)},' )'),'FontSize',10);
ylabel(strcat('G: log10 ( ',HeaderMtx{1}{ColorChoice(2)},' )'),'FontSize',10);


Hist3edge{1}=Bgrid;
Hist3edge{2}=Ggrid;

[BGHist3Data,BGHist3Pos]=hist3(horzcat(log10(XFP_RGB(:,3)),log10(XFP_RGB(:,2))),'Edges',Hist3edge);

pcolor(Bgrid,Ggrid,BGHist3Data')

set(gca,'Clim',[0,max(BGHist3Data(:))]);
colormap0white=colormap('jet');
colormap0white(1,:)=1;
colormap(colormap0white);
shading flat;

hold off


%== Fig2: Blue (x-axis) vs Red (y-axis)

scrszWscalefactor=500/scrsz(3);
scrszHscalefactor=500/scrsz(4);

figHandle2=figure('Position',[55 scrsz(4)-(scrsz(4)*scrszHscalefactor)-100 scrsz(3)*scrszWscalefactor scrsz(4)*scrszHscalefactor],'Color','w','PaperPositionMode', 'auto');

hold on
axis on
grid on

set(gca,'XColor',[0,0,0.5],'YColor',[0.5,0,0],'Position',[0.1 0.1 0.8 0.8]);
axis([log10(MinB) log10(MaxB) log10(MinR) log10(MaxR)]);
axis square;

xlabel(strcat('B: log10 ( ',HeaderMtx{1}{ColorChoice(3)},' )'),'FontSize',10);
ylabel(strcat('R: log10 ( ',HeaderMtx{1}{ColorChoice(1)},' )'),'FontSize',10);


Hist3edge{1}=Bgrid;
Hist3edge{2}=Rgrid;

[BRHist3Data,BRHist3Pos]=hist3(horzcat(log10(XFP_RGB(:,3)),log10(XFP_RGB(:,1))),'Edges',Hist3edge);

pcolor(Bgrid,Rgrid,BRHist3Data')

set(gca,'Clim',[0,max(BRHist3Data(:))]);
colormap(colormap0white);
shading flat;

hold off


%== Fig3:  Green (x-axis) vs Red (y-axis)

scrszWscalefactor=500/scrsz(3);
scrszHscalefactor=500/scrsz(4);

scrsz = get(0,'ScreenSize');
figHandle3=figure('Position',[105 scrsz(4)-(scrsz(4)*scrszHscalefactor)-100 scrsz(3)*scrszWscalefactor scrsz(4)*scrszHscalefactor],'Color','w','PaperPositionMode', 'auto');

hold on
axis on
grid on

set(gca,'XColor',[0,0.5,0],'YColor',[0.5,0,0],'Position',[0.1 0.1 0.8 0.8]);
axis([log10(MinB) log10(MaxB) log10(MinR) log10(MaxR)]);
axis square;

xlabel(strcat('G: log10 ( ',HeaderMtx{1}{ColorChoice(2)},' )'),'FontSize',10);
ylabel(strcat('R: log10 ( ',HeaderMtx{1}{ColorChoice(1)},' )'),'FontSize',10);

Hist3edge{1}=Ggrid;
Hist3edge{2}=Rgrid;

[GRHist3Data,GRHist3Pos]=hist3(horzcat(log10(XFP_RGB(:,2)),log10(XFP_RGB(:,1))),'Edges',Hist3edge);

pcolor(Ggrid,Rgrid,GRHist3Data')

set(gca,'Clim',[0,max(GRHist3Data(:))]);
colormap(colormap0white);
shading flat;
% colorbar('EastOutside');

hold off


%% GATE 1 (Auto)

GatingStartQuery=input('Gate Autofluorescence? Press "n" for no: ','s');

if GatingStartQuery=='n'
    error('Code terminated by user input. Goodbye!');
end;

FileName=strcat(FileName,' R',num2str(ColorChoice(1)),'G',num2str(ColorChoice(2)),'B',num2str(ColorChoice(3)));

if ~exist(OutputFolderString, 'dir')
            mkdir(OutputFolderString);
end;

NoiseGateTypeQuery=1;

R_LoNoiseGateCutoff=-6.0;
R_HiNoiseGateCutoff=6.0;
G_LoNoiseGateCutoff=-6.0;
G_HiNoiseGateCutoff=6.0;
B_LoNoiseGateCutoff=-6.0;
B_HiNoiseGateCutoff=6.0;

NoiseGate_RemoveRowIndex=zeros(1,size(XFP_RGB,1));

for j=1:size(XFP_RGB,1)
            
            if log10(XFP_RGB(j,1))<R_LoNoiseGateCutoff
                NoiseGate_RemoveRowIndex(1,j)=j;
            elseif log10(XFP_RGB(j,2))<G_LoNoiseGateCutoff
                NoiseGate_RemoveRowIndex(1,j)=j;
            elseif log10(XFP_RGB(j,3))<B_LoNoiseGateCutoff
                NoiseGate_RemoveRowIndex(1,j)=j;
                
            elseif log10(XFP_RGB(j,1))>R_HiNoiseGateCutoff
                NoiseGate_RemoveRowIndex(1,j)=j;
            elseif log10(XFP_RGB(j,2))>G_HiNoiseGateCutoff
                NoiseGate_RemoveRowIndex(1,j)=j;
            elseif log10(XFP_RGB(j,3))>B_HiNoiseGateCutoff
                NoiseGate_RemoveRowIndex(1,j)=j;
            end;
            
end;


%% USER-SPECIFIED GATE 2: GREEN-RED AUTOFLUORESCENCE (POLYGON) GATE

fprintf('\n\nSpecify Green-Red Autofluorescence polygon gate (1xAF).\n');

fprintf('Specify Autofluorescence Region in Fig4: Green vs Red Plot w/ polygon ROI. Double-click on polygon when done:\n\n')

Fig3Img=getframe(figHandle3,[50 50 400 400]);
figHandle4=figure('Position',[scrsz(3) 0 400 400],'OuterPosition',[scrsz(3) 0 400 400],'Resize','Off');
imshow(Fig3Img.cdata,'InitialMagnification',100);
Fig4ImgAxesHandle=gca;
axis off;
% set(Fig4ImgAxesHandle,'Position',[0 0 1 1],'xtick',[],'ytick',[]);

    
Fig4PolyROIHandle=impoly(Fig4ImgAxesHandle);
GRAutofluorPosition=wait(Fig4PolyROIHandle);
   
close(figHandle4);

GRAutofluorLogIntensityGate=[GRAutofluorPosition(:,1)*(log10(MaxG)-log10(MinG))/400,(400-GRAutofluorPosition(:,2))*(log10(MaxR)-log10(MinR))/400];


GRAutofluorGate_RemoveRowIndex=zeros(1,size(XFP_RGB,1));

for j=1:size(XFP_RGB,1)
            
            [InGRAutofluorGate OnGRAutofluorGate]=inpolygon(log10(XFP_RGB(j,2)),log10(XFP_RGB(j,1)),GRAutofluorLogIntensityGate(:,1),GRAutofluorLogIntensityGate(:,2));
            
            if (InGRAutofluorGate+OnGRAutofluorGate)<1
                GRAutofluorGate_RemoveRowIndex(1,j)=j;
            end;
            
end;


%% USER-SPECIFIED GATE 3: BLUE-RED AUTOFLUORESCENCE (POLYGON) GATE

fprintf('\nSpecify Blue-Red Autofluorescence (1xAF) polygon gate.\n');

fprintf('Specify Autofluorescence Region in Fig4: Blue vs Red Plot w/ polygon ROI. Double-click on polygon when done:\n\n')

Fig2Img=getframe(figHandle2,[50 50 400 400]);
figHandle4=figure('Position',[scrsz(3) 0 400 400],'OuterPosition',[scrsz(3) 0 400 400],'Resize','Off');
imshow(Fig2Img.cdata,'InitialMagnification',100);
Fig4ImgAxesHandle=gca;
axis off;
    
Fig4PolyROIHandle=impoly(Fig4ImgAxesHandle);
BRAutofluorPosition=wait(Fig4PolyROIHandle);
   
close(figHandle4);

BRAutofluorLogIntensityGate=[BRAutofluorPosition(:,1)*(log10(MaxB)-log10(MinB))/400,(400-BRAutofluorPosition(:,2))*(log10(MaxR)-log10(MinR))/400];

BRAutofluorGate_RemoveRowIndex=zeros(1,size(XFP_RGB,1));

for j=1:size(XFP_RGB,1)
            
            [InBRAutofluorGate OnBRAutofluorGate]=inpolygon(log10(XFP_RGB(j,3)),log10(XFP_RGB(j,1)),BRAutofluorLogIntensityGate(:,1),BRAutofluorLogIntensityGate(:,2));
            
            if (InBRAutofluorGate+OnBRAutofluorGate)<1
                BRAutofluorGate_RemoveRowIndex(1,j)=j;
            end;
            
end;


%% USER-SPECIFIED GATE 4: BLUE-GREEN AUTOFLUORESCENCE (POLYGON) GATE

fprintf('\nSpecify Blue-Green (1xAF) Autofluorescence polygon gate.\n');

fprintf('Specify Autofluorescence Region in Fig4: Blue vs Green Plot w/ polygon ROI. Double-click on polygon when done:\n\n')

Fig1Img=getframe(figHandle1,[50 50 400 400]);
figHandle4=figure('Position',[scrsz(3) 0 400 400],'OuterPosition',[scrsz(3) 0 400 400],'Resize','Off');
imshow(Fig1Img.cdata,'InitialMagnification',100);
Fig4ImgAxesHandle=gca;
axis off;
    
Fig4PolyROIHandle=impoly(Fig4ImgAxesHandle);
BGAutofluorPosition=wait(Fig4PolyROIHandle);
   
close(figHandle4);

BGAutofluorLogIntensityGate=[BGAutofluorPosition(:,1)*(log10(MaxB)-log10(MinB))/400,(400-BGAutofluorPosition(:,2))*(log10(MaxG)-log10(MinG))/400];

BGAutofluorGate_RemoveRowIndex=zeros(1,size(XFP_RGB,1));

for j=1:size(XFP_RGB,1)
            
            [InBGAutofluorGate OnBGAutofluorGate]=inpolygon(log10(XFP_RGB(j,3)),log10(XFP_RGB(j,2)),BGAutofluorLogIntensityGate(:,1),BGAutofluorLogIntensityGate(:,2));
            
            if (InBGAutofluorGate+OnBGAutofluorGate)<1
                BGAutofluorGate_RemoveRowIndex(1,j)=j;
            end;
            
end;

fprintf('\n\n');

close all;


%% 1xAF calculation 

AND_OR_AFQuery=1;
AutofluorGate_RemoveRowIndex=max([GRAutofluorGate_RemoveRowIndex;BRAutofluorGate_RemoveRowIndex;BGAutofluorGate_RemoveRowIndex]);

NoiseAFGate_RemoveRowIndex=max([NoiseGate_RemoveRowIndex;AutofluorGate_RemoveRowIndex]);
NoiseAFGate_RemoveRowIndex=nonzeros(NoiseAFGate_RemoveRowIndex');

XFP_RGB(NoiseAFGate_RemoveRowIndex,:)=[];
DataPtCt=size(XFP_RGB,1);

clearvars GRAutofluorGate_RemoveRowIndex BRAutofluorGate_RemoveRowIndex BGAutofluorGate_RemoveRowIndex NoiseGate_RemoveRowIndex;

RLogHistEdge=[0.00:0.1:6.00];
[RLogHistCount,RLogHistBin]=histc(log10(XFP_RGB(:,1)),RLogHistEdge);

BGLogHistEdge{1}=[0.00:0.1:6.00]; 
BGLogHistEdge{2}=[0.00:0.1:6.00]; 

BGHist3RLogStack=zeros([size(BGLogHistEdge{1},2),size(BGLogHistEdge{2},2),size(RLogHistEdge,2)]);
BGLogHist3Data_i=zeros([size(BGLogHistEdge{1},2),size(BGLogHistEdge{2},2)]);

for j=1:DataPtCt
    BGLogHist3Data_i=hist3([log10(XFP_RGB(j,3)),log10(XFP_RGB(j,2))],'Edges',BGLogHistEdge);  
    BGHist3RLogStack(:,:,RLogHistBin(j))=BGHist3RLogStack(:,:,RLogHistBin(j))+BGLogHist3Data_i;
end

BGHist3RLogStack_NaN=BGHist3RLogStack;
BGHist3RLogStack_NaN(BGHist3RLogStack_NaN==0)=NaN;

[MaxAFDensity MaxAFDensityIndexLoc]=max(BGHist3RLogStack(:));
[MaxAFDensity_Bsubscript,MaxAFDensity_Gsubscript,MaxAFDensity_Rsubscript]=ind2sub(size(BGHist3RLogStack),MaxAFDensityIndexLoc);

NoiseFree_AF_LogBmode=BGLogHistEdge{1}(MaxAFDensity_Bsubscript);
NoiseFree_AF_LogGmode=BGLogHistEdge{2}(MaxAFDensity_Gsubscript);
NoiseFree_AF_LogRmode=RLogHistEdge(MaxAFDensity_Rsubscript);

clearvars MaxAFDensity MaxAFDensityIndexLoc MaxAFDensity_Bsubscript MaxAFDensity_Gsubscript MaxAFDensity_Rsubscript;

fprintf(strcat('\nMaximum binned cell ct, autofluorescence:\t',num2str(max(BGHist3RLogStack(:))),'\n'));
AF_IsoSurf_Isovalue=input('Specify cell ct value to draw isosurface ("Enter" for default=0.02*max ct): ');

if isempty(AF_IsoSurf_Isovalue)
       AF_IsoSurf_Isovalue=0.02*max(BGHist3RLogStack(:));
end;

fprintf('\n\n\n');

BGHist3RLogStack_GBPermute_Smooth=smooth3(permute(BGHist3RLogStack,[2,1,3]),'gaussian',[3 3 3]);

AF_IsoSurf_VerticeFace=isosurface(BGLogHistEdge{1},BGLogHistEdge{2},RLogHistEdge,BGHist3RLogStack_GBPermute_Smooth,AF_IsoSurf_Isovalue);


%% FIGURE7 : PLOT 1xAF AUTOFLUORESCENCE 

AF_IsoSurf_ConvexHull_Index=convhull(AF_IsoSurf_VerticeFace.vertices(:,1),AF_IsoSurf_VerticeFace.vertices(:,2),AF_IsoSurf_VerticeFace.vertices(:,3));

% function inhull() from MATLAB file exchange. 
% Output: 1=inside convex hull. 2=outside.
XFP_RGB_IsoSurfConvHull_Query=inhull([log10(XFP_RGB(:,3)),log10(XFP_RGB(:,2)),log10(XFP_RGB(:,1))],[AF_IsoSurf_VerticeFace.vertices(:,1),AF_IsoSurf_VerticeFace.vertices(:,2),AF_IsoSurf_VerticeFace.vertices(:,3)],AF_IsoSurf_ConvexHull_Index);

scrszWscalefactor=1100/scrsz(3);
scrszHscalefactor=500/scrsz(4);

scrsz = get(0,'ScreenSize');
figHandle7=figure('Position',[150  scrsz(4)-(scrsz(4)*scrszHscalefactor)-100 scrsz(3)*scrszWscalefactor scrsz(4)*scrszHscalefactor],'Color','w','PaperPositionMode', 'auto');

set(gcf, 'Renderer', 'OpenGL');

subplot(1,2,1)
hold all;
axis on;
axis equal;
grid on;

axis([floor(min(min(log10(XFP_RGB(:,1:3))))) ceil(max(max(log10(XFP_RGB(:,1:3))))) floor(min(min(log10(XFP_RGB(:,1:3))))) ceil(max(max(log10(XFP_RGB(:,1:3))))) floor(min(min(log10(XFP_RGB(:,1:3))))) ceil(max(max(log10(XFP_RGB(:,1:3)))))]);
set(gca,'XColor',[0,0,0.5],'YColor',[0,0.5,0],'ZColor',[0.5,0,0]);

view([-30 30]); 

xlabel(strcat('B: log10 ( ',HeaderMtx{1}{ColorChoice(3)},' )'),'FontSize',10);
ylabel(strcat('G: log10 ( ',HeaderMtx{1}{ColorChoice(2)},' )'),'FontSize',10);
zlabel(strcat('R: log10 ( ',HeaderMtx{1}{ColorChoice(1)},' )'),'FontSize',10);
title({FileName;strcat('DENSITY ISOSURFACE (Isovalue=',num2str(AF_IsoSurf_Isovalue),') CONVEX HULL')});
titlehandle=get(gca,'title');
set(titlehandle, 'FontSize', 8);

Fig7A_ConvexHullHandle=trisurf(AF_IsoSurf_ConvexHull_Index,AF_IsoSurf_VerticeFace.vertices(:,1),AF_IsoSurf_VerticeFace.vertices(:,2),AF_IsoSurf_VerticeFace.vertices(:,3));

set(Fig7A_ConvexHullHandle,'FaceColor','red','EdgeColor','none');
camlight;
lighting phong;

hold off;


subplot(1,2,2)
hold all;
axis on;
axis equal;
grid on;

axis([floor(min(min(log10(XFP_RGB(:,1:3))))) ceil(max(max(log10(XFP_RGB(:,1:3))))) floor(min(min(log10(XFP_RGB(:,1:3))))) ceil(max(max(log10(XFP_RGB(:,1:3))))) floor(min(min(log10(XFP_RGB(:,1:3))))) ceil(max(max(log10(XFP_RGB(:,1:3)))))]);
set(gca,'XColor',[0,0,0.5],'YColor',[0,0.5,0],'ZColor',[0.5,0,0]);

view([-30 30]); 

xlabel(strcat('B: log10 ( ',HeaderMtx{1}{ColorChoice(3)},' )'),'FontSize',10);
ylabel(strcat('G: log10 ( ',HeaderMtx{1}{ColorChoice(2)},' )'),'FontSize',10);
zlabel(strcat('R: log10 ( ',HeaderMtx{1}{ColorChoice(1)},' )'),'FontSize',10);
title({FileName;strcat('DENSITY ISOSURFACE (Isovalue=',num2str(AF_IsoSurf_Isovalue),') CONVEX HULL')});
titlehandle=get(gca,'title');
set(titlehandle, 'FontSize', 8);

Fig7B_ConvexHullHandle=trisurf(AF_IsoSurf_ConvexHull_Index,AF_IsoSurf_VerticeFace.vertices(:,1),AF_IsoSurf_VerticeFace.vertices(:,2),AF_IsoSurf_VerticeFace.vertices(:,3));

set(Fig7B_ConvexHullHandle,'FaceColor','none','EdgeColor','red');
camlight;
lighting phong;

XFP_NoiseFreeAF_DataPtSize=3;
XFP_NoiseFreeAF_DataPtColor=[zeros(DataPtCt,2),abs(XFP_RGB_IsoSurfConvHull_Query-1)];

Fig7B_ScatterPlotHandle=scatter3(log10(XFP_RGB(:,3)),log10(XFP_RGB(:,2)),log10(XFP_RGB(:,1)),XFP_NoiseFreeAF_DataPtSize,XFP_NoiseFreeAF_DataPtColor,'filled');

hold off;

print(figHandle7,'-dtiffn','-r72',strcat(OutputFolderString,'/',FileName,' AFIsosurf',num2str(AF_IsoSurf_Isovalue),' 1xAFplot'));


%% Save 1xAF properties

AF_IsoSurfBGRLogVerticesFaces_ConvexHull=NaN(max(size(AF_IsoSurf_VerticeFace.vertices,1),size(AF_IsoSurf_VerticeFace.faces,1)),9);

AF_IsoSurfBGRLogVerticesFaces_ConvexHull(1:size(AF_IsoSurf_VerticeFace.vertices,1),1:3)=AF_IsoSurf_VerticeFace.vertices;
AF_IsoSurfBGRLogVerticesFaces_ConvexHull(1:size(AF_IsoSurf_VerticeFace.faces,1),4:6)=AF_IsoSurf_VerticeFace.faces;
AF_IsoSurfBGRLogVerticesFaces_ConvexHull(1:size(AF_IsoSurf_ConvexHull_Index,1),7:9)=AF_IsoSurf_ConvexHull_Index;

AF_IsoSurfBGRLogVerticesFaces_ConvexHull_Row1Centroid=horzcat(AF_IsoSurfBGRLogVerticesFaces_ConvexHull(1,:),NoiseFree_AF_LogBmode,NoiseFree_AF_LogGmode,NoiseFree_AF_LogRmode);

AF_IsoSurfConvexHull_HeaderRow={strcat('AFIsosurfVertices Log10(B): ',HeaderMtx{1}{ColorChoice(3)});strcat('AFIsosurfVertices Log10(G): ',HeaderMtx{1}{ColorChoice(2)});strcat('AFIsosurfVertices Log10(R): ',HeaderMtx{1}{ColorChoice(1)});
                                                                    strcat('AFIsosurfFaces: IsosurfVertexBGR1');strcat('AFIsosurfFaces: IsosurfVertexBGR2');strcat('AFIsoSurfFaces: IsosurfVertexBGR3');
                                                                    strcat('AFIsosurf ConvexHull: IsosurfVertexBGR1');strcat('AFIsosurf ConvexHull: IsosurfVertexBGR2');strcat('AFIsosurf ConvexHull: IsosurfVertexBGR3');
                                                                    strcat('AFMode Log10(B): ',HeaderMtx{1}{ColorChoice(3)});strcat('AFMode Log10(G): ',HeaderMtx{1}{ColorChoice(2)});strcat('AFMode Log10(R): ',HeaderMtx{1}{ColorChoice(1)})};
AF_IsoSurfConvexHull_FileNameString=strcat(OutputFolderString,'/',FileName,' 1xAF.txt');


fid=fopen(AF_IsoSurfConvexHull_FileNameString,'w');

for i=1:numel(AF_IsoSurfConvexHull_HeaderRow);
    fprintf(fid,'%s\t',AF_IsoSurfConvexHull_HeaderRow{i});
end

fprintf(fid,'\n');

dlmwrite(AF_IsoSurfConvexHull_FileNameString,AF_IsoSurfBGRLogVerticesFaces_ConvexHull_Row1Centroid,'delimiter','\t','-append');
dlmwrite(AF_IsoSurfConvexHull_FileNameString,AF_IsoSurfBGRLogVerticesFaces_ConvexHull(2:end,:),'delimiter','\t','-append');

fclose(fid);

clearvars AF_IsoSurfBGRLogVerticesFaces_ConvexHull AF_IsoSurfConvexHull_FileNameString AF_IsoSurfBGRLogVerticesFaces_ConvexHull_Row1Centroid  AF_IsoSurfConvexHull_FileNameString AF_IsoSurfConvexHull_HeaderRow;


 clc; close all hidden; clear all; clear java; fclose('all');

