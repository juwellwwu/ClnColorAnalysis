function [Rmode,Gmode,Bmode,Mode_THETA,Mode_PHI,Mode_RADIUS]=cmodspd(FileNameString,RchannelChoice,GchannelChoice,BchannelChoice,IsovalFrac,AFTieringQuery,AFMultiplierLoQuery,AFMultiplierHiQuery)
 
%% =====DESCRIPTION=====

% Calculate, save clonal chromatic mode and spread of one clone

% ==Usage: 
% User specifies variables in "USER INPUT" section.
% User specifies output file location for "*3DHistgm.tif",
% "*In2DChromSpread.txt", "*ChromModeSpdProperty.txt"

% ==Output file: "*3DHistgm.tif'"
% Plot 3D histogram of color data in RGB soace

% ==Output file: "*In2DChromSpread.txt"
% Save list of color data pts (in THETA-PHI deg) inside 2D chromatic spread 

% ==Output file: "*ChromModeSpdProperty.txt"
% Save 3D isosurface and 2D chromatic spread properties and chromatic mode values 
% Input for script "chromspdanalysis_batch.m" (prepare clone masks for clonal assignment, analyze single clone chromatic stability)


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUT

% 'y' to plot color data 3D histogram 
Fig4HistBGRPlotQuery='n';

% 'y' to plot chromatic spread on THETA-PHI grid and color data
Fig9CvxHull2DPlotQuery='y';
DataPtPlotFig=2.5E4; % Approx # data pts to plot, -1 to plot all

% 'y' to save list of color data pts inside 2D chromatic spread
SaveInHullRGB='n';

% 'y' to save 3D isosurface + 2D chromatic spread properties
SaveCvxHullProperty='y';
   

%% Load Color Data

[~,FileName,FileExt]=fileparts(FileNameString);

FileID=fopen(FileNameString);

HeaderRow=fgets(FileID);
HeaderMtx=textscan(HeaderRow,'%s','delimiter','\t');

DataPtCt=0;
while (fgets(FileID) ~= -1),
  DataPtCt=DataPtCt+1;
end

for i=1:numel(HeaderMtx{1})
    i_str=num2str(i);
    dispstring=[i_str,': ',HeaderMtx{1}{i}];
    disp(dispstring)
end

ColorChoice(1)=RchannelChoice;
ColorChoice(2)=GchannelChoice;
ColorChoice(3)=BchannelChoice;

fprintf('\n');

fprintf(strcat('R set to:\t',num2str(ColorChoice(1)),'\n'));
fprintf(strcat('G set to:\t',num2str(ColorChoice(2)),'\n'));
fprintf(strcat('B set to:\t',num2str(ColorChoice(3)),'\n'));

XFP_RGB=zeros(DataPtCt,4,'double');

for i=1:3
    if ColorChoice(i)<9999 
        XFP_RGB(:,i)=dlmread(FileNameString,'\t',[1,ColorChoice(i)-1,DataPtCt,ColorChoice(i)-1]);
    else
        XFP_RGB(:,i)=zeros(DataPtCt,1,'double');
        HeaderMtx{1}{ColorChoice(i)}=' ';
    end
end


%% Relative cell brightness, <1xAF exclusion

FlagLessThanHullQuery='';

if isempty(AFTieringQuery)
        AFTieringQuery='n';
end;

if AFTieringQuery=='y'

    Mode0_XFP_RGB_AFMultiplier=dlmread(FileNameString,'\t',[1,(numel(HeaderMtx{1})-3)-1,size(XFP_RGB,1),(numel(HeaderMtx{1})-3)-1]);

    FlagLessThanHullQuery='y';
    Flag_LessThanConvexHull=dlmread(FileNameString,'\t',[1,(numel(HeaderMtx{1})-2)-1,size(XFP_RGB,1),(numel(HeaderMtx{1})-2)-1]);
    XFP_RGB(find(Flag_LessThanConvexHull),:)=[];
    Mode0_XFP_RGB_AFMultiplier(find(Flag_LessThanConvexHull),:)=[];

    XFP_RGB(:,4)=Mode0_XFP_RGB_AFMultiplier;
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

clearvars Mode0_XFP_RGB_AFMultiplier AFMultiplier_Hi_GreaterThanIndex AFMultiplier_Lo_LessThanIndex  Flag_LessThanConvexHull AFMultiplier_RemoveIndex;


%% Remove data points near 0 or channel saturation (262144)

RGB_SatRemoveIndex=zeros(size(XFP_RGB,1),1);

RGB_SatRemoveIndex=bsxfun(@gt,XFP_RGB(:,1),2E5)+bsxfun(@lt,XFP_RGB(:,1),10)+bsxfun(@gt,XFP_RGB(:,2),2E5)+bsxfun(@lt,XFP_RGB(:,2),10)+bsxfun(@gt,XFP_RGB(:,3),2E5)+bsxfun(@lt,XFP_RGB(:,3),10);

XFP_RGB(find(RGB_SatRemoveIndex),:)=[];

fprintf(strcat('\n\n# Data points after near-0 or saturated R/G/B intensities removed: \t',num2str(size(XFP_RGB,1)),'\n\n\n')); 

DataPtCt=size(XFP_RGB,1);

clearvars RGB_SatRemoveIndex;


%% Prepare Color Data for RGB Binning 

RbinCtChoice=50;
GbinCtChoice=50;
BbinCtChoice=50;

GaussXFP_RGB=XFP_RGB.^(1/4)*100;
GaussRmin=min(GaussXFP_RGB(:,1));
GaussRmax=max(GaussXFP_RGB(:,1));
if round((GaussRmax-GaussRmin)/RbinCtChoice,3)<0.001
    GaussRbinEdge=GaussRmin:0.001:GaussRmax;
else
    GaussRbinEdge=GaussRmin:round((GaussRmax-GaussRmin)/RbinCtChoice,3):GaussRmax;
end;
[GaussRhistct,GaussRbinPos]=hist(GaussXFP_RGB(:,1),GaussRbinEdge);

GaussGmin=min(GaussXFP_RGB(:,2));
GaussGmax=max(GaussXFP_RGB(:,2));
if round((GaussGmax-GaussGmin)/GbinCtChoice,3)<0.001
    GuassGbinEdge=GaussGmin:0.001:GaussGmax;
else
    GaussGbinEdge=GaussGmin:floor((GaussGmax-GaussGmin)/GbinCtChoice):GaussGmax;
end;
[GaussGhistct,GaussGbinPos]=hist(GaussXFP_RGB(:,2),GaussGbinEdge);

GaussBmin=min(GaussXFP_RGB(:,3));
GaussBmax=max(GaussXFP_RGB(:,3));
if round((GaussBmax-GaussBmin)/BbinCtChoice,3)<0.001
    GaussBbinEdge=GaussBmin:0.001:GaussBmax;
else
    GaussBbinEdge=GaussBmin:floor((GaussBmax-GaussBmin)/BbinCtChoice):GaussBmax;
end;
[GaussBhistct,GaussBbinPos]=hist(GaussXFP_RGB(:,3),GaussBbinEdge);


GaussRhistGs1FitObj=fit(GaussRbinPos',GaussRhistct','gauss1');
GaussRhistGs1Fitcoeff=coeffvalues(GaussRhistGs1FitObj);

GaussGhistGs1FitObj=fit(GaussGbinPos',GaussGhistct','gauss1');
GaussGhistGs1Fitcoeff=coeffvalues(GaussGhistGs1FitObj);

GaussBhistGs1FitObj=fit(GaussBbinPos',GaussBhistct','gauss1');
GaussBhistGs1Fitcoeff=coeffvalues(GaussBhistGs1FitObj);

 MaxHistctGauss1D=max([max(GaussRhistct),max(GaussGhistct),max(GaussBhistct)]);

RGB_OutlierRemoveIndex=zeros(size(XFP_RGB,1),1);

RGB_OutlierRemoveIndex=bsxfun(@gt,GaussXFP_RGB(:,1),GaussRhistGs1Fitcoeff(1,2)+GaussRhistGs1Fitcoeff(1,3)*4/sqrt(2))+bsxfun(@lt,GaussXFP_RGB(:,1),GaussRhistGs1Fitcoeff(1,2)-GaussRhistGs1Fitcoeff(1,3)*4/sqrt(2))+bsxfun(@gt,GaussXFP_RGB(:,2),GaussGhistGs1Fitcoeff(1,2)+GaussGhistGs1Fitcoeff(1,3)*4/sqrt(2))+bsxfun(@lt,GaussXFP_RGB(:,2),GaussGhistGs1Fitcoeff(1,2)-GaussGhistGs1Fitcoeff(1,3)*4/sqrt(2))+bsxfun(@gt,GaussXFP_RGB(:,3),GaussBhistGs1Fitcoeff(1,2)+GaussBhistGs1Fitcoeff(1,3)*4/sqrt(2))+bsxfun(@lt,GaussXFP_RGB(:,3),GaussBhistGs1Fitcoeff(1,2)-GaussBhistGs1Fitcoeff(1,3)*4/sqrt(2));

XFP_RGB(find(RGB_OutlierRemoveIndex),:)=[];

DataPtCt=size(XFP_RGB,1);

clearvars RGB_OutlierRemoveIndex;
clearvars GaussXFP_RGB;


Rmin=min(XFP_RGB(:,1));
Rmax=max(XFP_RGB(:,1));
if round((Rmax-Rmin)/RbinCtChoice,1)<0.1
    RbinEdge=Rmin:0.1:Rmax;
else
    RbinEdge=Rmin:round((Rmax-Rmin)/RbinCtChoice,1):Rmax;
end;
[Rhistct,RbinPos]=hist(XFP_RGB(:,1),RbinEdge);

Gmin=min(XFP_RGB(:,2));
Gmax=max(XFP_RGB(:,2));
if round((Gmax-Gmin)/GbinCtChoice,1)<0.1
    GbinEdge=Gmin:0.1:Gmax;
else
    GbinEdge=Gmin:floor((Gmax-Gmin)/GbinCtChoice):Gmax;
end;
[Ghistct,GbinPos]=hist(XFP_RGB(:,2),GbinEdge);


Bmin=min(XFP_RGB(:,3));
Bmax=max(XFP_RGB(:,3));
if round((Bmax-Bmin)/BbinCtChoice,1)<0.1
    BbinEdge=Bmin:0.1:Bbmax;
else
    BbinEdge=Bmin:floor((Bmax-Bmin)/BbinCtChoice):Bmax;
end;
[Bhistct,BbinPos]=hist(XFP_RGB(:,3),BbinEdge);

 MaxHistct1D=max([max(Rhistct),max(Ghistct),max(Bhistct)]);

 RGBmin=min([Rmin,Gmin,Bmin]);
 RGBmax=max([Rmax,Gmax,Bmax]);


%% FILENAME Designate 

if AFTieringQuery=='y'
        FileNameSpec=strcat(FileName,' RCBLo',num2str(AFMultiplierLoQuery),'Hi',num2str(AFMultiplierHiQuery));
else
    FileNameSpec=strcat(FileName,' R',num2str(ColorChoice(1)),'G',num2str(ColorChoice(2)),'B',num2str(ColorChoice(3)));
end;


%% Generate 3D RGB Histogram for Chromatic mode, spread calculation

 Rstep=floor((Rmax-Rmin)/50);
 Gstep=floor((Gmax-Gmin)/50);
 Bstep=floor((Bmax-Bmin)/50);

Rmin=floor(Rmin/Rstep)*Rstep;
Rmax=ceil(Rmax/Rstep)*Rstep;

Gmin=floor(Gmin/Gstep)*Gstep;
Gmax=ceil(Gmax/Gstep)*Gstep;

Bmin=floor(Bmin/Bstep)*Bstep;
Bmax=ceil(Bmax/Bstep)*Bstep;

RHistEdge=[Rmin:Rstep:Rmax];
[RHistCount,RHistBin]=histc(XFP_RGB(:,1),RHistEdge);

BGHistEdge{1}=[Bmin:Bstep:Bmax]; 
BGHistEdge{2}=[Gmin:Gstep:Gmax];

if matlabpool('size')>0
    XFP_RGB_DataPtBlockNum=matlabpool('size');
else
    XFP_RGB_DataPtBlockNum=1;
end;


if XFP_RGB_DataPtBlockNum>1 && DataPtCt>5000

        XFP_RGB_DataPtBlockRowCt=floor(DataPtCt/XFP_RGB_DataPtBlockNum);
        for j=1:XFP_RGB_DataPtBlockNum
            if j<XFP_RGB_DataPtBlockNum
                XFP_RGB_DataPtBlockRowStartEnd(j,:)=[XFP_RGB_DataPtBlockRowCt*(j-1)+1,XFP_RGB_DataPtBlockRowCt*j];
            elseif j==XFP_RGB_DataPtBlockNum
                XFP_RGB_DataPtBlockRowStartEnd(j,:)=[XFP_RGB_DataPtBlockRowCt*(j-1)+1,DataPtCt];
            end;
        end;

         BGHist3RStack_Array=cell(XFP_RGB_DataPtBlockNum,1);
         BGHist3Data_i_Array=cell(XFP_RGB_DataPtBlockNum,1);

         for  j=1:XFP_RGB_DataPtBlockNum
             BGHist3RStack_Array{j}=zeros([size(BGHistEdge{1},2),size(BGHistEdge{2},2),size(RHistEdge,2)]);
             BGHist3Data_i_Array{j}=zeros([size(BGHistEdge{1},2),size(BGHistEdge{2},2)]);
         end;

        parfor j=1:XFP_RGB_DataPtBlockNum
            for jj=XFP_RGB_DataPtBlockRowStartEnd(j,1):XFP_RGB_DataPtBlockRowStartEnd(j,2)
                BGHist3Data_i_Array{j}=hist3([XFP_RGB(jj,3),XFP_RGB(jj,2)],'Edges',BGHistEdge);  
                BGHist3RStack_Array{j}(:,:,RHistBin(jj))=BGHist3RStack_Array{j}(:,:,RHistBin(jj))+BGHist3Data_i_Array{j};
            end;
        end;

        BGHist3RStack=zeros([size(BGHistEdge{1},2),size(BGHistEdge{2},2),size(RHistEdge,2)]);
        for j=1:XFP_RGB_DataPtBlockNum
            BGHist3RStack= BGHist3RStack+BGHist3RStack_Array{j};
        end;

        clearvars BGHist3Data_i_Array BGHist3RStack_Array RHistBin;

else

        BGHist3RStack=zeros([size(BGHistEdge{1},2),size(BGHistEdge{2},2),size(RHistEdge,2)]);
        BGHist3Data_i=zeros([size(BGHistEdge{1},2),size(BGHistEdge{2},2)]);

        for j=1:DataPtCt
            BGHist3Data_i=hist3([XFP_RGB(j,3),XFP_RGB(j,2)],'Edges',BGHistEdge);  
            BGHist3RStack(:,:,RHistBin(j))=BGHist3RStack(:,:,RHistBin(j))+BGHist3Data_i;
        end;

        clearvars BGHist3Data_i RHistBin;

end;


%% FIGURE4 : 3D RGB HISTOGRAM 

if isempty(Fig4HistBGRPlotQuery)
       Fig4HistBGRPlotQuery='n';
end;

if Fig4HistBGRPlotQuery=='y';

    BGHist3RStack_NaN=BGHist3RStack;
    BGHist3RStack_NaN(BGHist3RStack_NaN==0)=NaN;

    % Uncomment if want to plot bins > certain cell ct.
    % BGHist3RStack_NaN(BGHist3RStack_NaN<326)=NaN;

    scrsz = get(0,'ScreenSize');
    scrszWscalefactor=1200/scrsz(3);
    scrszHscalefactor=900/scrsz(4);

    figHandle4=figure('Position',[scrsz(3)-550  scrsz(4)-(scrsz(4)*scrszHscalefactor)-100 scrsz(3)*scrszWscalefactor scrsz(4)*scrszHscalefactor],'Color','w','PaperPositionMode', 'auto','Visible','off');

    hold on;
    axis on
    axis equal;

    grid off

    TargetColor=[1 0 0];                 
    TargetColorMapRes=100;                      
    TargetColorMap=zeros(3,TargetColorMapRes+1);
    TargetColorMap=horzcat(([TargetColorMapRes:-(TargetColorMapRes-(1-TargetColor(1,1))*TargetColorMapRes)/TargetColorMapRes:(1-TargetColor(1,1))*TargetColorMapRes]/TargetColorMapRes)',([TargetColorMapRes:-(TargetColorMapRes-(1-TargetColor(1,2))*TargetColorMapRes)/TargetColorMapRes:(1-TargetColor(1,2))*TargetColorMapRes]/TargetColorMapRes)',([TargetColorMapRes:-(TargetColorMapRes-(1-TargetColor(1,3))*TargetColorMapRes)/TargetColorMapRes:(1-TargetColor(1,3))*TargetColorMapRes]/TargetColorMapRes)');                       
    TargetColorMap=bsxfun(@minus,1,TargetColorMap);                    
    colormap(TargetColorMap);

    %  PATCH_3Darray from MATLAB File Exchange.
    [Fig4PATCHhandle,Fig4Colorbarhandle]=PATCH_3Darray_13Nov04Juwellmod(BGHist3RStack_NaN,BGHistEdge{1},BGHistEdge{2},RHistEdge,'col','cmap',TargetColorMap,'clim',[0,max(BGHist3RStack(:))]);  

    colorbar('EastOutside');

    axis([0 2E5 0 2E5 0 2E5]);

    set(gca,'XColor',[0,0,0.5],'YColor',[0,0.5,0],'ZColor',[0.5,0,0]);
   
    view([75 -15]); % 3D Cartesian view
   
    xlabel(strcat('B: ',HeaderMtx{1}{ColorChoice(3)}),'FontSize',10);
    ylabel(strcat('G: ',HeaderMtx{1}{ColorChoice(2)}),'FontSize',10);
    zlabel(strcat('R: ',HeaderMtx{1}{ColorChoice(1)}),'FontSize',10);
    title({strcat(FileName,'; Cell Ct: ',num2str(DataPtCt),' 3D RGB HISTOGRAM')});
    titlehandle=get(gca,'title');
    set(titlehandle, 'FontSize', 6);

    colorbar('EastOutside');

    print(figHandle4,'-dtiffn','-r150',strcat('ClnColorDEMO/cmodspd output/Isofrac',sprintf('%.2f',IsovalFrac),'/',FileNameSpec,' 3DHistgm.tif'));

    hold off

    close(figHandle4);

    hold off

   clearvars BGHist3RStack_NaN

end;


%% Locate Chromatic Mode 

[MaxDensity MaxDensityIndexLoc]=max(BGHist3RStack(:));
[MaxDensity_Bsubscript,MaxDensity_Gsubscript,MaxDensity_Rsubscript]=ind2sub(size(BGHist3RStack),MaxDensityIndexLoc);

Bmode=BGHistEdge{1}(MaxDensity_Bsubscript);
Gmode=BGHistEdge{2}(MaxDensity_Gsubscript);
Rmode=RHistEdge(MaxDensity_Rsubscript);

[Mode_THETA,Mode_PHI,Mode_RADIUS] = cart2sph(Bmode,Gmode,Rmode);

fprintf(strcat('\nMaximum binned cell ct:\t',num2str(max(BGHist3RStack(:))),'\n\n\n\n'));

clearvars MaxDensity MaxDensityIndexLoc MaxDensity_Bsubscript MaxDensity_Gsubscript MaxDensity_Rsubscript 



%% 2D Chromatic Spread: Calculate, plot 

IsoSurf_Isovalue=max(BGHist3RStack(:))*IsovalFrac;

BGHist3RStack_GBPermute_Smooth=smooth3(permute(BGHist3RStack,[2,1,3]),'gaussian',[3 3 3]);
IsoSurf_VerticeFace=isosurface(BGHistEdge{1},BGHistEdge{2},RHistEdge,BGHist3RStack_GBPermute_Smooth,IsoSurf_Isovalue);

IsoSurf_ConvexHull_Index=convhull(IsoSurf_VerticeFace.vertices(:,1),IsoSurf_VerticeFace.vertices(:,2),IsoSurf_VerticeFace.vertices(:,3));

IsoSurf_VerticeFace.verticesTHETArad=zeros(size(IsoSurf_VerticeFace.vertices(:,1)));
IsoSurf_VerticeFace.verticesTHETAdeg=zeros(size(IsoSurf_VerticeFace.vertices(:,1)));
IsoSurf_VerticeFace.verticesPHIrad=zeros(size(IsoSurf_VerticeFace.vertices(:,1)));
IsoSurf_VerticeFace.verticesPHIdeg=zeros(size(IsoSurf_VerticeFace.vertices(:,1)));

[IsoSurf_VerticeFace.verticesTHETArad,IsoSurf_VerticeFace.verticesPHIrad,IsoSurf_VerticeFace.verticesRADIUS]=cart2sph(IsoSurf_VerticeFace.vertices(:,1),IsoSurf_VerticeFace.vertices(:,2),IsoSurf_VerticeFace.vertices(:,3));

IsoSurf_VerticeFace.verticesTHETAdeg=180/pi().*(IsoSurf_VerticeFace.verticesTHETArad);
IsoSurf_VerticeFace.verticesPHIdeg=180/pi().*(IsoSurf_VerticeFace.verticesPHIrad);

IsoSurf_THETAPHIConvexHull_Index=convhull(IsoSurf_VerticeFace.verticesTHETAdeg,IsoSurf_VerticeFace.verticesPHIdeg);

XFP_RGB_IsoSurfTHETAPHIConvHull_Query=zeros(size(XFP_RGB,1),1);
XFP_THETArad=zeros(size(XFP_RGB,1),1);
XFP_PHIrad=zeros(size(XFP_RGB,1),1);
XFP_RADIUS=zeros(size(XFP_RGB,1),1);

[XFP_THETArad,XFP_PHIrad,XFP_RADIUS]=cart2sph(XFP_RGB(:,3),XFP_RGB(:,2),XFP_RGB(:,1));

XFP_RGB_IsoSurfTHETAPHIConvHull_Query=inpolygon(180/pi().*(XFP_THETArad),180/pi().*(XFP_PHIrad),IsoSurf_VerticeFace.verticesTHETAdeg(IsoSurf_THETAPHIConvexHull_Index),IsoSurf_VerticeFace.verticesPHIdeg(IsoSurf_THETAPHIConvexHull_Index));


% PLOT
if Fig9CvxHull2DPlotQuery=='y';

            scrsz = get(0,'ScreenSize');
            scrszWscalefactor=800/scrsz(3);
            scrszHscalefactor=600/scrsz(4);

            figHandle9=figure('Position',[25  scrsz(4)-(scrsz(4)*scrszHscalefactor)-100 scrsz(3)*scrszWscalefactor scrsz(4)*scrszHscalefactor],'Color','w','PaperPositionMode', 'auto');

            set(gcf, 'Renderer', 'OpenGL');


            % Figure 9
            hold all;
            axis on;
            axis square;
            grid on;

            axis([0 90 0 90]);
            set(gca,'XColor',[0,0,0],'YColor',[0,0,0]);

            xlabel(strcat('THETA'),'FontSize',8);
            ylabel(strcat('PHI'),'FontSize',8);

            Fig9_THETAPHIConvexHullHandle=plot(IsoSurf_VerticeFace.verticesTHETAdeg(IsoSurf_THETAPHIConvexHull_Index),IsoSurf_VerticeFace.verticesPHIdeg(IsoSurf_THETAPHIConvexHull_Index));
            set(Fig9_THETAPHIConvexHullHandle,'Color','red','LineWidth',2)

            % Color data within chromatic spread: black; outside: blue.
            XFP_DataPtSize=1;
            XFP_DataPtTPColor=[zeros(DataPtCt,2),abs(XFP_RGB_IsoSurfTHETAPHIConvHull_Query-1)];

            PrintDataPtSkip=max(1,round(DataPtCt/DataPtPlotFig));

            Fig9_ScatterPlotHandle=scatter(180/pi().*(XFP_THETArad(1:PrintDataPtSkip:DataPtCt)),180/pi().*(XFP_PHIrad(1:PrintDataPtSkip:DataPtCt)),XFP_DataPtSize,XFP_DataPtTPColor(1:PrintDataPtSkip:DataPtCt,:),'filled');

            title({FileNameSpec;strcat('#DataPt in 2D Chromatic Spread (Isovalue ',num2str(IsoSurf_Isovalue),'): ',num2str(sum(XFP_RGB_IsoSurfTHETAPHIConvHull_Query)));strcat('#DataPt printed/total: ',num2str(size(XFP_RGB(1:PrintDataPtSkip:DataPtCt,3),1)),'/',num2str(DataPtCt))});
            titlehandle=get(gca,'title');
            set(titlehandle, 'FontSize',7); 

            print(figHandle9,'-dtiffn',strcat('ClnColorDEMO/cmodspd output/Isofrac',sprintf('%.2f',IsovalFrac),'/',FileNameSpec,' Isoval',num2str(IsoSurf_Isovalue),' ChromSpread.tif'));

            hold off;

            close(figHandle9);
end;

clearvars XFP_DataPtTPColor;


close('all');
close all hidden;



%% EXPORT COLOR DATA INSIDE 2D CHROMATIC SPREAD

if SaveInHullRGB=='y'

             if AFTieringQuery=='y'
                XFP_RGB_InTPHull_HeaderRow={strcat('In2DChromSpread THETAdeg (B-G:',HeaderMtx{1}{ColorChoice(3)},'-',HeaderMtx{1}{ColorChoice(2)},')');strcat('In2DChromSpread PHIdeg (R:',HeaderMtx{1}{ColorChoice(1)},')');strcat('RelCellBright (xAF)')};
            else
                XFP_RGB_InTPHull_HeaderRow={strcat('In2DChromSpread THETAdeg (B-G:',HeaderMtx{1}{ColorChoice(3)},'-',HeaderMtx{1}{ColorChoice(2)},')');strcat('In2DChromSpread PHIdeg (R:',HeaderMtx{1}{ColorChoice(1)},')')};
            end;

            XFP_RGB_InTPHull_FileNameString=strcat('ClnColorDEMO/cmodspd output/Isofrac',sprintf('%.2f',IsovalFrac),'/',FileName,' Isoval',num2str(IsoSurf_Isovalue),'In2DChromSpread.txt');

            fid=fopen(XFP_RGB_InTPHull_FileNameString,'w');

            for i=1:numel(XFP_RGB_InTPHull_HeaderRow);
                fprintf(fid,'%s\t',XFP_RGB_InTPHull_HeaderRow{i});
            end

            fprintf(fid,'\n');

            if AFTieringQuery=='y'
                XFP_THETAPHIdeg=horzcat(180/pi().*(XFP_THETArad),180/pi().*(XFP_PHIrad),XFP_RGB(:,4));
            else
                 XFP_THETAPHIdeg=horzcat(180/pi().*(XFP_THETArad),180/pi().*(XFP_PHIrad));
            end;

            dlmwrite(XFP_RGB_InTPHull_FileNameString,XFP_THETAPHIdeg(find(XFP_RGB_IsoSurfTHETAPHIConvHull_Query),:),'delimiter','\t','-append','precision',16);

            fclose(fid);

            clearvars XFP_RGB XFP_THETArad XFP_PHIrad XFP_THETAPHIdeg 

end;


%% EXPORT 3D ISOSURFACE & 2D CHROMATIC SPREAD PROPERTIES 

if SaveCvxHullProperty=='y'

            IsoSurfBGRVerticesFaces_ConvexHull=NaN(max(size(IsoSurf_VerticeFace.vertices,1),size(IsoSurf_VerticeFace.faces,1)),11);

            IsoSurfBGRVerticesFaces_ConvexHull(1:size(IsoSurf_VerticeFace.vertices,1),1:3)=IsoSurf_VerticeFace.vertices;
            IsoSurfBGRVerticesFaces_ConvexHull(1:size(IsoSurf_VerticeFace.faces,1),4:6)=IsoSurf_VerticeFace.faces;
            IsoSurfBGRVerticesFaces_ConvexHull(1:size(IsoSurf_ConvexHull_Index,1),7:9)=IsoSurf_ConvexHull_Index;
            IsoSurfBGRVerticesFaces_ConvexHull(1:size(IsoSurf_THETAPHIConvexHull_Index,1),10)=IsoSurf_VerticeFace.verticesTHETAdeg(IsoSurf_THETAPHIConvexHull_Index);
            IsoSurfBGRVerticesFaces_ConvexHull(1:size(IsoSurf_THETAPHIConvexHull_Index,1),11)=IsoSurf_VerticeFace.verticesPHIdeg(IsoSurf_THETAPHIConvexHull_Index);

            IsoSurfBGRVerticesFaces_ConvexHull_Row1Centroid=horzcat(IsoSurfBGRVerticesFaces_ConvexHull(1,:),Bmode,Gmode,Rmode);


            IsoSurfConvexHull_HeaderRow={strcat('IsosurfVertices B: ',HeaderMtx{1}{ColorChoice(3)});strcat('IsosurfVertices G: ',HeaderMtx{1}{ColorChoice(2)});strcat('IsosurfVertices R: ',HeaderMtx{1}{ColorChoice(1)});
                                                                                strcat('IsosurfFaces: IsosurfVertexBGR1');strcat('IsosurfFaces: IsosurfVertexBGR2');strcat('IsoSurfFaces: IsosurfVertexBGR3');
                                                                                strcat('Isosurf ConvexHull: IsosurfVertexBGR1');strcat('Isosurf ConvexHull: IsosurfVertexBGR2');strcat('Isosurf ConvexHull: IsosurfVertexBGR3');
                                                                                strcat('2DChromSpread: THETAdeg');strcat('2DChromSpread: PHIdeg');
                                                                                strcat('Mode B: ',HeaderMtx{1}{ColorChoice(3)});strcat('Mode G: ',HeaderMtx{1}{ColorChoice(2)});strcat('Mode R: ',HeaderMtx{1}{ColorChoice(1)})};
            IsoSurfConvexHull_FileNameString=strcat('ClnColorDEMO/cmodspd output/Isofrac',sprintf('%.2f',IsovalFrac),'/',FileName,' Isoval',num2str(IsoSurf_Isovalue),' ChromModeSpdProperty.txt');


            fid=fopen(IsoSurfConvexHull_FileNameString,'w');

            for i=1:numel(IsoSurfConvexHull_HeaderRow);
                fprintf(fid,'%s\t',IsoSurfConvexHull_HeaderRow{i});
            end

            fprintf(fid,'\n');

            dlmwrite(IsoSurfConvexHull_FileNameString,IsoSurfBGRVerticesFaces_ConvexHull_Row1Centroid,'delimiter','\t','-append','precision',16);
            dlmwrite(IsoSurfConvexHull_FileNameString,IsoSurfBGRVerticesFaces_ConvexHull(2:end,:),'delimiter','\t','-append','precision',16);

            fclose(fid);
            % fclose('all');

end;

clearvars -except BatchInputFolder InputFiles NumInputFiles Rmode Gmode Bmode Mode_THETA Mode_PHI Mode_RADIUS

pause(2);
    









