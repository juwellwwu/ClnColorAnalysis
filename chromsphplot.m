function [] = chromsphplot(FileNameString,chromsphplotOutputFolder,XFP_CloneColorBGR,UfmAreaCorrFactorMatrix,DataPtCtChoice,RchannelChoice,GchannelChoice,BchannelChoice,CloneColorAssignChannelChoice,sphTHETAchoice,sphPHIchoice,SphMeshGrid,AFTieringQuery,AFMultiplierLoQuery,AFMultiplierHiQuery)

 %% =====DESCRIPTION=====

% Plot spherical scatter plot and spherical histogram of one population

% ==Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output files: "*SphScatterPlot.tif'
% Plot and save spherical scatter plots

% ==Output files: "*SphHistogram.tif'"
% Plot and save spherical histogram


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUT

% 'y' to plot spherical scatter plot
Fig2Query='y';

% Approx # data pts to plot in Fig2. -1 to plot all
DataPtPlotFig=2.5E4; 

% 'y' to plot spherical histogram
Fig4Query='y';

% 1 to plot Fig4 cts in linear scale, 2 in log10 scale
Fig4Lin1Log2Query=1;


%% Load Color Data

[~,FileName,FileExt]=fileparts(FileNameString);

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
ColorChoice(4)=CloneColorAssignChannelChoice;

fprintf(strcat('R set to:\t',num2str(ColorChoice(1)),'\n'));
fprintf(strcat('G set to:\t',num2str(ColorChoice(2)),'\n'));
fprintf(strcat('B set to:\t',num2str(ColorChoice(3)),'\n'));
fprintf(strcat('CloneID set to:\t',num2str(ColorChoice(4)),'\n\n\n'));

XFP_RGB=zeros(DataPtCt,4,'double');

parfor i=1:4
    if ColorChoice(i)<9999 
        XFP_RGB(:,i)=dlmread(FileNameString,'\t',[1,ColorChoice(i)-1,DataPtCt,ColorChoice(i)-1]);
    else
        XFP_RGB(:,i)=zeros(DataPtCt,1,'double');
    end;
end;

for i=1:3
    if ColorChoice(i)==9999 
        HeaderMtx{1}{ColorChoice(i)}='None Assigned ';
    end;
end;


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
        % FileNameSpec=strcat(FileName,' R',num2str(ColorChoice(1)),'G',num2str(ColorChoice(2)),'B',num2str(ColorChoice(3)),' T',num2str(sphTHETAchoice),' P',num2str(sphPHIchoice),' RCBLo',num2str(AFMultiplierLoQuery),'Hi',num2str(AFMultiplierHiQuery));
        FileNameSpec=strcat(FileName,' T',num2str(sphTHETAchoice),'P',num2str(sphPHIchoice),' RCBLo',num2str(AFMultiplierLoQuery),'Hi',num2str(AFMultiplierHiQuery));
else
        % FileNameSpec=strcat(FileName,' R',num2str(ColorChoice(1)),'G',num2str(ColorChoice(2)),'B',num2str(ColorChoice(3)),' T',num2str(sphTHETAchoice),'P',num2str(sphPHIchoice));
        FileNameSpec=strcat(FileName,' T',num2str(sphTHETAchoice),'P',num2str(sphPHIchoice));
end;


%% Convert to spherical coordinates

sphedge{1}=0:sphTHETAchoice:90;
sphedge{2}=0:sphPHIchoice:90;

XFP_THETArad=zeros([size(DataPtCt),1],'double');
XFP_PHIrad=zeros([size(DataPtCt),1],'double');
XFP_RADIUS=zeros([size(DataPtCt),1],'double');
XFP_SorttempRGB=zeros([size(DataPtCt),6],'double');
XFP_SortIndex=zeros([size(DataPtCt),1],'double');
XFP_THETAPHI=zeros([size(DataPtCt),2],'double');

[XFP_THETArad,XFP_PHIrad,XFP_RADIUS]=cart2sph(XFP_RGB(:,3),XFP_RGB(:,2),XFP_RGB(:,1));

[XFP_SorttempRGB,XFPSortIndex]=sortrows([XFP_RGB(:,1),XFP_RGB(:,2),XFP_RGB(:,3),XFP_THETArad,XFP_PHIrad,XFP_RADIUS,XFP_RGB(:,4)],-6);

XFP_RGB=[XFP_SorttempRGB(:,1),XFP_SorttempRGB(:,2),XFP_SorttempRGB(:,3)];
XFP_THETAPHI=[180/pi().*(XFP_SorttempRGB(:,4)),180/pi().*(XFP_SorttempRGB(:,5))];
XFP_RADIUS=[XFP_SorttempRGB(:,6)];

if CloneColorAssignChannelChoice<9999;
    XFP_RGB_CloneColorAssign=[XFP_SorttempRGB(:,7)];
else
    XFP_RGB_CloneColorAssign=[];
end;
    
NormXFP_RGB=bsxfun(@rdivide,XFP_RGB,XFP_RADIUS);

clear XFP_RGB_THETA_PHI_RADIUS_SortIndex XFP_SorttempRGB XFP_THETArad XFP_PHIrad XFP_THETA XFP_PHI XFPSortIndex SortFileHeaderRow SortFileNameString;


%% == FIGURE 2: SPHERICAL SCATTER PLOT

if Fig2Query=='y';
            
            scrsz = get(0,'ScreenSize');
            scrszWscalefactor=1200/scrsz(3);
            scrszHscalefactor=900/scrsz(4);

            figHandle2=figure('Position',[25 scrsz(4)-(scrsz(4)*scrszHscalefactor)-100 scrsz(3)*scrszWscalefactor scrsz(4)*scrszHscalefactor],'Color','w','PaperPositionMode','auto');
            
            hold all
            grid off;
            axis on;
            axis([0 1 0 1 0 1]);

            set(gca,'xtick',[0:0.25:1],'ytick',[0:0.25:1],'ztick',[0:0.25:1]);
            set(gca,'XColor',[0,0,0.5],'YColor',[0,0.5,0],'ZColor',[0.5,0,0]);
            
            xtickhandle=get(gca,'xtick');
            ytickhandle=get(gca,'ytick');
            ztickhandle=get(gca,'ztick');

            axis square; 
            
            set(gca,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on');
            set(gca,'XMinorTick','on','YMinorTick','on','ZMinorTick','on');

            view([135 45]);
            camtarget([0.5,0.5,0.5]);

            if numel(SphMeshGrid(:,:,1))<8300;                              
                    SphPtSize=1;                   
                    SphPtColor=horzcat(reshape(SphMeshGrid(:,:,3),numel(SphMeshGrid(:,:,1)),1),reshape(SphMeshGrid(:,:,2),numel(SphMeshGrid(:,:,1)),1),reshape(SphMeshGrid(:,:,1),numel(SphMeshGrid(:,:,1)),1));
                     
                    Fig2SphGridHandle=scatter3(reshape(SphMeshGrid(:,:,1),numel(SphMeshGrid(:,:,1)),1),reshape(SphMeshGrid(:,:,2),numel(SphMeshGrid(:,:,1)),1),reshape(SphMeshGrid(:,:,3),numel(SphMeshGrid(:,:,1)),1),SphPtSize,SphPtColor,'o','filled');
            end;

            XFPDataPtSize=10;

            if isempty(XFP_RGB_CloneColorAssign)==1
                XFP_NegClearRGB=(XFP_RGB>(0-realmin)).*XFP_RGB;
                XFPDataPtColor=bsxfun(@rdivide,XFP_NegClearRGB,XFP_RADIUS); 
            else
                XFPDataPtColor=fliplr(XFP_CloneColorBGR(XFP_RGB_CloneColorAssign,:));
            end;
            
            clear XFP_NegClearRGB;
            
            PrintDataPtSkip=max(1,round(DataPtCt/DataPtPlotFig));
            
            Fig2NormXFPPtHandle=scatter3(NormXFP_RGB(1:PrintDataPtSkip:DataPtCt,3),NormXFP_RGB(1:PrintDataPtSkip:DataPtCt,2),NormXFP_RGB(1:PrintDataPtSkip:DataPtCt,1),XFPDataPtSize,XFPDataPtColor(1:PrintDataPtSkip:DataPtCt,:),'filled');
            
            xlabelhandle=xlabel(strcat('Blue: ',HeaderMtx{1}{ColorChoice(3)}));
            ylabelhandle=ylabel(strcat('Green: ',HeaderMtx{1}{ColorChoice(2)}));
            zlabelhandle=zlabel(strcat('Red: ',HeaderMtx{1}{ColorChoice(1)}));
            set(xlabelhandle,'FontName','Arial','FontSize',10);
            set(ylabelhandle,'FontName','Arial','FontSize',10);
            set(zlabelhandle,'FontName','Arial','FontSize',10);

            title({FileNameSpec;strcat('# Data pts printed/total: ',num2str(size(XFP_RGB(1:PrintDataPtSkip:DataPtCt,3),1)),'/',num2str(DataPtCtAll))});
            titlehandle=get(gca,'title');
            set(titlehandle, 'FontSize', 6);

            hold off
            
            print(figHandle2,'-dtiffn','-r150',strcat(chromsphplotOutputFolder,'/',FileNameSpec,' SphScatterPlot.tif'));
            
            close(figHandle2);
            pause(1);
           
 end;

clear XFP_RGB XFPDataPtColor XFPDataPtSize XFP_RADIUS NormXFPDataPtSize SphPtColor SphPtSize TextDataLabels; 
 

%% == FIGURE 4: Spherical Histogram

if Fig4Query=='y'
    
            XFP_NearestSphGridPtSearchIndex=zeros(DataPtCt,1,'double');
            XFP_NearestSphGridPtSearchIndex=knnsearch([reshape(SphMeshGrid(:,:,1),numel(SphMeshGrid(:,:,1)),1),reshape(SphMeshGrid(:,:,2),numel(SphMeshGrid(:,:,1)),1),reshape(SphMeshGrid(:,:,3),numel(SphMeshGrid(:,:,1)),1)],fliplr(NormXFP_RGB));
            
            NormXFP_RGB_mean=[mean(NormXFP_RGB(:,1)),mean(NormXFP_RGB(:,2)),mean(NormXFP_RGB(:,3))];
            
            clear NormXFP_RGB
            
            SphGridBinCt=zeros(numel(SphMeshGrid(:,:,1)),1,'double');
            SphGridBinCt=histc(XFP_NearestSphGridPtSearchIndex,[1:1:numel(SphMeshGrid(:,:,1))]);
            
            clear XFP_NearestSphGridPtSearchIndex; 
            
            if isempty(UfmAreaCorrFactorMatrix)<1             
                    SphGridBinCt=SphGridBinCt./reshape(UfmAreaCorrFactorMatrix,numel(SphMeshGrid(:,:,1)),1);                 
            end;
  

            figHandle4=figure('Position',[50 scrsz(4)-800-100 1200 900],'Color','w','PaperPositionMode','auto');
            set(gcf,'Renderer','OpenGL');

            hold all
            grid on;
            axis on;
            
            axis([0 1 0 1 0 1]);
            
            set(gca,'xtick',[0:0.5:1],'ytick',[0:0.5:1],'ztick',[0:0.5:1],'XColor',[0,0,0.5],'YColor',[0,0.5,0],'ZColor',[0.5,0,0]);
            
            set(gca,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on');
            set(gca,'XMinorTick','on','YMinorTick','on','ZMinorTick','on');

            xtickhandle=get(gca,'xtick');
            ytickhandle=get(gca,'ytick');
            ztickhandle=get(gca,'ztick');

            axis square;

            view([135 45]);
            camtarget([0.5,0.5,0.5]);
            
            if Fig4Lin1Log2Query==1
                    set(gca,'Clim',[0,max(SphGridBinCt(:))]);
            elseif Fig4Lin1Log2Query==2
                    set(gca,'Clim',[0,log10(max(SphGridBinCt(:)))]);
            end;
            
            TargetColor=[];
            TargetColorMap=[];
            
            if CloneColorAssignChannelChoice<9999 && ~isempty(XFP_CloneColorBGR);
                if numel(unique(XFP_RGB_CloneColorAssign))==1 
                    TargetColor=fliplr(XFP_CloneColorBGR(XFP_RGB_CloneColorAssign(1,1),1:3));
                else
                    TargetColor=NormXFP_RGB_mean;
                end;
            else
                    TargetColor=NormXFP_RGB_mean;                
            end;
            
            TargetColorMapRes=100;                      
            TargetColorMap=zeros(3,TargetColorMapRes+1);
            TargetColorMap=horzcat(([TargetColorMapRes:-(TargetColorMapRes-(1-TargetColor(1,1))*TargetColorMapRes)/TargetColorMapRes:(1-TargetColor(1,1))*TargetColorMapRes]/TargetColorMapRes)',([TargetColorMapRes:-(TargetColorMapRes-(1-TargetColor(1,2))*TargetColorMapRes)/TargetColorMapRes:(1-TargetColor(1,2))*TargetColorMapRes]/TargetColorMapRes)',([TargetColorMapRes:-(TargetColorMapRes-(1-TargetColor(1,3))*TargetColorMapRes)/TargetColorMapRes:(1-TargetColor(1,3))*TargetColorMapRes]/TargetColorMapRes)');                       
            TargetColorMap=bsxfun(@minus,1,TargetColorMap);                    
            colormap(TargetColorMap);
                  
            shading flat;
            colorbar('EastOutside');
            
             if Fig4Lin1Log2Query==1
                    Fig4SurfDensityHandle=surf(SphMeshGrid(:,:,1),SphMeshGrid(:,:,2),SphMeshGrid(:,:,3),reshape(SphGridBinCt,size(SphMeshGrid(:,:,1),1),size(SphMeshGrid(:,:,1),2)),'EdgeAlpha',0,'FaceAlpha',1);
             elseif Fig4Lin1Log2Query==2
                    Fig4SurfDensityHandle=surf(SphMeshGrid(:,:,1),SphMeshGrid(:,:,2),SphMeshGrid(:,:,3),reshape(log10(SphGridBinCt),size(SphMeshGrid(:,:,1),1),size(SphMeshGrid(:,:,1),2)),'EdgeAlpha',0,'FaceAlpha',1);
             end;
             
            xlabelhandle=xlabel(strcat('Blue: ',HeaderMtx{1}{ColorChoice(3)}));
            ylabelhandle=ylabel(strcat('Green: ',HeaderMtx{1}{ColorChoice(2)}));
            zlabelhandle=zlabel(strcat('Red: ',HeaderMtx{1}{ColorChoice(1)}));
            set(xlabelhandle,'FontName','Arial','FontSize',10);
            set(ylabelhandle,'FontName','Arial','FontSize',10);
            set(zlabelhandle,'FontName','Arial','FontSize',10);

            title({strcat(FileNameSpec,'; # Data pts: ',num2str(DataPtCt),'/',num2str(DataPtCtAll),',MaxCt: ',num2str(max(SphGridBinCt(:))),' (',num2str(max(SphGridBinCt(:))/DataPtCt*100),'%)')});
            titlehandle=get(gca,'title');
            set(titlehandle, 'FontSize', 6);

            hold off
            
            if Fig4Lin1Log2Query==1
                print(figHandle4,'-dtiffn','-r150',strcat(chromsphplotOutputFolder,'/',FileNameSpec,' SphHistogram.tif'));
            elseif Fig4Lin1Log2Query==2
                print(figHandle4,'-dtiffn','-r150',strcat(chromsphplotOutputFolder,'/',FileNameSpec,' LOG10 SphHistogram.tif'));
            end;
                        
            close(figHandle4);
            pause(1);

end;


fprintf('Computation Complete.\n\n\n');

close all hidden;  fclose('all'); 

clear -except chromsphplotOutputFolder XFP_CloneColorBGR UfmAreaCorrFactorMatrix DataPtCtChoice RchannelChoice GchannelChoice BchannelChoice CloneColorAssignChannelChoice sphTHETAchoice sphPHIchoice SphMeshGrid AFTieringQuery AFMultiplierLoQuery AFMultiplierHiQuery




