 clc; close all; clear all; clear java;
 pause(2);

%% =====DESCRIPTION=====

% Calculate and save transformed THETA-PHI grid properties

% ==Usage: 
% User follows prompt for input

% ==Output files: "*tfmTHETA-PHIgrid.tif"
% Plot transformed THETA-PHI grid @ chosen THETA-PHI grid step

% ==Output files: transformed THETA-PHI grid property files:
% "*THETACTRxcorrgrid.txt"
% "*PHICTRycorrgrid.txt"
% "*IntensityRangeCorrFactorMatrix.txt"
% "*IntensityRangeCorrTHETAbinMatrix.txt"
% "*IntensityRangeCorrPHIbinMatrix.txt"
% "*IntensityRangeCorrOutputSize.txt"
% Save transformed THETA-PHI grid property files
% Stored in sub-directory "ClnColorDEMO/T*P*Step*tfmgrid"
% "cloneassign_batch.m"


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUT

fprintf('Choose THETA-PHI grid step:\n') 
fprintf('1: 0.01 deg for transformed THETA-PHI grid calculations (need 5.3GB space).\n') 
fprintf('2: 0.2 deg for adjusted spherical histograms (13.3MB).\n') 
fprintf('3: 5.0 deg for demo.\n') 
fprintf('Press "Enter" when done.\n\n') 

THETAPHIchoice=input('Grid step: ');

Fig3Query='';
 
if THETAPHIchoice==1
        Fig3Query='n';
        THETAchoice=0.01;
        PHIchoice=0.01;
        CorrTHETAedgeStepMultiple=20;
        CorrPHIedgeStepMultiple=20;
        if ~exist('ClnColorDEMO/T0.01P0.01Step20tfmgrid', 'dir')
            mkdir('ClnColorDEMO/T0.01P0.01Step20tfmgrid');
        end;
        Outputfoldername='ClnColorDEMO/T0.01P0.01Step20tfmgrid';
elseif THETAPHIchoice==2
        Fig3Query='n';
        THETAchoice=0.2;
        PHIchoice=0.2;
        CorrTHETAedgeStepMultiple=1;
        CorrPHIedgeStepMultiple=1;
        if ~exist('ClnColorDEMO/T0.2P0.2Step1tfmgrid', 'dir')
            mkdir('ClnColorDEMO/T0.2P0.2Step1tfmgrid');
        end;
        Outputfoldername='ClnColorDEMO/T0.2P0.2Step1tfmgrid';
elseif THETAPHIchoice==3
        Fig3Query='y';
        THETAchoice=5;
        PHIchoice=5;
        CorrTHETAedgeStepMultiple=1;
        CorrPHIedgeStepMultiple=1;
        if ~exist('ClnColorDEMO/T5.0P5.0Step1tfmgrid', 'dir')
            mkdir('ClnColorDEMO/T5.0P5.0Step1tfmgrid');
        end;
        Outputfoldername='ClnColorDEMO/T5.0P5.0Step1tfmgrid';
end;
 

 %% CALCULATIONS
 
THETAedge=0:THETAchoice:90; 
PHIedge=0:PHIchoice:90; 

THETACTRedge=(0+THETAchoice/2):THETAchoice:(90-THETAchoice/2); 
PHICTRedge=(0+PHIchoice/2):PHIchoice:(90-PHIchoice/2);

PHICTRVolSphElmtRatio=zeros(size(PHICTRedge,2),1);
ProductPHITHETABinCtrCART=zeros(size(PHIedge,2)-1,size(THETAedge,2)-1);
BinDeltaProduct_CART_XYZ=zeros(size(PHIedge,2)-1,size(THETAedge,2)-1);
THETAPHIgrid=zeros(size(THETACTRedge,2)*size(PHICTRedge,2),2);

Func_SphereVolElement=@(PHIvar) sin(pi/2-PHIvar);

 for i=1:size(PHICTRedge,2)
            PHICTRVolSphElmtRatio(i,1)=quadl(Func_SphereVolElement,(PHIedge(1,i)*pi/180),(PHIedge(1,i+1)*pi/180),1E-20);
 end;
 
[meshTHETACTRedge,meshPHICTRedge]=meshgrid(THETACTRedge*pi/180,PHICTRedge*pi/180);
ProductPHITHETABinCtrCART=cos(meshPHICTRedge).*cos(meshTHETACTRedge).*cos(meshPHICTRedge).*sin(meshTHETACTRedge).*sin(meshPHICTRedge);

BinDeltaProduct_CART_XYZ=(repmat(PHICTRVolSphElmtRatio,1,numel(THETACTRedge)))./(ProductPHITHETABinCtrCART);

THETAPHIgrid=[reshape(repmat(THETACTRedge,size(PHICTRedge,2),1),numel(THETACTRedge)*numel(PHICTRedge),1),repmat(PHICTRedge,1, size(THETACTRedge,2))'];

BinDeltaProduct_CART_XYZ=BinDeltaProduct_CART_XYZ./min(BinDeltaProduct_CART_XYZ(:));

THETACTRxcorrfactor=zeros(size(BinDeltaProduct_CART_XYZ));
PHICTRycorrfactor=zeros(size(BinDeltaProduct_CART_XYZ));
PHICTRycorrfactorCol=zeros(size(BinDeltaProduct_CART_XYZ,1),1);
 
THETACTRxcorrgrid=zeros(size(BinDeltaProduct_CART_XYZ));
PHICTRycorrgrid=zeros(size(BinDeltaProduct_CART_XYZ));
  
PHICTRycorrfactorCol(:,1)=sqrt(min(BinDeltaProduct_CART_XYZ,[],2));
PHICTRycorrfactorCol=PHICTRycorrfactorCol./min(PHICTRycorrfactorCol(:));
PHICTRycorrfactor=repmat(PHICTRycorrfactorCol,1,size(BinDeltaProduct_CART_XYZ,2));
 
THETACTRxcorrfactor=BinDeltaProduct_CART_XYZ./PHICTRycorrfactor;

CorrTHETAedgeStep=min(THETACTRxcorrfactor(:))*CorrTHETAedgeStepMultiple;
CorrPHIedgeStep=min(PHICTRycorrfactor(:))*CorrPHIedgeStepMultiple;
 
SumCol_THETACTRxcorrfactor=sum(THETACTRxcorrfactor,2);

minSumCol_THETACTRxcorrfactor=min(SumCol_THETACTRxcorrfactor(:));

THETAxcorrGridSize=max(SumCol_THETACTRxcorrfactor(:,1));

SumRow_PHICTRycorrfactor=sum(PHICTRycorrfactor,1);

minSumRow_PHICTRycorrfactor=min(SumRow_PHICTRycorrfactor(:));

PHIycorrGridSize=max(SumRow_PHICTRycorrfactor(1,:));

THETACTRxcorrgrid=cumsum(THETACTRxcorrfactor,2)+repmat((THETAxcorrGridSize/2-SumCol_THETACTRxcorrfactor/2),1,size(THETACTRxcorrfactor,2));

PHICTRycorrgrid=cumsum(PHICTRycorrfactor,1)+repmat((PHIycorrGridSize/2-SumRow_PHICTRycorrfactor/2),size(PHICTRycorrfactor,1),1);

PHICTRycorrgrid=PHICTRycorrgrid-max(PHICTRycorrgrid(:))/2;
THETACTRxcorrgrid=THETACTRxcorrgrid-(THETACTRxcorrgrid(1,floor(size(THETACTRxcorrgrid,2)/2))+THETACTRxcorrgrid(1,ceil(size(THETACTRxcorrgrid,2)/2)))/2;

clearvars THETACTRxcorrfactor PHICTRycorrfactor

CorrTHETAedgeEnd=max(abs(min(THETACTRxcorrgrid(:))-CorrTHETAedgeStep/2),abs(max(THETACTRxcorrgrid(:))+CorrTHETAedgeStep/2));

CorrTHETAedgeNeg=(-CorrTHETAedgeStep:-CorrTHETAedgeStep:-ceil(CorrTHETAedgeEnd/CorrTHETAedgeStep)*CorrTHETAedgeStep)';
CorrTHETAedge=[flipud(CorrTHETAedgeNeg);0;-(CorrTHETAedgeNeg)];

CorrTHETABinCtrHistIndMatrix=zeros(size(THETACTRxcorrgrid));

CorrTHETABinCtrHistIndTemp=zeros(size(THETACTRxcorrgrid,2),1);

for i=1:size(PHICTRycorrgrid,1)
    [CorrTHETABinCtrCt,CorrTHETABinCtrHistIndTemp]=histc(THETACTRxcorrgrid(i,:),CorrTHETAedge);
    CorrTHETABinCtrHistIndMatrix(i,:)=CorrTHETABinCtrHistIndTemp'; 
end;

clearvars CorrTHETABinCtrCt CorrTHETABinCtrHistIndTemp CorrTHETAedgeLef;

CorrPHIedgeEnd=max(abs(min(PHICTRycorrgrid(:))-CorrPHIedgeStep/2),abs(max(PHICTRycorrgrid(:))+CorrPHIedgeStep/2));

CorrPHIedgeNeg=(-CorrPHIedgeStep:-CorrPHIedgeStep:-ceil(CorrPHIedgeEnd/CorrPHIedgeStep)*CorrPHIedgeStep)';
CorrPHIedge=[flipud(CorrPHIedgeNeg);0;-(CorrPHIedgeNeg)];

CorrPHIBinCtrHistIndMatrix=zeros(size(PHICTRycorrgrid));

CorrPHIBinCtrHistIndTemp=zeros(size(PHICTRycorrgrid,2),1);

for i=1:size(THETACTRxcorrgrid,2)
    [CorrPHIBinCtrCt,CorrPHIBinCtrHistIndTemp]=histc(PHICTRycorrgrid(i,:),CorrPHIedge);
    CorrPHIBinCtrHistIndMatrix(i,:)=CorrPHIBinCtrHistIndTemp; 
end;

clearvars CorrPHIBinCtrCt CorrPHIBinCtrHistIndTemp CorrPHIedgeLeft;


%% FIGURE 3: Transformed THETA-PHI grid

if Fig3Query=='y'

    scrsz = get(0,'ScreenSize');

    figHandle3=figure('Position',[5 scrsz(4)-900-105 3360 900],'Color','w','PaperPositionMode','auto');

    subplot(1,1,1,'Position',[0.05 0.1 0.9, 0.8])
    hold all
    grid off

    axis on
    axis([min(THETACTRxcorrgrid(:)) max(THETACTRxcorrgrid(:)) min(PHICTRycorrgrid(:)) max(PHICTRycorrgrid(:))]);

    shading flat;
    colormap('jet');

    C=ones(size(THETACTRxcorrgrid));

    Fig3Handle=pcolor(THETACTRxcorrgrid,PHICTRycorrgrid,C);
    set(Fig3Handle,'FaceAlpha',0.05,'EdgeAlpha',1);

    xlabelhandle=xlabel(strcat('Tfm THETA'));
    ylabelhandle=ylabel(strcat('Tfm PHI'));

    set(xlabelhandle,'FontName','Arial','FontSize',8);
    set(ylabelhandle,'FontName','Arial','FontSize',8);
    
    set(gca,'TickDir','out');

    title({strcat('Transformed THETA-PHI grid     THETAstep: ',num2str(THETAchoice),'deg; PHIstep: ',num2str(PHIchoice),'deg')});
    titlehandle=get(gca,'title');
    set(titlehandle, 'FontSize', 7);

    hold off;

    print(figHandle3,'-dtiffn','-r72',strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' tfmTHETA-PHIgrid.tif'));
    
    close(figHandle3);
    
end;


%% Export Transformed THETA-PHI grid properties

%
THETACTRxcorrgrid_FileNameString=strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' THETACTRxcorrgrid.txt');

fid=fopen(THETACTRxcorrgrid_FileNameString,'w');

dlmwrite(THETACTRxcorrgrid_FileNameString,THETACTRxcorrgrid,'delimiter','\t','-append','precision',16);

fclose(fid);
fclose('all');

THETACTRxcorrgrid_sizeDim1=size(THETACTRxcorrgrid,1);
THETACTRxcorrgrid_sizeDim2=size(THETACTRxcorrgrid,2);

clearvars THETACTRxcorrgrid;


% 
PHICTRycorrgrid_FileNameString=strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' PHICTRycorrgrid.txt');

fid=fopen(PHICTRycorrgrid_FileNameString,'w');

dlmwrite(PHICTRycorrgrid_FileNameString,PHICTRycorrgrid,'delimiter','\t','-append','precision',16);

fclose(fid);
fclose('all');

PHICTRycorrgrid_sizeDim1=size(PHICTRycorrgrid,1);
PHICTRycorrgrid_sizeDim2=size(PHICTRycorrgrid,2);

clearvars PHICTRycorrgrid


% 
BinDeltaProduct_CART_XYZ_FileNameString=strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' IntensityRangeCorrFactorMatrix.txt');

fid=fopen(BinDeltaProduct_CART_XYZ_FileNameString,'w');

dlmwrite(BinDeltaProduct_CART_XYZ_FileNameString,BinDeltaProduct_CART_XYZ,'delimiter','\t','-append','precision',16);

fclose(fid);
fclose('all');

BinDeltaProduct_CART_XYZ_sizeDim1=size(BinDeltaProduct_CART_XYZ,1);
BinDeltaProduct_CART_XYZ_sizeDim2=size(BinDeltaProduct_CART_XYZ,2);

clearvars BinDeltaProduct_CART_XYZ


% 
CorrTHETABinCtrHistIndMatrix_FileNameString=strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' IntensityRangeCorrTHETAbinMatrix.txt');

fid=fopen(CorrTHETABinCtrHistIndMatrix_FileNameString,'w');

dlmwrite(CorrTHETABinCtrHistIndMatrix_FileNameString,CorrTHETABinCtrHistIndMatrix,'delimiter','\t','-append','precision',16);

fclose(fid);
fclose('all');

CorrTHETABinCtrHistIndMatrix_sizeDim1=size(CorrTHETABinCtrHistIndMatrix,1);
CorrTHETABinCtrHistIndMatrix_sizeDim2=size(CorrTHETABinCtrHistIndMatrix,2);

clearvars CorrTHETABinCtrHistIndMatrix 


% 
CorrPHIBinCtrHistIndMatrix_FileNameString=strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' IntensityRangeCorrPHIbinMatrix.txt');

fid=fopen(CorrPHIBinCtrHistIndMatrix_FileNameString,'w');

dlmwrite(CorrPHIBinCtrHistIndMatrix_FileNameString,CorrPHIBinCtrHistIndMatrix,'delimiter','\t','-append','precision',16);

fclose(fid);
fclose('all');

CorrPHIBinCtrHistIndMatrix_sizeDim1=size(CorrPHIBinCtrHistIndMatrix,1);
CorrPHIBinCtrHistIndMatrix_sizeDim2=size(CorrPHIBinCtrHistIndMatrix,2);

clearvars CorrPHIBinCtrHistIndMatrix 


%
CorrTHETAedge_FileNameString=strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' IntensityRangeCorrTHETAbinEdge.txt');

fid=fopen(CorrTHETAedge_FileNameString,'w');

dlmwrite(CorrTHETAedge_FileNameString,CorrTHETAedge,'delimiter','\t','-append','precision',16);

fclose(fid);
fclose('all');

CorrTHETAedge_sizeDim1=size(CorrTHETAedge,1);
CorrTHETAedge_sizeDim2=size(CorrTHETAedge,2);

clearvars CorrTHETAedge;


%
CorrPHIedge_FileNameString=strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' IntensityRangeCorrPHIbinEdge.txt');

fid=fopen(CorrPHIedge_FileNameString,'w');

dlmwrite(CorrPHIedge_FileNameString,CorrPHIedge,'delimiter','\t','-append','precision',16);

fclose(fid);
fclose('all');

CorrPHIedge_sizeDim1=size(CorrPHIedge,1);
CorrPHIedge_sizeDim2=size(CorrPHIedge,2);

clearvars CorrPHIedge;


%
CorrTHETAPHIFilesSize_FileNameString=strcat(Outputfoldername,'/THETA',num2str(THETAchoice),'PHI',num2str(PHIchoice),' IntensityRangeCorrOutputSize.txt');

CorrTHETAPHIFilesSize_HeaderRow={strcat('THETACTRxcorrgrid_sizeDim1');strcat('THETACTRxcorrgrid_sizeDim2');strcat('PHICTRycorrgrid_sizeDim1');strcat('PHICTRycorrgrid_sizeDim2');
                                                                    strcat('CorrFactorMatrix_sizeDim1');strcat('CorrFactorMatrix_sizeDim2');strcat('CorrTHETAbinMatrix_sizeDim1');strcat('CorrTHETAbinMatrix_sizeDim2');
                                                                    strcat('CorrPHIbinMatrix_sizeDim1');strcat('CorrPHIbinMatrix_sizeDim2');strcat('CorrTHETAedge_sizeDim1');strcat('CorrTHETAedge_sizeDim2');
                                                                    strcat('CorrPHIedge_sizeDim1');strcat('CorrPHIedge_sizeDim2');strcat('CorrTHETAminSpanStepping');strcat('CorrPHIStepping');
                                                                    strcat('CorrTHETA45degLength');strcat('CorrAllPHILength')};

fid=fopen(CorrTHETAPHIFilesSize_FileNameString,'w');

for i=1:numel(CorrTHETAPHIFilesSize_HeaderRow);
    fprintf(fid,'%s\t',CorrTHETAPHIFilesSize_HeaderRow{i});
end

fprintf(fid,'\n');

CorrTHETAPHIFilesSizeData=[THETACTRxcorrgrid_sizeDim1,THETACTRxcorrgrid_sizeDim2,PHICTRycorrgrid_sizeDim1,PHICTRycorrgrid_sizeDim2,BinDeltaProduct_CART_XYZ_sizeDim1,BinDeltaProduct_CART_XYZ_sizeDim2,CorrTHETABinCtrHistIndMatrix_sizeDim1,CorrTHETABinCtrHistIndMatrix_sizeDim2,CorrPHIBinCtrHistIndMatrix_sizeDim1,CorrPHIBinCtrHistIndMatrix_sizeDim2,CorrTHETAedge_sizeDim1,CorrTHETAedge_sizeDim2,CorrPHIedge_sizeDim1,CorrPHIedge_sizeDim2,CorrTHETAedgeStepMultiple,CorrPHIedgeStepMultiple,minSumCol_THETACTRxcorrfactor,minSumRow_PHICTRycorrfactor];
dlmwrite(CorrTHETAPHIFilesSize_FileNameString,CorrTHETAPHIFilesSizeData,'delimiter','\t','-append','precision',16);

fclose(fid);

fclose('all');

clear all; clear java;
