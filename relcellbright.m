function [] = relcellbright(AF_IsoSurf_VerticeFaceArray,FileNameString,OutputFolderString,DataPtCtChoice,RchannelChoice,GchannelChoice,BchannelChoice)

%% =====DESCRIPTION=====

% Determine relative cell brightness of datapts in one txt file

% ==Usage: 
% User specifies variables in "USER INPUT" section.
% User specifies output file location for "*RelCellBright Histgm.tif", "*RelCellBright.txt"

% ==Output file: "*RelCellBrightHistgm.tif'"
% Plot histogram of cell ct at each relative cell brightness bin

% ==Output file: "*RelCellBrightHistCt.txt"
% Save data for "*RelCellBright Histgm.tif'"

% ==Output file: "*RelCellBright.txt"
% Relative cell brightness info amended to color data file
% Input for script "cmodespread_batch.m" (single clone chromatic mode, chromatic spread analysis)
% Input for script "dataptct_batch.m" (single clone relative clonal brightness analysis)
% Input for script "cloneassign_batch.m" (polyclonal population clonal assignment)


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% USER INPUT

% 'y' to plot relative cell brightness histogram
Fig2Query='y';

% 'y' to save cell count in each relative cell brighness bin
AFMultiplierHistCtSave_Query='y';


%% Load Cell Color Data

[FilePath,FileName,FileExt]=fileparts(FileNameString);

FileID=fopen(FileNameString);
HeaderRow=fgets(FileID);
HeaderMtx=textscan(HeaderRow,'%s','delimiter','\t');

DataPtCt=0;
while (fgets(FileID) ~= -1),
  DataPtCt=DataPtCt+1;
end

DataPtCtAll=DataPtCt;

if DataPtCtChoice<0
    DataPtCtChoice=DataPtCt;
else
    DataPtCt=DataPtCtChoice;
end;

fprintf('\n');


ColorChoice(1)=RchannelChoice;
ColorChoice(2)=GchannelChoice;
ColorChoice(3)=BchannelChoice;

fprintf(strcat(FileNameString,'\n'));
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

fclose(FileID);
fclose('all');

FileName=strcat(FileName,' R',num2str(ColorChoice(1)),'G',num2str(ColorChoice(2)),'B',num2str(ColorChoice(3)));

if ~exist(OutputFolderString, 'dir')
            mkdir(OutputFolderString);
end;


%% Color Data Pt wrt AFMode

AFMode0_XFP_LinRGB=zeros(DataPtCt,4);
AFMode0_XFP_LinRGB=horzcat([bsxfun(@minus,XFP_RGB,fliplr(AF_IsoSurf_VerticeFaceArray.LinBGRMode)),[1:1:DataPtCt]']);
[AFMode0_XFP_LinRGB_THETArad,AFMode0_XFP_LinRGB_PHIrad,~] = cart2sph(AFMode0_XFP_LinRGB(:,3),AFMode0_XFP_LinRGB(:,2),AFMode0_XFP_LinRGB(:,1));

Flag_LessThanAFConvexHull=any(bsxfun(@lt,AFMode0_XFP_LinRGB(:,1:3),[min(AF_IsoSurf_VerticeFaceArray.vertices(:,3)),min(AF_IsoSurf_VerticeFaceArray.vertices(:,2)),min(AF_IsoSurf_VerticeFaceArray.vertices(:,1))]),2);


%% XFP_RGB Relative Cell Brightness in xAF Unit

AFMode0_XFP_LinRGB_AFMultiplier=zeros(DataPtCt,1);
XFP_AFIsoSurfConvHull_Query=zeros(DataPtCt,1);

XFP_AFIsoSurfConvHull_TierCt=20000;
DivisorStep=10;
AFScaleCLimMax=-9999;

DataPtCt_QueryTemp=DataPtCt;

waitbarhandle=waitbar(0,'Begin calculating relative cell brightness (in xAF units)...');

for i=1:XFP_AFIsoSurfConvHull_TierCt
    
    Divisor_i=DivisorStep+DivisorStep*(i-1);
    
    waitbar(i/XFP_AFIsoSurfConvHull_TierCt,waitbarhandle,sprintf('Testing %.0f x Autofluorescence. %.0f cells to be processed...',Divisor_i,DataPtCt_QueryTemp))

    XFP_AFIsoSurfConvHull_Query=inhull([AFMode0_XFP_LinRGB(:,3)/Divisor_i,AFMode0_XFP_LinRGB(:,2)/Divisor_i,AFMode0_XFP_LinRGB(:,1)/Divisor_i],[AF_IsoSurf_VerticeFaceArray.vertices(:,1),AF_IsoSurf_VerticeFaceArray.vertices(:,2),AF_IsoSurf_VerticeFaceArray.vertices(:,3)],AF_IsoSurf_VerticeFaceArray.index);
 
    AFMode0_XFP_LinRGB_AFMultiplier(AFMode0_XFP_LinRGB(find(XFP_AFIsoSurfConvHull_Query),4),1)=Divisor_i;

    AFMode0_XFP_LinRGB(find(XFP_AFIsoSurfConvHull_Query),:)=[];
  
    DataPtCt_QueryTemp=size( AFMode0_XFP_LinRGB,1);
    XFP_AFIsoSurfConvHull_Query=zeros(DataPtCt_QueryTemp,1);
    
    if nnz(AFMode0_XFP_LinRGB_AFMultiplier)>0.99*DataPtCt && AFScaleCLimMax==-9999
            AFScaleCLimMax=Divisor_i;
    end;
    
    if nnz(AFMode0_XFP_LinRGB_AFMultiplier)==DataPtCt
            AFMode0_XFP_LinRGB_AFMultiplier(:,i+1:end)=[];
            break;  
    end;
    
end;

close(waitbarhandle);

clearvars AFCentroid0_XFP_LinRGB AFCentroid0_XFP_LinRGB_RADIUS ;


%% Relative Cell Brightness Histogram: Plot and Save

AFMultiplier_HistogramCt=hist(AFMode0_XFP_LinRGB_AFMultiplier,[0:DivisorStep:max(AFMode0_XFP_LinRGB_AFMultiplier)]);

% Plot

if Fig2Query=='y'

            scrszWscalefactor=1000/1280;
            scrszHscalefactor=500/778;

            scrsz = get(0,'ScreenSize');
            figHandle2=figure('Position',[100 scrsz(4)-(scrsz(4)*scrszHscalefactor)-100 scrsz(3)*scrszWscalefactor scrsz(4)*scrszHscalefactor],'Color','w','PaperPositionMode', 'auto');

            set(gcf, 'Renderer', 'OpenGL');

            hold on;
            axis on;
            grid on;

            axis([0 max(AFMode0_XFP_LinRGB_AFMultiplier) 0 max(AFMultiplier_HistogramCt)]);

            xlabel(strcat('x Autofluorescence' ),'FontSize',10);
            ylabel(strcat('# XFP Data Pts' ),'FontSize',10);

            title({FileName;strcat(num2str(DataPtCt),' Data Pts, RelCellBright(xAF) Histogram');strcat('AF File: ',AF_IsoSurf_VerticeFaceArray.FileName)},'FontSize',10);
            titlehandle=get(gca,'title');
            set(titlehandle, 'FontSize', 8);

            Fig2HistogramHandle=bar([0:DivisorStep:max(AFMode0_XFP_LinRGB_AFMultiplier)],AFMultiplier_HistogramCt);

            export_fig(figHandle2,strcat(OutputFolderString,'/',FileName,' RelCellBrightHistgm.tif'),'-nocrop','-a1');

            hold off;

            close(figHandle2);

end;

% Save
if AFMultiplierHistCtSave_Query=='y'
    
            FileNameString_LinAFTHistCt=strcat(OutputFolderString,'/',FileName,' RelCellBrightHistCt.txt');

            LinAFTHistCt_HeaderMtx{1}=strcat('RelCellBright(xAF) HistBinCenter');
            LinAFTHistCt_HeaderMtx{2}=strcat('RelCellBright(xAF) HistBinCenter CellCt');

            FileID=fopen(FileNameString_LinAFTHistCt, 'w');

            for i=1:numel(LinAFTHistCt_HeaderMtx);
                fprintf(FileID,'%s\t',LinAFTHistCt_HeaderMtx{i});        
            end;

            fprintf(FileID,'\n');

            dlmwrite(FileNameString_LinAFTHistCt,horzcat([0:DivisorStep:max(AFMode0_XFP_LinRGB_AFMultiplier)]',AFMultiplier_HistogramCt'),'delimiter','\t','-append');

            fclose(FileID);
            fclose('all');

            clear FileNameString_LinAFTHistCt LinAFTHistCt_HeaderMtx;
    
end;


%% Save Relative Cell Brightness (xAF) info

% Append to original Color Data File
if DataPtCt==DataPtCtAll
    
    FileNameString_TierIndexAppend=strcat(OutputFolderString,'/',FileName,' RelCellBright.txt');

    HeaderMtx{1}{numel(HeaderMtx{1})+1}=strcat(AF_IsoSurf_VerticeFaceArray.FileName,' RelCellBright (xAF)');
    HeaderMtx{1}{numel(HeaderMtx{1})+1}=strcat('FLAG 1 if < 1xAF');
    HeaderMtx{1}{numel(HeaderMtx{1})+1}=strcat('THETA(deg) wrt AFMode');
    HeaderMtx{1}{numel(HeaderMtx{1})+1}=strcat('PHI(deg) wrt AFMode');

    FileID=fopen(FileNameString);
    FCSData=dlmread(FileNameString,'\t',1,0);
    fclose(FileID);
    fclose('all');

    FCSData=horzcat(FCSData,AFMode0_XFP_LinRGB_AFMultiplier,Flag_LessThanAFConvexHull,180/pi().*([AFMode0_XFP_LinRGB_THETArad,AFMode0_XFP_LinRGB_PHIrad]));

    FileID=fopen(FileNameString_TierIndexAppend, 'w');

    HeaderMtx_TierIndexAppend=HeaderMtx{1}; 

    for i=1:numel(HeaderMtx_TierIndexAppend);
        fprintf(FileID,'%s\t',HeaderMtx_TierIndexAppend{i});
    end
    fprintf(FileID,'\n');

    dlmwrite(FileNameString_TierIndexAppend,FCSData,'delimiter','\t','-append');

    fclose(FileID);
    fclose('all');

    clearvars FCS_Data AFMode0_XFP_LinRGB_THETArad AFMode0_XFP_LinRGB_PHIrad;    
    
end;


end

