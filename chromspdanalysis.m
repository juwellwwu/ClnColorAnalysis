function[ClnMaskPixelArea,UniqueClnMaskPixelArea_OrAndMin_MaskCt] = chromspdanalysis(UniqueClnNameList,Batch_CvxHull_VerticeFace_CellArray,BatchCvxHullClonalMaskOutputFolder,THETAchoice,PHIchoice,CorrTHETAbinMatrix_FileNameString,CorrTHETAbinEdge_FileNameString,CorrPHIbinMatrix_FileNameString,CorrPHIbinEdge_FileNameString,THETACTRxcorrgrid_FileNameString,PHICTRycorrgrid_FileNameString,IntensityRangeCorrOutputSize,TargetImgCorrTHETASpan,TargetImgCorrPHISpan,TargetImgCorrTHETAPxWidth,TargetImgCorrPHIPxHeight)

%% =====DESCRIPTION=====

% Create clonal masks of multiple clones for clonal assignment on partial THETA-PHI grid
% Calculate chromatic stability for single clones

% ==Usage: 
% Call from function "chromspdanalysis_batch.m"

% ==Output files:  "*CloneMask.tif'"
% Save clone masks for clonal assignments. One for each %nmax chromatic
% spead for each clone
% Stored in sub-directory "ClnColorDEMO/cmodspd output/Isofrac*/CloneMasks"

% ==Output files: "*ORChromSpdArea.tif'"
% Save union of clone masks (of specific %nmax chromatic spread) of each unique single in working directory
% Stored in sub-directory "Working_directory_name CloneMasks" 

% ==Output files: "*ANDChromSpdArea.tif'"
% Save intersection of clone masks (of specific %nmax chromatic spread) of each unique single in working directory
% Stored in sub-directory "Working_directory_name CloneMasks" 

% ==Output files: "*MINChromSpdArea.tif'"
% Save smallest clone masks (of specific %nmax chromatic spread) of each unique single in working directory
% Stored in sub-directory "Working_directory_name CloneMasks"


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% Convert chromatic spreads to transformed THETA-PHI coordinates

ClnMaskPixelArea=zeros(length(UniqueClnNameList),1);

UniqueClnMaskOr_Stack=zeros(TargetImgCorrPHIPxHeight,TargetImgCorrTHETAPxWidth,length(UniqueClnNameList));
UniqueClnMaskAnd_Stack=ones(TargetImgCorrPHIPxHeight,TargetImgCorrTHETAPxWidth,length(UniqueClnNameList));
UniqueClnMaskMin_Stack=ones(TargetImgCorrPHIPxHeight,TargetImgCorrTHETAPxWidth,length(UniqueClnNameList));
UniqueClnMaskPixelArea_OrAndMin_MaskCt=zeros(length(UniqueClnNameList),4);

for j=1:size(Batch_CvxHull_VerticeFace_CellArray,2)
    
            XFP_RGB=zeros(size(Batch_CvxHull_VerticeFace_CellArray{j}.THETAPHIvertices,1),3);
            [XFP_RGB(:,3),XFP_RGB(:,2),XFP_RGB(:,1)]=sph2cart(Batch_CvxHull_VerticeFace_CellArray{j}.THETAPHIvertices(:,1),Batch_CvxHull_VerticeFace_CellArray{j}.THETAPHIvertices(:,2),1.00000);
            DataPtCt=size(XFP_RGB,1);

            sphTHETAchoice=THETAchoice;
            sphPHIchoice=PHIchoice;

            sphTHETAgrid=0:sphTHETAchoice:90;
            sphPHIgrid=0:sphPHIchoice:90;
            
            OrgPHIDataPtHistInd=zeros(size(DataPtCt,1),'uint32');
            OrgTHETADataPtHistInd=zeros(size(DataPtCt,1),'uint32');
            CorrPHIDataPtHistInd=zeros(size(DataPtCt,1),'uint32');
            CorrTHETADataPtHistInd=zeros(size(DataPtCt,1),'uint32');

            [OrgPHIDataPtCt,OrgPHIDataPtHistInd]=histc(Batch_CvxHull_VerticeFace_CellArray{j}.THETAPHIvertices(:,2),sphPHIgrid-PHIchoice/2);
            [OrgTHETADataPtCt,OrgTHETADataPtHistInd]=histc(Batch_CvxHull_VerticeFace_CellArray{j}.THETAPHIvertices(:,1),sphTHETAgrid-THETAchoice/2);

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

            parfor k=1:size(CorrPHIDataPtHistInd,1)
                CorrPHIDataPtbinEdge(k,1)=dlmread(CorrPHIbinEdge_FileNameString,'\t',[CorrPHIDataPtHistInd(k,1)-1 0 CorrPHIDataPtHistInd(k,1)-1 0]);
            end;

            fclose('all');

            CorrPHITHETAEdge{1}=CorrPHIbinEdge(1:end); 

            ClnCorrPHImin=CorrPHIbinEdge(2,1);
            ClnCorrPHImax=CorrPHIbinEdge(end-1,1);

            clearvars CorrPHIbinEdge;

            CorrTHETAbinEdge=zeros([maxCorrTHETADataPtHistInd-minCorrTHETADataPtHistInd+2,1]);
            CorrTHETAbinEdge_FileID=fopen(CorrTHETAbinEdge_FileNameString);
            CorrTHETAbinEdge=dlmread(CorrTHETAbinEdge_FileNameString,'\t',[((minCorrTHETADataPtHistInd-1)-1) 0 ((maxCorrTHETADataPtHistInd+1)-1) 0]);

            parfor k=1:size(CorrTHETADataPtHistInd,1)
                CorrTHETADataPtbinEdge(k,1)=dlmread(CorrTHETAbinEdge_FileNameString,'\t',[CorrTHETADataPtHistInd(k,1)-1 0 CorrTHETADataPtHistInd(k,1)-1 0]);
            end;

            fclose('all');

            CorrPHITHETAEdge{2}=CorrTHETAbinEdge(1:end);

            ClnCorrTHETAmin=CorrTHETAbinEdge(2,1);
            ClnCorrTHETAmax=CorrTHETAbinEdge(end-1,1);

            clearvars CorrTHETAbinEdge;

            CorrTHETADataPtbinEdge_ImgShift=round((CorrTHETADataPtbinEdge+TargetImgCorrTHETASpan/2)*TargetImgCorrTHETAPxWidth/(TargetImgCorrTHETASpan))';
            CorrPHIDataPtbinEdge_ImgShift=TargetImgCorrPHIPxHeight-round((CorrPHIDataPtbinEdge+TargetImgCorrPHISpan/2)*TargetImgCorrPHIPxHeight/(TargetImgCorrPHISpan))';

            [ImgMeshCorrTHETA,ImgMeshCorrPHI]=meshgrid([1:1:TargetImgCorrTHETAPxWidth],[1:1:TargetImgCorrPHIPxHeight]);
            CloneMask=inpolygon(ImgMeshCorrTHETA,ImgMeshCorrPHI,CorrTHETADataPtbinEdge_ImgShift,CorrPHIDataPtbinEdge_ImgShift);
            
            UniqueClnMaskOr_Stack(:,:,Batch_CvxHull_VerticeFace_CellArray{j}.UniqueClnID)=bsxfun(@or,UniqueClnMaskOr_Stack(:,:,Batch_CvxHull_VerticeFace_CellArray{j}.UniqueClnID),CloneMask);
            UniqueClnMaskAnd_Stack(:,:,Batch_CvxHull_VerticeFace_CellArray{j}.UniqueClnID)=bsxfun(@and,UniqueClnMaskAnd_Stack(:,:,Batch_CvxHull_VerticeFace_CellArray{j}.UniqueClnID),CloneMask);
                
            ClnMaskPixelArea(j,1)=numel(find(CloneMask));
            if numel(find(UniqueClnMaskMin_Stack(:,:,Batch_CvxHull_VerticeFace_CellArray{j}.UniqueClnID))) > ClnMaskPixelArea(j,1)
                UniqueClnMaskMin_Stack(:,:,Batch_CvxHull_VerticeFace_CellArray{j}.UniqueClnID) = CloneMask;
            end;
            
            UniqueClnMaskPixelArea_OrAndMin_MaskCt(Batch_CvxHull_VerticeFace_CellArray{j}.UniqueClnID,4)=UniqueClnMaskPixelArea_OrAndMin_MaskCt(Batch_CvxHull_VerticeFace_CellArray{j}.UniqueClnID,4)+1;
            
            % Save chromatic spreads on transformed THETA-PHI grid
            imwrite(CloneMask,strcat(BatchCvxHullClonalMaskOutputFolder,'/',Batch_CvxHull_VerticeFace_CellArray{j}.FileName,' CloneMask.tif'), 'Compression','none');
            
            fclose('all');
            close all hidden 
            clearvars -except ClnMaskPixelArea UniqueClnMaskPixelArea_OrAndMin_MaskCt UniqueClnMaskOr_Stack UniqueClnMaskAnd_Stack UniqueClnMaskMin_Stack UniqueClnNameList Batch_CvxHull_VerticeFace_CellArray BatchCvxHullClonalMaskOutputFolder THETAchoice PHIchoice CorrTHETAbinMatrix_FileNameString CorrTHETAbinEdge_FileNameString CorrPHIbinMatrix_FileNameString CorrPHIbinEdge_FileNameString THETACTRxcorrgrid_FileNameString PHICTRycorrgrid_FileNameString IntensityRangeCorrOutputSize Fig1Query Fig3Query Fig5Query TargetImgCorrTHETASpan TargetImgCorrPHISpan TargetImgCorrTHETAPxWidth TargetImgCorrPHIPxHeight;  
            
end;


% Chromatic spread calc (for chromatic stability)
UniqueClnMaskOr_Stack(UniqueClnMaskOr_Stack>0)=1;

for jj=1:length(UniqueClnNameList)
    imwrite(UniqueClnMaskOr_Stack(:,:,jj),char(strcat(BatchCvxHullClonalMaskOutputFolder,'/',UniqueClnNameList(jj),' ORChromSpdArea.tif')), 'Compression', 'none');
    imwrite(UniqueClnMaskAnd_Stack(:,:,jj),char(strcat(BatchCvxHullClonalMaskOutputFolder,'/',UniqueClnNameList(jj),' ANDChromSpdArea.tif')), 'Compression', 'none');
    imwrite(UniqueClnMaskMin_Stack(:,:,jj),char(strcat(BatchCvxHullClonalMaskOutputFolder,'/',UniqueClnNameList(jj),' MINChromSpdArea.tif')), 'Compression', 'none');
    UniqueClnMaskPixelArea_OrAndMin_MaskCt(jj,1)=numel(find(UniqueClnMaskOr_Stack(:,:,jj)));
    UniqueClnMaskPixelArea_OrAndMin_MaskCt(jj,2)=numel(find(UniqueClnMaskAnd_Stack(:,:,jj)));
    UniqueClnMaskPixelArea_OrAndMin_MaskCt(jj,3)=numel(find(UniqueClnMaskMin_Stack(:,:,jj)));
    
end;


 close all hidden; fclose('all'); clearvars -except ClnMaskPixelArea UniqueClnMaskPixelArea_OrAndMin_MaskCt