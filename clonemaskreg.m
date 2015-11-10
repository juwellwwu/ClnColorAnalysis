function[]=clonemaskreg(ClonalMask_FileNameString,XRegGridCoordinate,YRegGridCoordinate,RegImg_Width,RegImg_Height)

%% =====DESCRIPTION=====

% Register clone mask using *inverse_transf RAW.txt" from bUnwarpJ

% ==Usage: 
% User specifies variables in "USER INPUT" section.
% User specifies location of tfmTHETA-PHI grid property files
% (in section "%% Transformed THETA-PHI Grid Info"

% ==Output files: "Working_directory/Fsk*/Result/*CloneCellCount.tif'"
% Save cell ct assigned to each clone for each polyclonal population


%%  =====DO NOT REMOVE=====

% Supplementary software code for Wu et al. "Defining Clonal Color in Fluorescent Multi-Clonal Tracking"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: Nov-2015


%% Warp clone mask using bUnwarpJ data

ClonalMask_PreRegImg=zeros(RegImg_Height,RegImg_Width);
ClonalMask_RegShiftImg=zeros(RegImg_Height,RegImg_Width);

ClonalMask_FileNameString=char(ClonalMask_FileNameString);
[ClonalMask_FilePath,ClonalMask_FileName,ClonalMask_FileExt]=fileparts(ClonalMask_FileNameString);
ClonalMask_FileName=regexprep(ClonalMask_FileName,' CloneMask', '');

ClonalMask_PreRegImg=imread(ClonalMask_FileNameString);

% Determine the register-shifted coordinates of non-zero (in clone)
% clonal mask points. Find gives LINEAR indices of non-zero values

RegClonalMask_YX=zeros(size(find(ClonalMask_PreRegImg),1),2); 

RegClonalMask_YX(:,1)=round(YRegGridCoordinate(find(ClonalMask_PreRegImg)));
RegClonalMask_YX(:,2)=round(XRegGridCoordinate(find(ClonalMask_PreRegImg)));

clearvars ClonalMask_PreRegImg;


% Make sure X coordinates stay within limits
RegClonalMask_OutofBoundsIndex=bsxfun(@gt,RegClonalMask_YX(:,2),RegImg_Width);
RegClonalMask_YX(RegClonalMask_OutofBoundsIndex,2)=RegImg_Width;

RegClonalMask_OutofBoundsIndex=bsxfun(@lt,RegClonalMask_YX(:,2),1);
RegClonalMask_YX(RegClonalMask_OutofBoundsIndex,2)=1;

% Make sure Y coordinates stay within limits
RegClonalMask_OutofBoundsIndex=bsxfun(@gt,RegClonalMask_YX(:,1),RegImg_Height);
RegClonalMask_YX(RegClonalMask_OutofBoundsIndex,1)=RegImg_Height;

RegClonalMask_OutofBoundsIndex=bsxfun(@lt,RegClonalMask_YX(:,1),1);
RegClonalMask_YX(RegClonalMask_OutofBoundsIndex,1)=1;

% Write "1" values at register shifted positions.  
% Y-axis values are row subscripts, X-axis values are col.
ClonalMask_RegShiftImg=accumarray([RegClonalMask_YX(:,1),RegClonalMask_YX(:,2)],1,[RegImg_Height,RegImg_Width]);

% Dilate silghtly to close out any small holes as result of grid stretching
% Image Dilation: Overlap region given a dilution of 2x2 diamond
 ClonalMask_RegShiftImg_DiluteStructureElement=strel('diamond',2);
 ClonalMask_RegShiftImg=imdilate(ClonalMask_RegShiftImg,ClonalMask_RegShiftImg_DiluteStructureElement);

% bwboundaries gives (X,Y) coordinates of boundary of mask. Mask needs
% inversion. cell2mat combines coordinates of all "objects"
RegShiftClonalMask_BoundaryYX=bwboundaries(ClonalMask_RegShiftImg,'noholes');
RegShiftClonalMask_BoundaryYX=cell2mat(RegShiftClonalMask_BoundaryYX);

% Meshgrid for inpolygon later. 
RegImgMeshX=zeros(RegImg_Height,RegImg_Width);
RegImgMeshY=zeros(RegImg_Height,RegImg_Width);
[RegImgMeshX,RegImgMeshY]=meshgrid([1:1:RegImg_Width],[1:1:RegImg_Height]);

ClonalMask_RegShiftImg=inpolygon(RegImgMeshY,RegImgMeshX,RegShiftClonalMask_BoundaryYX(:,1),RegShiftClonalMask_BoundaryYX(:,2));

mkdir(ClonalMask_FilePath,'REG');
imwrite(ClonalMask_RegShiftImg,strcat(ClonalMask_FilePath,'/REG/',ClonalMask_FileName,' RegCloneMask.tif'), 'Compression','none');
pause(1);

clearvars RegImgMeshX RegImgMeshY RegShiftClonalMask_BoundaryYX RegClonalMask_YX ClonalMask_RegShiftImg
