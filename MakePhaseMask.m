function out=MakePhaseMask(mainfold)
%% Maskes phase mask images for movies for use in SMALL LABS




% PHASEMASK PARAMETERS.
paramsPhase.nDilation = 1;
paramsPhase.lowThresh = 0;
paramsPhase.highThresh = 1;
paramsPhase.autofill = 1;
paramsPhase.minArea = 100;
paramsPhase.maxArea = 1e4;
fieldsPhase = fieldnames(paramsPhase);



%% find all the movie files
[datalist,dataloc,~]=uigetfile([mainfold filesep '*.nd2;*.tif*;*.bin'],...
    'Select the movie files','multiselect','on');

if ~dataloc
    display('no data selected')
    return
end

if ~iscell(datalist); datalist={datalist}; end
datalist = cellfun(@(x)[dataloc, x],datalist,'uniformoutput',false);
[dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);


%% get the phase images
[phaselist,phaselistloc,~]=uigetfile([mainfold filesep...
    '*.nd2;*.tif*;*.bin'],'Select the phase images','multiselect','on');

if ~phaselistloc
    display('no data selected')
    return
end

if ~iscell(phaselist); phaselist={phaselist}; end
phaselist = cellfun(@(x)[phaselistloc, x],phaselist,'uniformoutput',false);
[plocs,pnames,pexts]=cellfun(@fileparts,phaselist,'uniformoutput',false);

%% load the phase images

for ii=1:numel(pnames)
    m=matfile([fullfile(dlocs{ii},dnames{ii}),'_PhaseMask.mat'],'Writable',true);
    
    if strcmp(pexts{ii},'.nd2')
        vidid = bfGetReader(phaselist{ii});
        phaseImg = bfGetPlane(vidid,1);
    elseif strcmp(pexts{ii},'.bin')
        phaseImg=binGetFrames2(phaselist{ii},[]);
    else
        phaseImg=imread(fullfile(plocs{ii},[pnames{ii},pexts{ii}]));
    end
    m.phaseImg=phaseImg;
end

%% construct the masks

for ii=1:numel(dnames)
    m=matfile([fullfile(dlocs{ii},dnames{ii}),'_PhaseMask.mat'],'Writable',true);
    phaseImg=m.phaseImg;
    going = 1;
    oldParams = [paramsPhase.nDilation,paramsPhase.lowThresh,paramsPhase.highThresh,...
        paramsPhase.autofill,paramsPhase.minArea,paramsPhase.maxArea];
    phaseMask = valley2(phaseImg,paramsPhase);
    while going
        subplot(2,1,1)
        imshow(phaseMask)
        subplot(2,1,2)
        imshow(phaseImg,[])
        subplot(2,1,1)
        title('Press 1 to change Params, + to add cell by hand, - to remove,0 to end');
        [~,~,Button]=ginput;
        if Button==49
            % prompt user for parameter changes
            dValues = cellfun(@num2str,num2cell(oldParams),'uniformoutput',0);
            newParams=cellfun(@str2double,...
                inputdlg(fieldsPhase,'Phasemask Parameters',...
                numDlgLines,dValues,opts))';
            paramsPhase.nDilation = newParams(1);
            paramsPhase.lowThresh = newParams(2);
            paramsPhase.highThresh = newParams(3);
            paramsPhase.autofill = newParams(4);
            paramsPhase.minArea = newParams(5);
            paramsPhase.maxArea = newParams(6);
            phaseMask = valley2(phaseImg,paramsPhase);
        elseif Button==45
            subplot(2,1,1)
            [xc,yc]=ginput;
            subCells=diag(phaseMask(round(yc),round(xc)));
            for ii=1:length(subCells)
                phaseMask(phaseMask==subCells(ii))=0;
            end
        elseif any(Button==[43,61])
            subplot(2,1,2)
            going2=1;
            while going2
                AddCell=roipoly;
                AddCell=AddCell.*(max(max(phaseMask))+1);
                phaseMask=phaseMask+AddCell;
                subplot(2,1,1)
                imshow(phaseMask)
                title('Continue to add cells? Y/1 or N/0 then enter')
                [~,~,cont]=ginput;
                if cont==110||cont==48
                    going2=0;
                elseif cont==121||cont==49
                    going2=1;
                    subplot(2,1,2)
                end
            end
        elseif Button==48
            going=0;
        end
    end
    phaseMask = selectCells(phaseMask,phaseImg);
    m.PhaseMask = phaseMask;
    m.paramsPhase = paramsPhase;
end
close all
end
function phaseMask=selectCells(phaseMask,img)
% Let user click on the phase mask of cell images to decide which cell to
% analyze subsequently

%#ok<*AGROW>

h=figure;
c=onCleanup(@()close(h));
vSize = size(phaseMask);

subplot(121);
imshow(img,[]);
subplot(122);
imshow(phaseMask~=0,[]);
title('Hit Tab to enter a user defined region');
hold all

% highlight the outlines of the cells
rProp=regionprops(phaseMask,'Convexhull');
for i=1:numel(rProp)
    plot(rProp(i,1).ConvexHull(:,1),rProp(i,1).ConvexHull(:,2),'c-','linewidth',2)
end

set(gcf,'NextPlot','add');
axes
h1 = title('Click the cells to be analyzed. Enter to Proceed.');
set(gca,'Visible','off');
set(h1,'Visible','on');

% Let the user pick cells with good shapes
keyInput=0; clicksX=[]; clicksY=[];
while keyInput~=121
    [clickY,clickX,keyInput]=ginput(1);
    if keyInput==9
        hold on
        title('Select the Phase region')
        subplot(121)
        phaseMask=roipoly;
        break
    end
    plot(clickY,clickX,'m*','markersize',12);
    clicksX=cat(1,clicksX,clickX);
    clicksY=cat(1,clicksY,clickY);
end

% If the user did click on something
if ~isempty(clicksX)
    pmID = findPmID(phaseMask,vSize,cat(2,clicksX,clicksY));
    phaseMask(~ismember(phaseMask,pmID)) = 0;
end
end
function pmID = findPmID(phaseMask,imSize,points)
pmID = zeros(size(points,1),1);

if numel(pmID)>0
    p1 = points(:,1) + 1 < imSize(1);
    p2 = points(:,2) + 1 < imSize(2);
    p3 = points(:,1)>0|points(:,2)>0;
    
    g1 = p1&p2&p3;
    
    points = ceil(points);
    
    pmID(g1) = phaseMask(sub2ind(imSize,points(g1,1),points(g1,2)));
end
end