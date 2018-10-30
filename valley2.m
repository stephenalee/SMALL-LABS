function PhaseMask=valley2(f,params)
% clear all, close all;clc
% Update History:

% 6/27/2013 YL: Added the "autofill" function. Now it is OK to use lower
% values for "higher_thresh" and/or higher values for "lower_thresh"
% without worrying too much about over-segmentation (i.e., cell masks
% having holes inside). Also, the connectivity criterion for the binary
% labelling has been changed from 8(default) to 4 in the autofill section.

% 6/172013 YL: Updated Comments.

% 5/3/2016 DR: parameters are now fields of a structure

%% ------------------------------------------------------------------------
%  User-Defined Parameters:
%  ------------------------------------------------------------------------

% Change the to zero only if you still have under-segmentation regardless
% of threshold.
dilate_factor=params.nDilation;

% Can be negative, but usually use [-1 0.5].
lower_thresh=params.lowThresh;

% Usually use positive values [0 3] , should be larger than
% lower_threshold. *Decrease* this value first when cells are
% UNDER-segmented (e.g., 2 or more cells grouped as 1 region).
higher_thresh=params.highThresh;% .001;

autofill=params.autofill;

% Note: If cells are under-segmented, increase the lower_thresh and/or
% decrease the higher_thresh (recommended); If cells are over-segmented,
% decrease the lower_thresh and increase the higher_thresh.

% minimum allowable area (in px) for a single segmented region. Regions
% with area smaller than this will be discarded.
min_area=params.minArea;

% maximum allowable area (in px) for a single segmented region. Regions
% with area larger than this will be discarded.
max_area=params.maxArea;

%%

f=uint16(f);
f=imcomplement(f); % Invert intensity

% Laplacian of Gaussian filtering
[g,~]=edge(f,'log', 0);

f = imfilter(f, fspecial('average', 5), 'replicate');
% figure, imshow(f, []); title('average- or gaussian- filtered image');

%-------------------------------------------------------------------------%
% Define Valley filters
V = zeros(3,3); V(2,2) = -1;
A1 = V; A1(3,1) = 1;
A2 = V; A2(2,1) = 1;
A3 = V; A3(1,1) = 1;
B1 = V; B1(3,2) = 1;
B3 = V; B3(1,2) = 1;
C1 = V; C1(3,3) = 1;
C2 = V; C2(2,3) = 1;
C3 = V; C3(1,3) = 1;
%-------------------------------------------------------------------------%

A1 = imfilter(f, A1, 'corr', 'replicate', 'same');
C3 = imfilter(f, C3, 'corr', 'replicate', 'same');
A2 = imfilter(f, A2, 'corr', 'replicate', 'same');
C2 = imfilter(f, C2, 'corr', 'replicate', 'same');
A3 = imfilter(f, A3, 'corr', 'replicate', 'same');
C1 = imfilter(f, C1, 'corr', 'replicate', 'same');
B3 = imfilter(f, B3, 'corr', 'replicate', 'same');
B1 = imfilter(f, B1, 'corr', 'replicate', 'same');

V1 = min(A1, C3); V2 = min(A2, C2); V3 = min(A3, C1); V4 = min(B3, B1);
V1 = max(V1, V2); V2 = max(V3, V4); V = max(V1, V2);

% figure, imshow(V,[]), title('valley values')
% figure, imshow(V~=0,[]), title('non-zero valley values')

% 8-connectivity correlation
mean_V = mean(V(:)); std_V = std2(V);
T = [mean_V+lower_thresh*std_V, mean_V + higher_thresh*std_V];
V2 = uint8(V>= T(2))*10; V2(V2~=10)=1; % Strong threshold
V1 = uint8(V>= T(1)); % Weak threshold
V1 = imfilter(V1.*V2, [1 1 1; 1 1 1; 1 1 1], 'corr', 'replicate', 'same')>=10;

% figure, imshow(V2,[]); title('strong threshold')
% figure, imshow(~V1,[]); title('strong and weak threshold')

% Thresholding
f2 = im2bw(f, graythresh(f));
% figure, imshow(f2, []); title('thresholding')

f3 = f2 & ~V1 & ~g;
% figure, imshow(f3, []); title('threshold and edge combined')

f4 = bwmorph(f3,'close',inf);
% figure, imshow(f4,[]); title('close')

f5 = imfill(f4,'holes');
% figure, imshow(f5,[]); title('imfill')

f6 = bwmorph(f5, 'majority', inf');
% figure,imshow(f6,[]); title('majority')

f7 = bwlabel(f6, 4);
% figure, imshow(f7, []); title('bwlabel')

seg_area = regionprops(f7, 'area');

f7(ismember(f7,find([seg_area.Area]<min_area|[seg_area.Area]>max_area)))=0;

f7 = bwmorph(f7, 'dilate', dilate_factor);
% figure, imshow(f7, []); title('final dilation by 1')
% f7 = bwlabel(f7 ~= 0, 4);
% figure, imshow(f7, []); title('bwlabel without small regions and dilated by 1')
% % Reassign region label after getting rid of small regions

f7 = bwlabel(~bwmorph(~f7,'diag',5));
% figure, imshow(f7, []); title('Unfilled Phase Mask')

PhaseMask = f7;

if autofill == 1
    PhaseMask = bwconvhull(PhaseMask~=0,'objects',4);
    PhaseMask = bwlabel(PhaseMask, 4);
    % figure, imshow(PhaseMask, []); title('Filled Phase Mask')
end
end