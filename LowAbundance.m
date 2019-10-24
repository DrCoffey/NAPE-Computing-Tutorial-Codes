%% RNA Scope Image analysis for Low Abundance Genes
clear all
close all
fname="C:\Users\Kevin\Desktop\NAPE Computing\Low Abundance Raw\2.jpg";

%% Read Image
im = (imread(fname));
[~,imName] = fileparts(fname);

%% Original Unedited Figure
figure('Color','k')
f = imshow(im);

%% Extract Blue Channel for Dapi Segmentation
imgblue = im(:,:,3); % blue, nuclei only signal
f = imshow(imgblue);

% Subtract the uneven background from the image
se = strel('disk',60);
imgblue = imtophat(imgblue,se);
f = imshow(imgblue);

% Binarize Image
bw_blue = imbinarize(imgblue);
f = imshow(bw_blue);

% Find Connected Aread of White
bw_blue = bwareaopen(bw_blue, 50);
f = imshow(bw_blue);

% % Remove Clumping
bw_blue = bwareafilt(bw_blue,[50 2500]);
f = imshow(bw_blue);

% Look for connected components
cc_blue = bwconncomp(bw_blue, 8);

% Create Dapi Centroids
radius = 30;
s = regionprops(cc_blue,'centroid');
centroids = cat(1, s.Centroid);

% Remove centroids touching the border (hard to look inside the partial circles)
NotAtBorder = centroids(:,1) < size(im,2) - (radius + 1) &...
centroids(:,1) > radius + 1 &...
centroids(:,2) < size(im,1) - (radius + 1) &...
centroids(:,2) > radius + 1;
centroids = centroids(NotAtBorder,:);  
cell_radious(1:length(centroids),1)=30;

%Show Cell Identification
f = imshow(im)
hold on
viscircles(centroids,cell_radious,'LineStyle','-','LineWidth',.1,'Color','cyan');


%% Extract Green Channel for Signal Segmentation
imggreen = im(:,:,2);
f=imshow(imggreen);

% Subtract the uneven background from the image
imggreen = imtophat(imggreen,se);
figure;
f=imshow(imggreen);

% Binarize image 
bw_green = imbinarize(imggreen,.15);
f=imshow(bw_green);

% filter clumps
bw_green = bwareafilt(bw_green,[1 50]);
f=imshow(bw_green);

% Look for connected components
cc_green = bwconncomp(bw_green, 8);
s = regionprops(cc_green,'centroid'); 
green_centroids = cat(1, s.Centroid);

%Show Cell Identification
f = imshow(im)
hold on
% scatter(centroids(:,1),centroids(:,2), 'c','Marker','.')
scatter(green_centroids(:,1),green_centroids(:,2), 'g','Marker','.')


%% Extract Red Channel for Signal Segmentation
imgred = im(:,:,1);
f = imshow(imgred)

% Subtract the uneven background from the image
imgred = imtophat(imgred,se);
f=imshow(imgred);

% Binarize image 
bw_red = imbinarize(imgred,.15);
f=imshow(bw_red);

% filter clumps
bw_red = bwareafilt(bw_red,[1 50]);
f=imshow(bw_red);

% imshow(bw_red);
cc_red = bwconncomp(bw_red, 8);
s = regionprops(cc_red,'centroid');
red_centroids = cat(1, s.Centroid);

%Show Cell Identification
f = imshow(im)
hold on
% scatter(centroids(:,1),centroids(:,2), 'c','Marker','.')
scatter(red_centroids(:,1),red_centroids(:,2), 'r','Marker','.')

%% Final Analasis Figure
figure;
imshow(im)
hold on
scatter(green_centroids(:,1),green_centroids(:,2), 'g','Marker','.')
scatter(red_centroids(:,1),red_centroids(:,2), 'r','Marker','.')
viscircles(centroids,cell_radious,'LineStyle','-','LineWidth',.1,'Color','cyan');

%% Create fake BW images with 1 white pixel at the centroids for cell quant
a=size(bw_red);
red_dot_matrix=zeros(a(1),a(2));
for y=1:length(red_centroids(:,1))
red_dot_matrix(floor(red_centroids(y,2)),floor(red_centroids(y,1)))=1;
end
green_dot_matrix=zeros(a(1),a(2));
for y=1:length(green_centroids(:,1))
green_dot_matrix(floor(green_centroids(y,2)),floor(green_centroids(y,1)))=1;
end

%% Green & Red per Cell (clever trick to use the strel function to
% approximate a circle)
SE = strel('disk',radius,0);
figure;
f=imshow(SE.Neighborhood);

%% create indicies of the circle in reference to the centroid
circle = find(SE.Neighborhood);
[I,J] = ind2sub(size(SE.Neighborhood),circle);
I = I - radius - 1;
J = J - radius - 1;
dist = sqrt(J.^2 + I.^2);
[dist,distIdx] = sort(dist,'ascend');
I = I(distIdx);
J = J(distIdx);

%% Create index of each cell within the full image
for k=1:length(centroids)
% Index of cell    
cellIndex = sub2ind(size(imgred),round(centroids(k,2)) + I,round(centroids(k,1)) + J);

% Amount of red/green dots per cell
red_count(k,1) = sum(red_dot_matrix(cellIndex));
green_count(k,1)=sum(green_dot_matrix(cellIndex));

% Total RNA Scope Signal/ Cell
redSignal(k,:) = sum(imgred(cellIndex));
greenSignal(k,:)= sum(imggreen(cellIndex));
end

%% Another way Insert Circles
I=insertShape(im,'circle',[centroids(green_count>=1,:) cell_radious(green_count>=1)],'Color','green','LineWidth',2);
I=insertShape(I,'circle',[centroids(red_count>=1,:) cell_radious(red_count>=1)],'Color','red','LineWidth',2);
I=insertShape(I,'circle',[centroids(green_count>=1 & red_count>=1,:) cell_radious(green_count>=1 & red_count>=1)],'Color','yellow','LineWidth',2);
figure;
f=imshow(I);
hold on
scatter(green_centroids(:,1),green_centroids(:,2), 'g','Marker','.')
scatter(red_centroids(:,1),red_centroids(:,2), 'r','Marker','.')


%% Data Section
totalCells=length(centroids);
totalGreenDots=length(green_centroids);
totalRedDots=length(red_centroids);
totalGreenCells=sum(green_count~=0);
totalRedCells=sum(red_count~=0);
totalBothCells=sum(green_count~=0 & red_count~=0);
greenPerGreenCell=mean(green_count(green_count~=0));
redPerRedCell=mean(red_count(red_count~=0));
proportionGreenCells=totalGreenCells/totalCells;
proportionRedCells=totalRedCells/totalCells;
proportionBothCells=totalBothCells/totalCells;
