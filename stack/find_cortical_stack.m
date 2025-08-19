function [cortical] = find_cortical_stack(image_path, growth_plate_data)
% Modified version of find_cortical.m to work with the stack workflow
% This function identifies cortical regions of bone using an image and growth plate data
% 
% Input:
% image_path - Path to the image file (NoDAPI overlay from Registration folder)
% growth_plate_data - Growth plate data structure with top and bottom coordinates
%
% Output:
% cortical - Binary mask of cortical regions

% Load the image
try
    img = imread(image_path);
    disp(['Loaded image: ', image_path]);
catch ME
    error(['Error loading image: ', ME.message]);
end

% Convert to grayscale if it's a color image
if size(img, 3) > 1
    img_gray = rgb2gray(img);
    disp('Converted color image to grayscale');
else
    img_gray = img;
end

% Extract growth plate coordinates from data
if isstruct(growth_plate_data) && isfield(growth_plate_data, 'bt1')
    top_y = growth_plate_data.bt1.top_y;
    bottom_y = growth_plate_data.bt1.bottom_y;
    disp(['Using growth plate coordinates from structure: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
elseif iscell(growth_plate_data) && numel(growth_plate_data) >= 1
    try
        top_y = growth_plate_data{1, 2, 1};
        bottom_y = growth_plate_data{1, 3, 1};
        disp(['Using growth plate coordinates from cell array: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
    catch
        % Use default values based on image height
        [height, ~] = size(img_gray);
        top_y = round(height * 0.25);
        bottom_y = round(height * 0.75);
        disp(['Using default growth plate coordinates: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
    end
else
    % Use default values based on image height
    [height, ~] = size(img_gray);
    top_y = round(height * 0.25);
    bottom_y = round(height * 0.75);
    disp(['Using default growth plate coordinates: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
end

% Apply thresholding to segment the bone
level = graythresh(img_gray);
bone_binary = imbinarize(img_gray, level * 0.5); % Use an even lower threshold to capture more of the bone

% Clean up the binary image
bone_binary = imclose(bone_binary, strel('disk', 10));
bone_binary = imfill(bone_binary, 'holes');
bone_binary = imdilate(bone_binary, strel('disk', 3));

% Display thresholding result
figure('visible', 'off');
subplot(1,2,1); imshow(img_gray); title('Original Grayscale');
subplot(1,2,2); imshow(bone_binary); title('Thresholded');
[img_dir, img_name, ~] = fileparts(image_path);
viz_dir = fullfile(img_dir, '..', 'Analyzed', 'Cortical');
if ~exist(viz_dir, 'dir')
    mkdir(viz_dir);
end
saveas(gcf, fullfile(viz_dir, [img_name, '_thresholding.png']));
disp(['Saved thresholding visualization to: ', fullfile(viz_dir, [img_name, '_thresholding.png'])]);
close(gcf);

% Create ROI mask using growth plate boundaries
[height, width] = size(bone_binary);
roi = bone_binary;
roi(1:top_y-1, :) = false;  % Clear above top boundary
roi(bottom_y+1:end, :) = false;  % Clear below bottom boundary

% Find edges of the bone using erosion
template = bone_binary - imerode(bone_binary, strel('disk', 1));
template(1:top_y-1, :) = false;
template(bottom_y+1:end, :) = false;

% Label connected components
[L, n] = bwlabel(template);

% Try to increase contrast if not enough regions found
if n < 2
    disp('Not enough cortical regions detected, trying with enhanced edge detection...');
    % Try Canny edge detection instead
    edges = edge(img_gray, 'canny');
    edges = edges & roi; % Only keep edges within ROI
    
    % Display edge detection result
    figure('visible', 'off');
    imshow(edges);
    title('Canny Edge Detection');
    saveas(gcf, fullfile(viz_dir, [img_name, '_edges.png']));
    disp(['Saved edge detection visualization to: ', fullfile(viz_dir, [img_name, '_edges.png'])]);
    close(gcf);
    
    % Try with edge detection
    [L, n] = bwlabel(edges);
    
    if n < 2
        warning('Could not find two distinct cortical regions even after edge enhancement.');
        cortical = false(size(roi));
        
        % Create empty visualization for consistency
        figure('visible', 'off', 'Position', [100, 100, width, height*2]);
        subplot(2,2,1); imshow(img); title('Original Image');
        subplot(2,2,2); imshow(bone_binary); title('Bone Binary Mask');
        subplot(2,2,3); imshow(roi); title('ROI');
        subplot(2,2,4); imshow(cortical); title('No Cortical Regions Detected');
        saveas(gcf, fullfile(viz_dir, [img_name, '_cortical_analysis.png']));
        disp(['Saved empty cortical visualization to: ', fullfile(viz_dir, [img_name, '_cortical_analysis.png'])]);
        close(gcf);
        
        % Save empty cortical mask
        cortical_mask = cortical;
        save(fullfile(viz_dir, 'cortical_mask.mat'), 'cortical_mask');
        disp(['Saved empty cortical mask to: ', fullfile(viz_dir, 'cortical_mask.mat')]);
        disp('Cortical region analysis completed with no regions detected.');
        return;
    else
        disp(['Edge detection found ', num2str(n), ' regions.']);
        template = edges;
    end
end

stat = regionprops(L, 'BoundingBox', 'Image', 'Area', 'Centroid');
clear L;

% Find cortical regions on left and right sides
% Instead of finding the two longest regions, find the leftmost and rightmost large regions

% Get areas and x-coordinates of centroids
areas = [stat.Area];
centroids_x = zeros(1, n);
for i = 1:n
    centroids_x(i) = stat(i).Centroid(1);
end

% Filter regions by size (keep only reasonably large ones)
min_area = max(areas) * 0.05; % At least 5% of the largest region
valid_regions = find(areas > min_area);

if length(valid_regions) < 2
    warning('Could not find two distinct cortical regions after size filtering.');
    cortical = false(size(roi));
    return;
end

% Among valid regions, find the leftmost and rightmost ones
[~, left_idx_in_valid] = min(centroids_x(valid_regions));
[~, right_idx_in_valid] = max(centroids_x(valid_regions));

% Get the actual indices
idx1 = valid_regions(left_idx_in_valid);
idx2 = valid_regions(right_idx_in_valid);

% If they're the same region (shouldn't happen but just in case)
if idx1 == idx2
    warning('Could not distinguish left and right cortical regions.');
    cortical = false(size(roi));
    return;
end

% Determine which is left and which is right
if stat(idx1).BoundingBox(1) > stat(idx2).BoundingBox(1)
    left = idx2;
    right = idx1;
else
    left = idx1;
    right = idx2;
end

% Initialize cortical mask
cortical = false(size(roi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   left cortex
temp = false(size(roi));
st_x = round(stat(left).BoundingBox(2));
ed_x = st_x + stat(left).BoundingBox(4) - 1;
st_y = round(stat(left).BoundingBox(1));
ed_y = st_y + stat(left).BoundingBox(3) - 1;
temp(st_x:ed_x, st_y:ed_y) = stat(left).Image;
temp_cortical = temp;

cont = 1;
count = 0;
s(1) = length(find(temp & bone_binary));
while cont
    count = count + 1;
    s(count+1) = length(find(temp(:,1:end-count+1) & bone_binary(:,count:end)));
    if s(count+1) / s(1) < 0.75 || count > 500
        break
    end
    temp_cortical(:,count:end) = temp_cortical(:,count:end) | temp(:,1:end-count+1); 
end
cortical = cortical | (temp_cortical & bone_binary);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   right cortex
temp = false(size(roi));
st_x = round(stat(right).BoundingBox(2));
ed_x = st_x + stat(right).BoundingBox(4) - 1;
st_y = round(stat(right).BoundingBox(1));
ed_y = st_y + stat(right).BoundingBox(3) - 1;
temp(st_x:ed_x, st_y:ed_y) = stat(right).Image;
temp_cortical = temp;
clear s

cont = 1;
count = 0;
s(1) = length(find(temp & bone_binary));
while cont
    count = count + 1;
    s(count+1) = length(find(temp(:,count:end) & bone_binary(:,1:end-count+1)));
    if s(count+1) / s(1) < 0.75 || count > 500
        break
    end
    temp_cortical(:,1:end-count+1) = temp_cortical(:,1:end-count+1) | temp(:,count:end); 
end
cortical = cortical | (temp_cortical & bone_binary);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create visualization to check results
figure('visible', 'off', 'Position', [100, 100, width, height*2]);

% Plot original image
subplot(2, 2, 1);
imshow(img);
title('Original Image', 'FontSize', 12);

% Plot bone binary
subplot(2, 2, 2);
imshow(bone_binary);
title('Bone Binary Mask', 'FontSize', 12);

% Plot ROI mask
subplot(2, 2, 3);
imshow(roi);
hold on;
line([1, width], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
line([1, width], [bottom_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
title('ROI with Growth Plate Boundaries', 'FontSize', 12);

% Plot cortical regions
subplot(2, 2, 4);
imshow(img);
hold on;
h = imshow(cat(3, zeros(size(cortical)), cortical, zeros(size(cortical))));
set(h, 'AlphaData', cortical * 0.5);
title('Cortical Regions', 'FontSize', 12);

% Save the visualization
[img_dir, img_name, ~] = fileparts(image_path);
viz_dir = fullfile(img_dir, '..', 'Analyzed', 'Cortical');
if ~exist(viz_dir, 'dir')
    mkdir(viz_dir);
end
viz_path = fullfile(viz_dir, [img_name, '_cortical_analysis.png']);
saveas(gcf, viz_path);
disp(['Saved cortical region visualization to: ', viz_path]);
close(gcf);

% Save the cortical mask to MAT file
cortical_mask = cortical;
save(fullfile(viz_dir, 'cortical_mask.mat'), 'cortical_mask');
disp(['Saved cortical mask to: ', fullfile(viz_dir, 'cortical_mask.mat')]);

disp('Cortical region analysis completed successfully.');
end
