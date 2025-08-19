function find_growth_plate_stack_for_dia(dia_img_path)
% Function to detect growth plate from 2_Dia.png for cortical analysis
% This is a specialized version of find_growth_plate_stack that works with a single image
% 
% Input:
% dia_img_path - Full path to the 2_Dia.png image

% Check if input image exists
if ~exist(dia_img_path, 'file')
    error(['Image file not found: ', dia_img_path]);
end

% Load the image
try
    img = imread(dia_img_path);
    disp(['Loaded image: ', dia_img_path]);
catch ME
    error(['Error loading image: ', ME.message]);
end

% Look for mineral image (CCC_E01_hF_FL1_M_s1c2.png) in the same folder
[img_dir, ~, ~] = fileparts(dia_img_path);
mineral_img_path = fullfile(img_dir, 'CCC_E01_hF_FL1_M_s1c2.png');

% Check if mineral image exists
has_mineral_image = false;
if exist(mineral_img_path, 'file')
    try
        mineral_img = imread(mineral_img_path);
        disp(['Found mineral image: ', mineral_img_path]);
        has_mineral_image = true;
    catch ME
        disp(['Warning: Error loading mineral image: ', ME.message]);
        disp('Will use Dia image for all boundaries');
    end
else
    disp('Mineral image not found. Will use Dia image for all boundaries.');
end

% Convert to grayscale if color image
if size(img, 3) > 1
    img_gray = rgb2gray(img);
    disp('Converted color image to grayscale');
else
    img_gray = img;
end

% Get image dimensions
[height, width] = size(img_gray);
disp(['Image dimensions: ', num2str(height), 'x', num2str(width)]);

% Apply thresholding to segment the bone
level = graythresh(img_gray);
bone_binary = imbinarize(img_gray, level * 0.7); % Use a lower threshold to capture more of the bone

% Clean up the binary image
bone_binary = imclose(bone_binary, strel('disk', 5));
bone_binary = imfill(bone_binary, 'holes');

% Find growth plate (estimate top and bottom of growth plate)
% First, find the bone boundaries
col_sums = sum(bone_binary, 2);
col_sums_normalized = col_sums / max(col_sums);

% Use gradient to find transition points that might indicate growth plate
col_grad = gradient(col_sums_normalized);
col_grad_smoothed = movmean(col_grad, height/40);

% Find potential growth plate boundaries based on gradient
[~, max_grad_idx] = findpeaks(col_grad_smoothed, 'MinPeakHeight', 0.002, 'MinPeakDistance', height/10);
[~, min_grad_idx] = findpeaks(-col_grad_smoothed, 'MinPeakHeight', 0.002, 'MinPeakDistance', height/10);

% Combine and sort all transition points
all_transitions = sort([max_grad_idx; min_grad_idx]);

% If we have enough transitions, use them to determine growth plate
if length(all_transitions) >= 2
    % Use the first quarter of transitions for top
    top_candidates = all_transitions(all_transitions < height/2);
    if ~isempty(top_candidates)
        top_y = top_candidates(end);
    else
        top_y = round(height * 0.25); % Default if no good candidates
    end
    
    % Use the last quarter of transitions for bottom
    bottom_candidates = all_transitions(all_transitions > height/2);
    if ~isempty(bottom_candidates)
        bottom_y = bottom_candidates(1);
    else
        bottom_y = round(height * 0.75); % Default if no good candidates
    end
else
    % Use default positions if not enough transitions found
    disp('Not enough transition points found. Using default growth plate positions.');
    top_y = round(height * 0.25);
    bottom_y = round(height * 0.75);
end

disp(['Detected growth plate boundaries: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);

% Initialize default horizontal boundaries
left_x = 1;
right_x = width;

% Use mineral image to determine left and right boundaries if available
if has_mineral_image
    % Convert mineral image to grayscale if needed
    if size(mineral_img, 3) > 1
        mineral_gray = rgb2gray(mineral_img);
    else
        mineral_gray = mineral_img;
    end
    
    % Threshold the mineral image to segment bone
    mineral_level = graythresh(mineral_gray);
    mineral_binary = imbinarize(mineral_gray, mineral_level * 0.7);
    
    % Clean up the binary image
    mineral_binary = imclose(mineral_binary, strel('disk', 5));
    mineral_binary = imfill(mineral_binary, 'holes');
    
    % Find the largest connected component (the main bone)
    [labeled, ~] = bwlabel(mineral_binary);
    stats = regionprops(labeled, 'Area', 'BoundingBox', 'PixelIdxList');
    if ~isempty(stats)
        [~, idx] = max([stats.Area]);
        
        % Create bone mask from largest component
        bone_mask = false(size(mineral_binary));
        bone_mask(stats(idx).PixelIdxList) = true;
        
        % Try to use find_periosteum_4 to get a more precise bone boundary
        try
            disp('Applying find_periosteum_4 to refine bone boundary...');
            [bone_mask_refined, ~] = find_periosteum_4(bone_mask);
            
            % Use the refined mask if it's not empty
            if any(bone_mask_refined(:))
                bone_mask = bone_mask_refined;
                disp('Successfully refined bone boundary using find_periosteum_4');
            else
                disp('Warning: Refined bone mask is empty. Using original mask.');
            end
        catch ME
            disp(['Warning: Error in find_periosteum_4: ', ME.message]);
            disp('Using original bone mask.');
        end
        
        % Get bounding box of the bone mask
        stats_refined = regionprops(bone_mask, 'BoundingBox');
        if ~isempty(stats_refined)
            bbox = stats_refined.BoundingBox;
            left_x = max(1, round(bbox(1)));
            right_x = min(width, left_x + round(bbox(3)) - 1);
            disp(['Detected left/right boundaries from mineral image: left_x=', num2str(left_x), ', right_x=', num2str(right_x)]);
        else
            disp('Could not get bounding box from refined mask. Using original.');
            bbox = stats(idx).BoundingBox;
            left_x = max(1, round(bbox(1)));
            right_x = min(width, left_x + round(bbox(3)) - 1);
        end
    else
        disp('No bone regions found in mineral image. Using full width.');
    end
else
    % Use dia image for left and right boundaries if mineral image not available
    [labeled, ~] = bwlabel(bone_binary);
    stats = regionprops(labeled, 'Area', 'BoundingBox');
    
    if ~isempty(stats)
        [~, idx] = max([stats.Area]);
        bbox = stats(idx).BoundingBox;
        left_x = max(1, round(bbox(1)));
        right_x = min(width, left_x + round(bbox(3)) - 1);
        disp(['Detected left/right boundaries from dia image: left_x=', num2str(left_x), ', right_x=', num2str(right_x)]);
    else
        disp('No bone regions found. Using full width.');
    end
end

% Create a visualization
figure('visible', 'off', 'Position', [100, 100, width, height]);
imshow(img);
hold on;

% Draw growth plate boundaries
line([1, width], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
line([1, width], [bottom_y, bottom_y], 'Color', 'g', 'LineWidth', 2);

% Draw left/right boundaries
line([left_x, left_x], [1, height], 'Color', 'b', 'LineWidth', 2);
line([right_x, right_x], [1, height], 'Color', 'b', 'LineWidth', 2);

% Add labels
text(10, top_y - 20, 'Top Growth Plate', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
text(10, bottom_y + 20, 'Bottom Growth Plate', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
text(left_x + 10, 30, 'Left Boundary', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
text(right_x - 150, 30, 'Right Boundary', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);

title('Growth Plate Detection for Cortical Analysis', 'FontSize', 14);

% Save the visualization
[img_dir, img_name, ~] = fileparts(dia_img_path);
output_dir = fullfile(img_dir, '..', 'CCC_E01_hF_FL1_stack_output', 'Analyzed');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
viz_path = fullfile(output_dir, 'dia_growth_plate_visualization.png');
saveas(gcf, viz_path);
disp(['Saved growth plate visualization to: ', viz_path]);
close(gcf);

% Create ROI visualization on shift2_NoDAPI.jpg
registration_dir = fullfile(img_dir, '..', 'CCC_E01_hF_FL1_stack_output', 'Registration');
if exist(registration_dir, 'dir')
    nodapi_path = fullfile(registration_dir, 'CCC_E01_hF_FL1_shift2_NoDAPI.jpg');
    if exist(nodapi_path, 'file')
        try
            % Load NoDAPI image
            nodapi_img = imread(nodapi_path);
            
            % Create ROI mask
            roi_mask = false(size(img_gray));
            roi_mask(top_y:bottom_y, left_x:right_x) = true;
            
            % Create visualization with ROI overlay
            figure('visible', 'off', 'Position', [100, 100, width, height]);
            imshow(nodapi_img);
            hold on;
            
            % Add semi-transparent ROI overlay
            h = imshow(cat(3, roi_mask*0.8, zeros(size(roi_mask)), zeros(size(roi_mask))));
            set(h, 'AlphaData', roi_mask * 0.3);
            
            % Draw boundaries
            line([1, width], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
            line([1, width], [bottom_y, bottom_y], 'Color', 'g', 'LineWidth', 2);
            line([left_x, left_x], [1, height], 'Color', 'b', 'LineWidth', 2);
            line([right_x, right_x], [1, height], 'Color', 'b', 'LineWidth', 2);
            
            % Add labels
            text(10, top_y - 20, 'Top', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
            text(10, bottom_y + 20, 'Bottom', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
            text(left_x + 10, 30, 'Left', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
            text(right_x - 50, 30, 'Right', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
            
            title('ROI Mask on NoDAPI Image', 'FontSize', 14);
            
            % Save visualization
            roi_viz_dir = fullfile(output_dir, 'ROI_Visualizations');
            if ~exist(roi_viz_dir, 'dir')
                mkdir(roi_viz_dir);
            end
            roi_viz_path = fullfile(roi_viz_dir, 'roi_on_nodapi.jpg');
            saveas(gcf, roi_viz_path);
            disp(['Saved ROI visualization to: ', roi_viz_path]);
            close(gcf);
        catch ME
            disp(['Warning: Error creating ROI visualization: ', ME.message]);
        end
    else
        disp(['NoDAPI image not found: ', nodapi_path]);
    end
else
    disp(['Registration directory not found: ', registration_dir]);
end

% Create growth_plate1 structure
bt1 = struct();
bt1.top_y = top_y;
bt1.bottom_y = bottom_y;
bt1.left_x = left_x;
bt1.right_x = right_x;

% Put in cell array format compatible with existing code
growth_plate1{1,1,1} = 'growth_plate_data';
growth_plate1{1,2,1} = top_y;
growth_plate1{1,3,1} = bottom_y;
growth_plate1{1,4,1} = left_x; % Add left_x to cell array
growth_plate1{1,5,1} = right_x; % Add right_x to cell array

% Save the growth plate data
save(fullfile(output_dir, 'growth_plate1.mat'), 'growth_plate1', 'bt1');
disp(['Saved growth plate data to: ', fullfile(output_dir, 'growth_plate1.mat')]);

disp('Growth plate detection completed successfully.');
end
