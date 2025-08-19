function analyze_stack()
% Modified version of c_analysis.m to work with the output from find_growth_plate_stack
% This function analyzes the growth plate data for the CCC_E01_hF_FL1_stack folder

% Get current directory
disp('Starting analysis...');
currDir = pwd;
disp(['Current directory is: ', currDir]);

% Get the home directory from dir_info.txt
try
    a = textread('dir_info.txt', '%s', 2);
    home_dir = a{1,:};
    phr_temp = a{2,:};
    disp(['Using home directory: ', home_dir]);
    disp(['Using phrase template: ', phr_temp]);
catch ME
    disp(['Error reading dir_info.txt: ', ME.message]);
    disp('Creating default dir_info.txt');
    home_dir = fileparts(currDir);
    if home_dir(end) ~= '/'
        home_dir = [home_dir, '/'];
    end
    phr_temp = 'CCC_E01_hF';
    
    % Create dir_info.txt
    fid = fopen('dir_info.txt', 'w');
    fprintf(fid, '%s\n', home_dir);
    fprintf(fid, '%s\n', phr_temp);
    fclose(fid);
end

% Determine delimiter based on platform
if ispc()
    delimeter = '\';
else
    delimeter = '/';
end

% Set paths
stack_dir = [home_dir, 'CCC_E01_hF_FL1_stack', delimeter];
output_dir = [home_dir, 'CCC_E01_hF_FL1_stack_output', delimeter];
disp(['Stack directory: ', stack_dir]);
disp(['Output directory: ', output_dir]);

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    disp(['Created output directory: ', output_dir]);
end

% Create analysis directory
analysis_dir = [output_dir, 'Analyzed', delimeter];
if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
    disp(['Created analysis directory: ', analysis_dir]);
end

% Try to load growth plate data
try
    growth_plate1_path = [output_dir, 'growth_plate1.mat'];
    growth_plate2_path = [output_dir, 'growth_plate2.mat'];
    
    if exist(growth_plate1_path, 'file')
        load(growth_plate1_path, 'growth_plate1');
        disp('Loaded growth_plate1.mat');
    else
        disp('growth_plate1.mat not found. Using default values.');
        growth_plate1 = {};
    end
    
    if exist(growth_plate2_path, 'file')
        load(growth_plate2_path, 'growth_plate2');
        disp('Loaded growth_plate2.mat');
    else
        disp('growth_plate2.mat not found. Using default values.');
        growth_plate2 = {};
    end
catch ME
    disp(['Error loading growth plate data: ', ME.message]);
    disp('Using default values.');
    growth_plate1 = {};
    growth_plate2 = {};
end

% Process the images in the stack directory
stack_images = dir([stack_dir, '*.png']);
disp(['Found ', num2str(length(stack_images)), ' PNG images in the stack directory']);

% Create a simple output structure
output_file = [analysis_dir, 'analysis_results.txt'];
fid = fopen(output_file, 'w');
fprintf(fid, 'Growth Plate Analysis Results\n');
fprintf(fid, '===========================\n\n');
fprintf(fid, 'Analysis performed on: %s\n\n', datestr(now));

% Create a montage directory for combined visualizations
montage_dir = [analysis_dir, 'ROI_Visualizations', delimeter];
if ~exist(montage_dir, 'dir')
    mkdir(montage_dir);
    disp(['Created montage directory: ', montage_dir]);
end

% First, try to find the mineral image (CCC_E01_hF_FL1_M_s1c2.png) to determine ROI coordinates
mineral_image_idx = -1;
dia_image_idx = -1;
top_y = -1;
bottom_y = -1;
bone_mask = [];

% Look for mineral image first, then Dia as fallback
for i = 1:length(stack_images)
    if contains(stack_images(i).name, 'CCC_E01_hF_FL1_M_s1c2') || contains(stack_images(i).name, '_M_s1c2')
        mineral_image_idx = i;
        break;
    elseif strcmp(stack_images(i).name, '2_Dia.png')
        dia_image_idx = i;
    end
end

% Use the mineral image if found, otherwise try the Dia image
if mineral_image_idx > 0
    disp(['Found mineral image: ', stack_images(mineral_image_idx).name, '. Using it to determine ROI coordinates for all images.']);
    
    % Load the mineral image
    img_path = [stack_dir, stack_images(mineral_image_idx).name];
    mineral_img = imread(img_path);
    [height, width, channels] = size(mineral_img);
    
    % Convert to grayscale if color image
    if channels > 1
        mineral_gray = rgb2gray(mineral_img);
    else
        mineral_gray = mineral_img;
    end
    
    % Apply thresholding to segment the bone
    level = graythresh(mineral_gray);
    bone_binary = imbinarize(mineral_gray, level * 0.7); % Use a lower threshold to capture more of the bone
    
    % Clean up the binary image
    bone_binary = imclose(bone_binary, strel('disk', 5));
    bone_binary = imfill(bone_binary, 'holes');
    
    % Find the largest connected component (the bone)
    cc = bwconncomp(bone_binary);
    stats = regionprops(cc, 'Area', 'BoundingBox', 'Image', 'PixelIdxList');
    
    if ~isempty(stats)
        [~, idx] = max([stats.Area]);
        
        % Create initial bone mask
        bone_mask = false(size(bone_binary));
        bone_mask(stats(idx).PixelIdxList) = true;
        
        % Use find_periosteum_4 to get a more precise bone boundary that follows the periosteum
        try
            disp('Applying find_periosteum_4 to refine bone boundary...');
            [bone_mask_refined, bone_filled] = find_periosteum_4(bone_mask);
            disp('Successfully refined bone boundary using find_periosteum_4');
            
            % Use the refined mask if it's not empty
            if any(bone_mask_refined(:))
                bone_mask = bone_mask_refined;
                % Also try the filled version for visualization purposes
                bone_filled_viz = bone_filled;
            else
                disp('Warning: Refined bone mask is empty. Using original mask.');
                bone_filled_viz = bone_mask;
            end
        catch ME
            disp(['Warning: Error in find_periosteum_4: ', ME.message]);
            disp('Using original bone mask.');
            bone_filled_viz = bone_mask;
        end
        
        % Get bounding box of the refined mask
        stats_refined = regionprops(bone_mask, 'BoundingBox');
        if ~isempty(stats_refined)
            bbox = stats_refined.BoundingBox;
            left_x = max(1, round(bbox(1)));
            top_bbox = max(1, round(bbox(2)));
            bone_width = round(bbox(3));
            bone_height = round(bbox(4));
            right_x = min(width, left_x + bone_width - 1);
            bottom_bbox = min(height, top_bbox + bone_height - 1);
            
            disp(['Refined bone bounding box: left=', num2str(left_x), ', right=', num2str(right_x), ...
                  ', top=', num2str(top_bbox), ', bottom=', num2str(bottom_bbox)]);
        else
            disp('Could not get bounding box of refined mask. Using original.');
            bbox = stats(idx).BoundingBox;
            left_x = max(1, round(bbox(1)));
            top_bbox = max(1, round(bbox(2)));
            bone_width = round(bbox(3));
            bone_height = round(bbox(4));
            right_x = min(width, left_x + bone_width - 1);
            bottom_bbox = min(height, top_bbox + bone_height - 1);
        end
        
        % Try to identify cortical bone regions
        try
            disp('Applying find_cortical to identify cortical bone regions...');
            temp_roi = bone_mask;
            cortical_mask = find_cortical(bone_binary, bone_mask, temp_roi);
            disp('Successfully identified cortical regions.');
            
            % Use cortical mask for visualization
            if any(cortical_mask(:))
                cortical_viz = cortical_mask;
            else
                disp('Warning: Cortical mask is empty.');
                cortical_viz = [];
            end
        catch ME
            disp(['Warning: Error in find_cortical: ', ME.message]);
            cortical_viz = [];
        end
    else
        disp('Could not segment bone. Using full image.');
        bone_mask = true(size(bone_binary));
        bone_filled_viz = bone_mask;
        left_x = 1;
        right_x = width;
        cortical_viz = [];
    end
    
    % Get growth plate coordinates from growth_plate1 (if available)
    if ~isempty(growth_plate1) && isfield(growth_plate1, 'bt1')
        top_y = growth_plate1.bt1.top_y;
        bottom_y = growth_plate1.bt1.bottom_y;
        
        % Check for left_x and right_x in bt1
        if isfield(growth_plate1.bt1, 'left_x') && isfield(growth_plate1.bt1, 'right_x')
            left_x = growth_plate1.bt1.left_x;
            right_x = growth_plate1.bt1.right_x;
            disp(['Using all growth plate coordinates from structure: top_y=', num2str(top_y), ...
                  ', bottom_y=', num2str(bottom_y), ', left_x=', num2str(left_x), ...
                  ', right_x=', num2str(right_x)]);
        else
            disp(['Using vertical growth plate coordinates from structure: top_y=', num2str(top_y), ...
                  ', bottom_y=', num2str(bottom_y)]);
            disp('Horizontal boundaries not found in structure.');
        end
    else
        % Use the coordinates from find_growth_plate_stack if available
        try
            disp('Attempting to access growth_plate1 data...');
            disp(['growth_plate1 type: ', class(growth_plate1)]);
            
            if iscell(growth_plate1) && numel(growth_plate1) >= 1
                disp(['growth_plate1 size: ', num2str(size(growth_plate1))]);
                top_y = growth_plate1{1, 2, 1};
                bottom_y = growth_plate1{1, 3, 1};
                
                % Check for left_x and right_x in cell array
                if size(growth_plate1, 2) >= 5
                    left_x = growth_plate1{1, 4, 1};
                    right_x = growth_plate1{1, 5, 1};
                    disp(['Found all growth plate coordinates: top_y=', num2str(top_y), ...
                          ', bottom_y=', num2str(bottom_y), ', left_x=', num2str(left_x), ...
                          ', right_x=', num2str(right_x)]);
                else
                    disp(['Found vertical growth plate coordinates: top_y=', num2str(top_y), ...
                          ', bottom_y=', num2str(bottom_y)]);
                    disp('Horizontal boundaries not found in cell array.');
                end
            else
                disp('growth_plate1 is not a cell array or is empty');
                % Use default values based on image height
                top_y = round(height * 0.25);
                bottom_y = round(height * 0.75);
                disp(['Using default growth plate coordinates: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
            end
        catch ME
            disp(['Error accessing growth plate data: ', ME.message]);
            % Use default values based on image height
            top_y = round(height * 0.25);
            bottom_y = round(height * 0.75);
            disp(['Using default growth plate coordinates: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
        end
    end
    
    % Create ROI mask using both vertical and horizontal boundaries
    roi_mask = bone_mask;
    
    % Apply vertical boundaries
    roi_mask(1:top_y-1, :) = false;  % Clear above top boundary
    roi_mask(bottom_y+1:end, :) = false;  % Clear below bottom boundary
    
    % Apply horizontal boundaries if available
    if exist('left_x', 'var') && exist('right_x', 'var')
        roi_mask(:, 1:left_x-1) = false;  % Clear left of left boundary
        roi_mask(:, right_x+1:end) = false;  % Clear right of right boundary
        disp('Applied both vertical and horizontal boundaries to ROI mask.');
    else
        disp('Applied only vertical boundaries to ROI mask. Using bone mask for horizontal boundaries.');
    end
    
    % Create a visualization of the bone segmentation
    figure('visible', 'off', 'Position', [100, 100, width, height]);
    
    % Create a subfigure for the original image with bone boundary
    subplot(2, 2, 1);
    imshow(mineral_img);
    hold on;
    
    % Get the boundary of the bone mask
    [B, ~] = bwboundaries(bone_mask, 'noholes');
    
    % Plot each boundary
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
    end
    
    % Draw horizontal lines for the growth plate
    line([1, width], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
    line([1, width], [bottom_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
    
    % Draw vertical lines for left/right boundaries if available
    if exist('left_x', 'var') && exist('right_x', 'var')
        line([left_x, left_x], [1, height], 'Color', 'g', 'LineWidth', 2);
        line([right_x, right_x], [1, height], 'Color', 'g', 'LineWidth', 2);
    end
    
    title('Bone Boundary', 'FontSize', 12);
    
    % Create a subfigure for the binary bone mask
    subplot(2, 2, 2);
    imshow(bone_mask);
    title('Bone Mask', 'FontSize', 12);
    
    % Create a subfigure for the ROI mask
    subplot(2, 2, 3);
    imshow(roi_mask);
    hold on;
    line([1, width], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
    line([1, width], [bottom_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
    
    % Draw vertical lines for left/right boundaries if available
    if exist('left_x', 'var') && exist('right_x', 'var')
        line([left_x, left_x], [1, height], 'Color', 'g', 'LineWidth', 2);
        line([right_x, right_x], [1, height], 'Color', 'g', 'LineWidth', 2);
        title('ROI Mask (All Boundaries)', 'FontSize', 12);
    else
        title('ROI Mask', 'FontSize', 12);
    end
    
    % Create a subfigure for the mineral image with ROI overlay
    subplot(2, 2, 4);
    imshow(mineral_img);
    hold on;
    
    % Create a semi-transparent overlay of the ROI mask
    h = imshow(cat(3, zeros(size(roi_mask)), ones(size(roi_mask)), zeros(size(roi_mask))));
    set(h, 'AlphaData', roi_mask * 0.3);
    
    % Get the boundary of the ROI mask
    [B_roi, ~] = bwboundaries(roi_mask, 'noholes');
    
    % Plot each ROI boundary
    for k = 1:length(B_roi)
        boundary = B_roi{k};
        plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 2);
    end
    
    % Draw horizontal lines for the growth plate
    line([1, width], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
    line([1, width], [bottom_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
    
    % Draw vertical lines for left/right boundaries if available
    if exist('left_x', 'var') && exist('right_x', 'var')
        line([left_x, left_x], [1, height], 'Color', 'g', 'LineWidth', 2);
        line([right_x, right_x], [1, height], 'Color', 'g', 'LineWidth', 2);
    end
    
    title('ROI Overlay', 'FontSize', 12);
    
    % Add legend to the right of the subplot
    annotation('textbox', [0.75, 0.02, 0.2, 0.1], 'String', 'Bone Boundary (Green)', 'FitBoxToText', 'on', 'EdgeColor', 'none');
    annotation('textbox', [0.75, 0.06, 0.2, 0.1], 'String', 'ROI Boundary (Yellow)', 'FitBoxToText', 'on', 'EdgeColor', 'none');
    annotation('textbox', [0.75, 0.10, 0.2, 0.1], 'String', 'Top Growth Plate (Red)', 'FitBoxToText', 'on', 'EdgeColor', 'none');
    annotation('textbox', [0.75, 0.14, 0.2, 0.1], 'String', 'Bottom Growth Plate (Blue)', 'FitBoxToText', 'on', 'EdgeColor', 'none');
    
    if exist('left_x', 'var') && exist('right_x', 'var')
        annotation('textbox', [0.75, 0.18, 0.2, 0.1], 'String', 'Left/Right Boundaries (Green)', 'FitBoxToText', 'on', 'EdgeColor', 'none');
    end
    
    % Add title to the entire figure
    sgtitle('Bone Segmentation and ROI Boundaries', 'FontSize', 14);
    
    % Save the visualization
    saveas(gcf, [montage_dir, 'bone_segmentation.png']);
    disp(['Saved bone segmentation visualization to: ', montage_dir, 'bone_segmentation.png']);
    close(gcf);
    
    % Create HTML report for better visualization
    html_content = ['<!DOCTYPE html>\n<html>\n<head>\n<title>ROI Analysis Report</title>\n' ...
                   '<style>\n' ...
                   'body { font-family: Arial, sans-serif; margin: 20px; }\n' ...
                   'h1 { color: #333366; }\n' ...
                   'img { max-width: 100%; border: 1px solid #ddd; margin: 10px 0; }\n' ...
                   '.info { background-color: #f0f0f0; padding: 10px; border-radius: 5px; }\n' ...
                   '</style>\n</head>\n<body>\n' ...
                   '<h1>ROI Analysis Report</h1>\n' ...
                   '<div class="info">\n' ...
                   '<p><b>Stack:</b> ', stack_dir, '</p>\n' ...
                   '<p><b>Mineral Image:</b> ', stack_images(mineral_image_idx).name, '</p>\n' ...
                   '<p><b>ROI Coordinates:</b></p>\n' ...
                   '<ul>\n' ...
                   '<li>Left: ', num2str(left_x), '</li>\n' ...
                   '<li>Right: ', num2str(right_x), '</li>\n' ...
                   '<li>Top: ', num2str(top_y), '</li>\n' ...
                   '<li>Bottom: ', num2str(bottom_y), '</li>\n' ...
                   '</ul>\n' ...
                   '</div>\n' ...
                   '<h2>Bone Segmentation and ROI</h2>\n' ...
                   '<img src="bone_segmentation.png" alt="Bone Segmentation">\n' ...
                   '<p>This visualization shows:</p>\n' ...
                   '<ul>\n' ...
                   '<li><b>Top Left:</b> Original mineral image with bone boundary (green) and growth plate lines (red/blue)</li>\n' ...
                   '<li><b>Top Right:</b> Binary bone mask</li>\n' ...
                   '<li><b>Bottom Left:</b> ROI mask (bone mask + growth plate boundaries)</li>\n' ...
                   '<li><b>Bottom Right:</b> Original image with semi-transparent ROI overlay</li>\n' ...
                   '</ul>\n' ...
                   '</body>\n</html>'];
    
    % Write HTML file
    html_file = [montage_dir, 'roi_report.html'];
    fid = fopen(html_file, 'w');
    fprintf(fid, '%s', html_content);
    fclose(fid);
    disp(['Created HTML report: ', html_file]);
    
    % Open HTML report in browser if on macOS
    if ismac
        try
            system(['open ', html_file]);
        catch
            disp('Could not automatically open HTML report. Please open it manually.');
        end
    end
elseif dia_image_idx > 0
    disp('Found 2_Dia.png. Using it to determine ROI coordinates for all images.');
    
    % Load the Dia image
    img_path = [stack_dir, stack_images(dia_image_idx).name];
    reference_img = imread(img_path);
    [height, width, ~] = size(reference_img);
    
    % Use default full width
    left_x = 1;
    right_x = width;
    
    % Try to detect horizontal boundaries from Dia image too
    try
        % Convert to grayscale if it's a color image
        if size(reference_img, 3) > 1
            gray_img = rgb2gray(reference_img);
        else
            gray_img = reference_img;
        end
        
        % Basic edge detection to find bone boundaries
        col_sums = sum(gray_img, 1);
        col_sums = col_sums / max(col_sums);
        col_sums_smooth = movmean(col_sums, width/20);
        
        % Find edges
        threshold = 0.1;
        left_indices = find(col_sums_smooth > threshold);
        right_indices = find(col_sums_smooth > threshold);
        
        if ~isempty(left_indices) && ~isempty(right_indices)
            left_x = left_indices(1);
            right_x = right_indices(end);
            
            % Add margins
            margin = round(width * 0.02);
            left_x = max(1, left_x - margin);
            right_x = min(width, right_x + margin);
            
            disp(['Found horizontal boundaries from Dia image: left_x=', num2str(left_x), ', right_x=', num2str(right_x)]);
        end
    catch ME
        disp(['Could not detect horizontal boundaries from Dia image: ', ME.message]);
    end
else
    disp('Warning: Neither mineral image nor 2_Dia.png found. Will determine ROI coordinates for each image separately.');
    left_x = -1;
    right_x = -1;
end

% Process each image
for i = 1:length(stack_images)
    img_name = stack_images(i).name;
    disp(['Processing image: ', img_name]);
    
    % Create a folder for this image
    img_dir = [analysis_dir, img_name(1:end-4), delimeter];
    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    
    % Load the image
    img_path = [stack_dir, img_name];
    img = imread(img_path);
    disp(['Loaded image: ', img_path]);
    
    % Get image dimensions
    [height, width, channels] = size(img);
    disp(['Image dimensions: ', num2str(height), 'x', num2str(width), 'x', num2str(channels)]);
    
    % Use pre-determined coordinates from 2_Dia.png if available, otherwise determine them for this image
    if top_y > 0 && bottom_y > 0
        disp(['Using pre-determined growth plate coordinates: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
    else
        % Determine growth plate coordinates for this specific image
        if ~isempty(growth_plate1) && isfield(growth_plate1, 'bt1')
            top_y = growth_plate1.bt1.top_y;
            bottom_y = growth_plate1.bt1.bottom_y;
            disp(['Using growth plate coordinates from file: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
        else
            % Use the coordinates from find_growth_plate_stack if available
            try
                disp('Attempting to access growth_plate1 data...');
                disp(['growth_plate1 type: ', class(growth_plate1)]);
                
                if iscell(growth_plate1) && numel(growth_plate1) >= 1
                    disp(['growth_plate1 size: ', num2str(size(growth_plate1))]);
                    top_y = growth_plate1{1, 2, 1};
                    bottom_y = growth_plate1{1, 3, 1};
                    disp(['Found growth plate coordinates: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
                else
                    disp('growth_plate1 is not a cell array or is empty');
                    % Use default values based on image height
                    top_y = round(height * 0.25);
                    bottom_y = round(height * 0.75);
                    disp(['Using default growth plate coordinates: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
                end
            catch ME
                disp(['Error accessing growth plate data: ', ME.message]);
                % Use default values based on image height
                top_y = round(height * 0.25);
                bottom_y = round(height * 0.75);
                disp(['Using default growth plate coordinates: top_y=', num2str(top_y), ', bottom_y=', num2str(bottom_y)]);
            end
        end
    end
    
    % Determine horizontal boundaries for this image if not already set
    if left_x <= 0 || right_x <= 0
        left_x = 1;
        right_x = width;
        
        % Try to detect from current image if it's a mineral or Dia image
        if contains(img_name, '_M_s1c2') || strcmp(img_name, '2_Dia.png')
            try
                % Convert to grayscale if it's a color image
                if size(img, 3) > 1
                    gray_img = rgb2gray(img);
                else
                    gray_img = img;
                end
                
                % Find edges
                col_sums = sum(gray_img, 1);
                col_sums = col_sums / max(col_sums);
                col_sums_smooth = movmean(col_sums, width/20);
                
                % Find boundaries
                threshold = 0.1;
                left_indices = find(col_sums_smooth > threshold);
                right_indices = find(col_sums_smooth > threshold);
                
                if ~isempty(left_indices) && ~isempty(right_indices)
                    left_x = left_indices(1);
                    right_x = right_indices(end);
                    
                    % Add margins
                    margin = round(width * 0.02);
                    left_x = max(1, left_x - margin);
                    right_x = min(width, right_x + margin);
                    
                    disp(['Detected horizontal boundaries for ', img_name, ': left_x=', num2str(left_x), ', right_x=', num2str(right_x)]);
                end
            catch
                disp(['Could not detect horizontal boundaries for ', img_name]);
            end
        end
    end
    
    % Save a copy of the image with growth plate markers
    figure('visible', 'off', 'Position', [100, 100, width, height]);
    imshow(img);
    hold on;
    
    % Draw top boundary with text label
    line([left_x, right_x], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
    text(left_x + 20, top_y - 20, 'Top Boundary', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
    
    % Draw bottom boundary with text label
    line([left_x, right_x], [bottom_y, bottom_y], 'Color', 'g', 'LineWidth', 2);
    text(left_x + 20, bottom_y + 20, 'Bottom Boundary', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
    
    % Draw vertical lines at the edges of the ROI
    line([left_x, left_x], [top_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
    text(left_x + 5, (top_y + bottom_y)/2, 'Left', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7], 'Rotation', 90);
    
    line([right_x, right_x], [top_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
    text(right_x - 20, (top_y + bottom_y)/2, 'Right', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7], 'Rotation', 90);
    
    % Add a title with image name and ROI info
    roi_width = right_x - left_x;
    roi_height = bottom_y - top_y;
    title_str = sprintf('%s (ROI: %dx%d pixels)', img_name, roi_width, roi_height);
    title(title_str, 'FontSize', 14);
    
    % Save the visualization to the image folder
    saveas(gcf, [img_dir, 'growth_plate_visualization.png']);
    disp(['Saved growth plate visualization to: ', img_dir, 'growth_plate_visualization.png']);
    
    % Also save to the montage directory
    saveas(gcf, [montage_dir, img_name(1:end-4), '_boundaries.png']);
    
    % Close the current figure
    close(gcf);
    
    % Create a highlighted version to better show the ROI
    figure('visible', 'off', 'Position', [100, 100, width, height]);
    imshow(img);
    hold on;
    
    % Create a semi-transparent overlay for the ROI area
    roi_overlay = zeros(height, width, 'uint8');
    for y = top_y:bottom_y
        for x = left_x:right_x
            roi_overlay(y, x) = 1;
        end
    end
    
    % Apply color overlay to highlight ROI
    colored_overlay = cat(3, zeros(size(roi_overlay)), roi_overlay, zeros(size(roi_overlay)));
    h = imshow(colored_overlay);
    set(h, 'AlphaData', 0.3 * roi_overlay);
    
    % Draw boundary lines
    line([left_x, right_x], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
    line([left_x, right_x], [bottom_y, bottom_y], 'Color', 'g', 'LineWidth', 2);
    line([left_x, left_x], [top_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
    line([right_x, right_x], [top_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
    
    % Add a title
    title(['ROI Highlight for ', img_name], 'FontSize', 14);
    
    % Save the highlighted visualization
    saveas(gcf, [montage_dir, img_name(1:end-4), '_roi_highlight.png']);
    disp(['Saved highlighted ROI visualization to: ', montage_dir, img_name(1:end-4), '_roi_highlight.png']);
    
    % Close the highlighted figure
    close(gcf);
    
    % Extract ROI using both horizontal and vertical boundaries
    if top_y > 0 && bottom_y <= height && top_y < bottom_y && left_x > 0 && right_x <= width && left_x < right_x
        roi_img = img(top_y:bottom_y, left_x:right_x, :);
        roi_height = bottom_y - top_y;
        roi_width = right_x - left_x;
        roi_area = roi_height * roi_width;
        
        % Save ROI image
        figure('visible', 'off', 'Position', [100, 100, roi_width, roi_height]);
        imshow(roi_img);
        title(['ROI for ', img_name], 'FontSize', 14);
        saveas(gcf, [img_dir, 'roi_visualization.png']);
        disp(['Saved ROI visualization to: ', img_dir, 'roi_visualization.png']);
        
        % Also save to the montage directory
        saveas(gcf, [montage_dir, img_name(1:end-4), '_roi_only.png']);
        close(gcf);
        
        % Create a side-by-side comparison (original vs ROI)
        figure('visible', 'off', 'Position', [100, 100, width*2, height]);
        
        % First subplot - original image with ROI boundaries
        subplot(1, 2, 1);
        imshow(img);
        hold on;
        rectangle('Position', [left_x, top_y, right_x-left_x, bottom_y-top_y], ...
                 'EdgeColor', 'y', 'LineWidth', 2, 'LineStyle', '--');
        line([left_x, right_x], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
        line([left_x, right_x], [bottom_y, bottom_y], 'Color', 'g', 'LineWidth', 2);
        line([left_x, left_x], [top_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
        line([right_x, right_x], [top_y, bottom_y], 'Color', 'b', 'LineWidth', 2);
        title('Original Image with ROI Boundaries', 'FontSize', 12);
        
        % Second subplot - extracted ROI
        subplot(1, 2, 2);
        imshow(roi_img);
        title('Extracted ROI', 'FontSize', 12);
        
        % Add a common title
        sgtitle([img_name, ' - Region of Interest Analysis'], 'FontSize', 14);
        
        % Save the side-by-side comparison
        saveas(gcf, [montage_dir, img_name(1:end-4), '_comparison.png']);
        disp(['Saved comparison visualization to: ', montage_dir, img_name(1:end-4), '_comparison.png']);
        close(gcf);
        
        % Calculate basic statistics for each channel
        mean_values = zeros(1, channels);
        std_values = zeros(1, channels);
        
        for c = 1:channels
            channel_data = double(roi_img(:,:,c));
            mean_values(c) = mean(channel_data(:));
            std_values(c) = std(double(channel_data(:)));
        end
        
        % Save results to text file
        result_file = [img_dir, 'analysis_results.txt'];
        fid2 = fopen(result_file, 'w');
        fprintf(fid2, 'Analysis Results for %s\n', img_name);
        fprintf(fid2, '===========================\n\n');
        fprintf(fid2, 'Image Dimensions: %d x %d x %d\n', height, width, channels);
        fprintf(fid2, 'Growth Plate Coordinates:\n');
        fprintf(fid2, '  Top Y: %d\n', top_y);
        fprintf(fid2, '  Bottom Y: %d\n', bottom_y);
        fprintf(fid2, 'ROI Dimensions: %d x %d\n', roi_height, roi_width);
        fprintf(fid2, 'ROI Area: %d pixels\n\n', roi_area);
        fprintf(fid2, 'Channel Statistics:\n');
        
        for c = 1:channels
            fprintf(fid2, '  Channel %d:\n', c);
            fprintf(fid2, '    Mean Value: %.2f\n', mean_values(c));
            fprintf(fid2, '    Standard Deviation: %.2f\n', std_values(c));
        end
        
        fclose(fid2);
        disp(['Saved analysis results to: ', result_file]);
        
        % Add summary to main output file
        fprintf(fid, 'Image: %s\n', img_name);
        fprintf(fid, '  ROI Height: %d pixels\n', roi_height);
        fprintf(fid, '  ROI Area: %d pixels\n', roi_area);
        for c = 1:channels
            fprintf(fid, '  Mean Value Channel %d: %.2f\n', c, mean_values(c));
        end
        fprintf(fid, '\n');
    else
        disp('Invalid growth plate coordinates. Skipping ROI extraction.');
        fprintf(fid, 'Image: %s - Invalid growth plate coordinates\n\n', img_name);
    end
end

fclose(fid);
disp(['Saved overall analysis results to: ', output_file]);

% Create a summary visualization of all ROIs
disp('Creating summary visualization of all ROIs...');
try
    % Get all comparison images
    comparison_images = dir([montage_dir, '*_comparison.png']);
    
    if length(comparison_images) > 0
        % Create a figure for the summary
        figure('visible', 'off', 'Position', [100, 100, 1200, 800]);
        
        % Determine grid size
        num_images = length(comparison_images);
        grid_cols = min(3, num_images);
        grid_rows = ceil(num_images / grid_cols);
        
        % Plot each image in the grid
        for i = 1:num_images
            img_path = [montage_dir, comparison_images(i).name];
            img = imread(img_path);
            
            subplot(grid_rows, grid_cols, i);
            imshow(img);
            title(comparison_images(i).name(1:end-14), 'FontSize', 10);
        end
        
        % Add overall title
        sgtitle('Summary of ROI Analysis for All Images', 'FontSize', 16);
        
        % Save the summary visualization
        summary_path = [analysis_dir, 'roi_analysis_summary.png'];
        saveas(gcf, summary_path);
        disp(['Saved ROI analysis summary to: ', summary_path]);
        close(gcf);
        
        % Create an HTML report for easy viewing
        html_file = [analysis_dir, 'roi_analysis_report.html'];
        fid_html = fopen(html_file, 'w');
        
        fprintf(fid_html, '<!DOCTYPE html>\n');
        fprintf(fid_html, '<html>\n');
        fprintf(fid_html, '<head>\n');
        fprintf(fid_html, '  <title>ROI Analysis Report</title>\n');
        fprintf(fid_html, '  <style>\n');
        fprintf(fid_html, '    body { font-family: Arial, sans-serif; margin: 20px; }\n');
        fprintf(fid_html, '    h1 { color: #333366; }\n');
        fprintf(fid_html, '    h2 { color: #666699; margin-top: 30px; }\n');
        fprintf(fid_html, '    .image-container { margin: 20px 0; }\n');
        fprintf(fid_html, '    img { max-width: 100%%; border: 1px solid #ddd; }\n');
        fprintf(fid_html, '  </style>\n');
        fprintf(fid_html, '</head>\n');
        fprintf(fid_html, '<body>\n');
        fprintf(fid_html, '  <h1>ROI Analysis Report</h1>\n');
        fprintf(fid_html, '  <p>Analysis performed on: %s</p>\n', datestr(now));
        
        fprintf(fid_html, '  <h2>Summary Visualization</h2>\n');
        fprintf(fid_html, '  <div class="image-container">\n');
        fprintf(fid_html, '    <img src="roi_analysis_summary.png" alt="ROI Analysis Summary">\n');
        fprintf(fid_html, '  </div>\n');
        
        fprintf(fid_html, '  <h2>Individual Image Analysis</h2>\n');
        
        % Add all individual images
        all_images = dir([montage_dir, '*.png']);
        image_names = {all_images.name};
        
        % Group images by their base name
        base_names = unique(cellfun(@(x) x(1:end-18), {comparison_images.name}, 'UniformOutput', false));
        
        for i = 1:length(base_names)
            base_name = base_names{i};
            fprintf(fid_html, '  <h3>%s</h3>\n', base_name);
            
            % Find all visualizations for this image
            boundary_img = [base_name, '_boundaries.png'];
            highlight_img = [base_name, '_roi_highlight.png'];
            roi_img = [base_name, '_roi_only.png'];
            comparison_img = [base_name, '_comparison.png'];
            
            % Add each image if it exists
            if any(strcmp(image_names, boundary_img))
                fprintf(fid_html, '  <div class="image-container">\n');
                fprintf(fid_html, '    <h4>Boundaries</h4>\n');
                fprintf(fid_html, '    <img src="ROI_Visualizations/%s" alt="%s">\n', boundary_img, boundary_img);
                fprintf(fid_html, '  </div>\n');
            end
            
            if any(strcmp(image_names, highlight_img))
                fprintf(fid_html, '  <div class="image-container">\n');
                fprintf(fid_html, '    <h4>ROI Highlight</h4>\n');
                fprintf(fid_html, '    <img src="ROI_Visualizations/%s" alt="%s">\n', highlight_img, highlight_img);
                fprintf(fid_html, '  </div>\n');
            end
            
            if any(strcmp(image_names, comparison_img))
                fprintf(fid_html, '  <div class="image-container">\n');
                fprintf(fid_html, '    <h4>Side-by-Side Comparison</h4>\n');
                fprintf(fid_html, '    <img src="ROI_Visualizations/%s" alt="%s">\n', comparison_img, comparison_img);
                fprintf(fid_html, '  </div>\n');
            end
            
            if any(strcmp(image_names, roi_img))
                fprintf(fid_html, '  <div class="image-container">\n');
                fprintf(fid_html, '    <h4>Extracted ROI</h4>\n');
                fprintf(fid_html, '    <img src="ROI_Visualizations/%s" alt="%s">\n', roi_img, roi_img);
                fprintf(fid_html, '  </div>\n');
            end
        end
        
        fprintf(fid_html, '</body>\n');
        fprintf(fid_html, '</html>\n');
        
        fclose(fid_html);
        disp(['Created HTML report: ', html_file]);
    else
        disp('No comparison images found for summary visualization.');
    end
catch ME
    disp(['Error creating summary visualization: ', ME.message]);
end

disp('Analysis completed successfully.');
end
