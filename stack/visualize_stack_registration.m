function visualize_stack_registration()
% This function visualizes the registration of CCC_E01_hF_FL1_stack
% It assumes the images are already aligned and generates various visualizations:
%   1. Combined registration visualization of all images
%   2. Individual registration visualizations with grid overlay
%   3. RGB composite of the first 3 images
%   4. Multi-color overlay of all stack images
%   5. Edge-based overlay showing alignment of features
%
% All visualizations are saved in the CCC_E01_hF_FL1_stack_output/Registration directory
%
% Usage:
%   visualize_stack_registration() - Run from the KOMP_matlab directory
%
% No input parameters required

% Determine delimiter based on platform
if ispc()
    delimiter = '\';
else
    delimiter = '/';
end

% Get the home directory
home_dir = pwd;
% Check if we're in the src directory and need to move up one level
if contains(home_dir, [delimiter, 'src']) || endsWith(home_dir, 'src')
    [home_dir, ~, ~] = fileparts(home_dir);
    disp('Detected running from src directory, moving up one level');
end

if home_dir(end) ~= delimiter
    home_dir = [home_dir, delimiter];
end

% Create dir_info.txt file for compatibility with other scripts
fid = fopen('dir_info.txt', 'w');
fprintf(fid, '%s\n', home_dir);
fprintf(fid, '%s\n', 'CCC_E01_hF');
fclose(fid);
disp(['Created dir_info.txt in ', pwd]);

% Set paths
stack_dir = [home_dir, 'CCC_E01_hF_FL1_stack', delimiter];
output_dir = [home_dir, 'CCC_E01_hF_FL1_stack_output', delimiter];
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    disp(['Created output directory: ', output_dir]);
end

% Create registration visualization directory
reg_vis_dir = [output_dir, 'Registration', delimiter];
if ~exist(reg_vis_dir, 'dir')
    mkdir(reg_vis_dir);
    disp(['Created registration visualization directory: ', reg_vis_dir]);
end

% Get a list of all PNG images in the stack directory
stack_images = dir([stack_dir, '*.png']);
disp(['Found ', num2str(length(stack_images)), ' PNG images in the stack directory']);

if isempty(stack_images)
    error('No images found in the stack directory');
end

% Create a combined image showing all the aligned images
figure('visible', 'off', 'Position', [100, 100, 1200, 800]);

% Process each image
num_images = length(stack_images);
num_rows = ceil(sqrt(num_images));
num_cols = ceil(num_images / num_rows);

for i = 1:length(stack_images)
    img_name = stack_images(i).name;
    disp(['Processing image: ', img_name]);
    
    % Load the image
    img_path = [stack_dir, img_name];
    img = imread(img_path);
    
    % Create subplot
    subplot(num_rows, num_cols, i);
    imshow(img);
    title(img_name, 'Interpreter', 'none');
    
    % Save individual visualization
    figure('visible', 'off');
    imshow(img);
    title(['Registration visualization: ', img_name], 'Interpreter', 'none');
    
    % Add a grid to show alignment
    hold on;
    
    % Add horizontal grid lines
    for y = 500:500:size(img, 1)
        line([1, size(img, 2)], [y, y], 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
    end
    
    % Add vertical grid lines
    for x = 500:500:size(img, 2)
        line([x, x], [1, size(img, 1)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
    end
    
    % Save individual image with grid
    saveas(gcf, [reg_vis_dir, img_name(1:end-4), '_registration.png']);
    close(gcf);
end

% Save the combined visualization
sgtitle('Stack Registration Visualization', 'FontSize', 16);
saveas(gcf, [reg_vis_dir, 'combined_registration.png']);
close(gcf);

% Create an RGB composite visualization using the first 3 images (if available)
if length(stack_images) >= 3
    % Create RGB composite
    figure('visible', 'off', 'Position', [100, 100, 1200, 800]);
    
    % Load the first 3 images
    img1 = imread([stack_dir, stack_images(1).name]);
    img2 = imread([stack_dir, stack_images(2).name]);
    img3 = imread([stack_dir, stack_images(3).name]);
    
    % Convert to grayscale if needed
    if size(img1, 3) > 1
        img1 = rgb2gray(img1);
    end
    if size(img2, 3) > 1
        img2 = rgb2gray(img2);
    end
    if size(img3, 3) > 1
        img3 = rgb2gray(img3);
    end
    
    % Normalize each channel
    img1 = double(img1) / double(max(img1(:)));
    img2 = double(img2) / double(max(img2(:)));
    img3 = double(img3) / double(max(img3(:)));
    
    % Create RGB composite
    rgb_img = zeros(size(img1, 1), size(img1, 2), 3);
    rgb_img(:,:,1) = img1;  % Red channel
    rgb_img(:,:,2) = img2;  % Green channel
    rgb_img(:,:,3) = img3;  % Blue channel
    
    % Display and save
    imshow(rgb_img);
    title(['RGB Composite: ', stack_images(1).name, ' (R), ', stack_images(2).name, ' (G), ', stack_images(3).name, ' (B)'], 'Interpreter', 'none');
    saveas(gcf, [reg_vis_dir, 'rgb_composite.png']);
    close(gcf);
end

% Create a multi-color overlay of all images
figure('visible', 'off', 'Position', [100, 100, 1200, 800]);

% First get the size of the first image to initialize the overlay
first_img = imread([stack_dir, stack_images(1).name]);
if size(first_img, 3) > 1
    first_img = rgb2gray(first_img);
end

% Initialize overlay image
overlay_img = zeros(size(first_img));
overlay_color = zeros(size(first_img, 1), size(first_img, 2), 3);

% Define a set of colors to use for the overlay
colors = [
    1.0, 0.0, 0.0;  % Red
    0.0, 1.0, 0.0;  % Green
    0.0, 0.0, 1.0;  % Blue
    1.0, 1.0, 0.0;  % Yellow
    1.0, 0.0, 1.0;  % Magenta
    0.0, 1.0, 1.0;  % Cyan
    0.5, 0.5, 0.5;  % Gray
    1.0, 0.5, 0.0;  % Orange
    0.5, 0.0, 0.5;  % Purple
    0.0, 0.5, 0.5;  % Teal
];

% Process each image
for i = 1:min(length(stack_images), size(colors, 1))
    img_name = stack_images(i).name;
    
    % Load the image
    img = imread([stack_dir, img_name]);
    if size(img, 3) > 1
        img = rgb2gray(img);
    end
    
    % Normalize
    img = double(img) / double(max(img(:)));
    
    % Add to overlay with assigned color
    for c = 1:3
        overlay_color(:,:,c) = overlay_color(:,:,c) + img * colors(i, c);
    end
end

% Normalize the overlay image
max_val = max(overlay_color(:));
if max_val > 0
    overlay_color = overlay_color / max_val;
end

% Display and save
imshow(overlay_color);
title('Multi-Color Overlay of All Stack Images', 'FontSize', 16);

% Add a legend for colors
legend_text = cell(min(length(stack_images), size(colors, 1)), 1);
for i = 1:min(length(stack_images), size(colors, 1))
    legend_text{i} = stack_images(i).name;
end

% Add a text-based legend in the top-left corner
text(10, 20, 'Color Legend:', 'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold');
for i = 1:min(length(stack_images), size(colors, 1))
    text(10, 20 + i*20, [num2str(i), ': ', legend_text{i}], 'Color', colors(i,:), 'FontSize', 12);
end

% Save the overlay
saveas(gcf, [reg_vis_dir, 'multicolor_overlay.png']);
close(gcf);

% Edge-based registration visualization
figure('visible', 'off', 'Position', [100, 100, 1200, 800]);

% Create an edge image for each stack image and overlay them
edge_overlay = zeros(size(first_img, 1), size(first_img, 2), 3);

for i = 1:min(length(stack_images), size(colors, 1))
    img_name = stack_images(i).name;
    
    % Load the image
    img = imread([stack_dir, img_name]);
    if size(img, 3) > 1
        img = rgb2gray(img);
    end
    
    % Detect edges
    edges = edge(img, 'canny');
    
    % Add to overlay with assigned color
    for c = 1:3
        edge_overlay(:,:,c) = edge_overlay(:,:,c) + double(edges) * colors(i, c);
    end
end

% Normalize the edge overlay
max_val = max(edge_overlay(:));
if max_val > 0
    edge_overlay = edge_overlay / max_val;
end

% Display and save
imshow(edge_overlay);
title('Edge-Based Registration Visualization', 'FontSize', 16);

% Add the same legend as before
text(10, 20, 'Edge Color Legend:', 'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold');
for i = 1:min(length(stack_images), size(colors, 1))
    text(10, 20 + i*20, [num2str(i), ': ', legend_text{i}], 'Color', colors(i,:), 'FontSize', 12);
end

% Save the edge overlay
saveas(gcf, [reg_vis_dir, 'edge_overlay.png']);
close(gcf);

% Generate shifted images as in the original registration
disp('Generating shifted images...');

% Use the first non-Epi/Dia image as the reference
reference_idx = -1;
for i = 1:length(stack_images)
    if ~contains(stack_images(i).name, 'Epi') && ~contains(stack_images(i).name, 'Dia')
        reference_idx = i;
        break;
    end
end

if reference_idx == -1
    disp('Warning: Could not find a suitable reference image.');
    return;
end

reference_img = imread([stack_dir, stack_images(reference_idx).name]);

% Process each image to create shift versions
for i = 1:length(stack_images)
    % Skip reference image
    if i == reference_idx
        continue;
    end
    
    img_name = stack_images(i).name;
    % Skip Epi/Dia images
    if contains(img_name, 'Epi') || contains(img_name, 'Dia')
        continue;
    end
    
    disp(['  Creating shift image for: ', img_name]);
    
    % Load the image
    img = imread([stack_dir, img_name]);
    
    % Since the images are already aligned, we don't need to calculate shifts
    % Just save a copy with _shift suffix to maintain compatibility
    [~, name, ~] = fileparts(img_name);
    shift_img_name = [name, '_shift.jpg'];
    
    % Save the shift image with original colors preserved
    imwrite(img, [reg_vis_dir, shift_img_name]);
end

% Create a single NoDAPI overlay image
disp('Creating shift2_NoDAPI overlay...');

% Find the DAPI channel image
dapi_idx = -1;
for i = 1:length(stack_images)
    if contains(stack_images(i).name, 'CCC_E01_hF_FL1_A_s1c2')
        dapi_idx = i;
        break;
    end
end

if dapi_idx == -1
    disp('Warning: Could not find DAPI channel (CCC_E01_hF_FL1_A_s1c2). Using all images for overlay.');
end

% Initialize a blank RGB image for the overlay
sample_img = imread([stack_dir, stack_images(1).name]);
overlay_img = zeros(size(sample_img));

% Count of non-DAPI images used in overlay
overlay_count = 0;

% Add all non-DAPI, non-Epi/Dia images to the overlay
for i = 1:length(stack_images)
    % Skip DAPI, Epi, and Dia images
    if i == dapi_idx || contains(stack_images(i).name, 'Epi') || contains(stack_images(i).name, 'Dia')
        continue;
    end
    
    img = imread([stack_dir, stack_images(i).name]);
    
    % Add to overlay
    overlay_img = overlay_img + double(img);
    overlay_count = overlay_count + 1;
end

% Normalize the overlay
if overlay_count > 0
    overlay_img = uint8(overlay_img / overlay_count);
end

% Save the NoDAPI overlay
shift2_nodapi_name = 'CCC_E01_hF_FL1_shift2_NoDAPI.jpg';
imwrite(overlay_img, [reg_vis_dir, shift2_nodapi_name]);
disp(['  Created ', shift2_nodapi_name, ' (overlay of ', num2str(overlay_count), ' non-DAPI images)']);

disp('Registration visualization completed successfully.');
end
