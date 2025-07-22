function c_analysis_stack()
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
    home_dir = currDir;
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
    
    % Get growth plate coordinates
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
    
    % Save a copy of the image with growth plate markers
    figure('visible', 'off');
    imshow(img);
    hold on;
    line([1, width], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
    line([1, width], [bottom_y, bottom_y], 'Color', 'g', 'LineWidth', 2);
    title(['Growth plate for ', img_name]);
    saveas(gcf, [img_dir, 'growth_plate_visualization.png']);
    disp(['Saved growth plate visualization to: ', img_dir, 'growth_plate_visualization.png']);
    close(gcf);
    
    % Extract ROI
    if top_y > 0 && bottom_y <= height && top_y < bottom_y
        roi_img = img(top_y:bottom_y, :, :);
        roi_height = bottom_y - top_y;
        roi_width = width;
        roi_area = roi_height * roi_width;
        
        % Save ROI image
        figure('visible', 'off');
        imshow(roi_img);
        title(['ROI for ', img_name]);
        saveas(gcf, [img_dir, 'roi_visualization.png']);
        disp(['Saved ROI visualization to: ', img_dir, 'roi_visualization.png']);
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
disp('Analysis completed successfully.');
end
