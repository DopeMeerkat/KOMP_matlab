function find_growth_plate_stack()
% Modified version of find_growth_plate_2.m to work with CCC_E01_hF_FL1_stack folder
% Instead of using GUI to select coordinates, uses uppermost and lowermost pixel
% coordinates from 2_Dia.png

% Set paths and parameters
home_dir = fullfile(pwd, filesep);
stack_dir = [home_dir, 'CCC_E01_hF_FL1_stack', filesep];
phr_temp = 'CCC_E01_hF';
exp_type = 'F';

% Create dir_info.txt file for compatibility with other scripts
fid = fopen('dir_info.txt', 'w');
fprintf(fid, '%s\n', home_dir);
fprintf(fid, '%s\n', phr_temp);
fclose(fid);
disp(['Created dir_info.txt in ', pwd]);

% Determine delimiter based on platform
if ispc()
    delimeter = '\';
else
    delimeter = '/';
end

% Setup experiment names
exp_name = {[phr_temp,'_h',exp_type,'_F'];[phr_temp,'_h',exp_type,'_M']};  % gene name, exp_name, (female, male)

% Create output directory for storing results
output_dir = [home_dir, 'CCC_E01_hF_FL1_stack_output', delimeter];
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    disp(['Created output directory: ', output_dir]);
end

% Process the 2_Dia.png image to get uppermost and lowermost pixel coordinates
disp('Processing 2_Dia.png to find coordinates...');
img_path = fullfile(stack_dir, '2_Dia.png');

% Check if the image file exists
if ~exist(img_path, 'file')
    error(['Image file not found: ', img_path]);
end

try
    a = imread(img_path);
    disp('Image loaded successfully.');
catch ME
    error(['Error loading image: ', ME.message]);
end

% Find non-zero pixels (object of interest)
try
    gray_img = rgb2gray(a);
catch
    % If conversion fails, the image might already be grayscale
    gray_img = a;
end

[rows, cols] = find(gray_img > 0);
if isempty(rows)
    error('No object found in the image. Check if the image is valid.');
end

% Get the uppermost and lowermost pixels
top_y = min(rows);
bottom_y = max(rows);

disp(['Found coordinates: top_y = ', num2str(top_y), ', bottom_y = ', num2str(bottom_y)]);

% Setup structure similar to original function
no_bone_type = length(exp_name);
bone_type = 1:no_bone_type;
bone_number = {'1'};  % We're only processing one bone/sample in this case
sample_per_bone = 1;  % We're only processing one section in this case

% Create output directory if it doesn't exist
direct = output_dir;

% Process each bone type (though we're just working with one in this case)
for i = 1:no_bone_type
    bt = bone_type(i);
    phr1 = exp_name{bt,:};
    
    % Initialize growth_plate structure
    growth_plate = cell(1, 3);
    growth_plate{1, 1} = bone_number{1};
    
    % Instead of using ginput, use the coordinates we calculated
    growth_plate{1, 2, 1} = top_y;     % top of growth plate
    growth_plate{1, 3, 1} = bottom_y;  % bottom of growth plate
    
    % Save the growth plate data
    eval(['growth_plate', num2str(bt), ' = growth_plate;']);
    save_path = fullfile(direct, ['growth_plate', num2str(bt), '.mat']);
    eval(['save ''', save_path, ''' growth_plate', num2str(bt)]);
    disp(['Saved growth plate data to: ', save_path]);
    
    % Display the result
    figure;
    imshow(a);
    hold on;
    % Draw lines at the detected positions
    line([1, size(a, 2)], [top_y, top_y], 'Color', 'r', 'LineWidth', 2);
    line([1, size(a, 2)], [bottom_y, bottom_y], 'Color', 'g', 'LineWidth', 2);
    title(['Growth plate detection for ', phr1]);
    saveas(gcf, fullfile(direct, ['growth_plate_detection_', num2str(bt), '.png']));
    disp(['Saved visualization to: ', fullfile(direct, ['growth_plate_detection_', num2str(bt), '.png'])]);
end

disp('Growth plate detection completed.');
disp(['Results saved to ', direct]);
end
