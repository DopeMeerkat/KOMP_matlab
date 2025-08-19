function control_data_write_stack()
% Modified version of control_data_write_2.m to work with the stack directory structure
% This function writes data from analysis results to an Excel file for the CCC_E01_hF_FL1_stack folder

% Get current directory
disp('Starting data writing process...');
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
stack_dir = [home_dir(1:end-4), 'CCC_E01_hF_FL1_stack', delimeter];
output_dir = [home_dir(1:end-4), 'CCC_E01_hF_FL1_stack_output', delimeter];
analyzed_dir = [output_dir, 'Analyzed', delimeter];
disp(['Stack directory: ', stack_dir]);
disp(['Output directory: ', output_dir]);
disp(['Analyzed directory: ', analyzed_dir]);

% Extract bone type from the directory name (F for femur)
bone_type = 'F';
exp_name = 'CCC_E01';
common = [exp_name, '_h', bone_type, '_'];

% Define phrases for female/male
phrase2 = {'Female'; 'Male'};

% Create the phrases for the Excel file
phrase1 = {[common, 'F']; [common, 'M']};

% Write data to Excel
disp('Writing data to Excel...');
write_data6_cell_stack(analyzed_dir, phrase1, phrase2, bone_type, delimeter, home_dir);

% Add a note about ROI visualizations
roi_vis_dir = [analyzed_dir, 'ROI_Visualizations', delimeter];
if exist(roi_vis_dir, 'dir')
    disp('ROI visualization directory found.');
    % HTML report generation has been removed as per requirements
end

disp('Data writing process completed.');
end

function write_data6_cell_stack(direct, ~, ~, ~, delimeter, home_dir)
% Modified version of write_data6_cell for the stack workflow
% Parameters phrase1, phrase2, and bone_type are kept for compatibility but not used
disp('Using write_data6_cell_stack function...');

% Excel template path
template_path = [fileparts(home_dir), '/Excel_template_withGraphs.xlsx'];
disp(['Using Excel template: ', template_path]);

% Define output files
output_excel = [direct, 'analysis_results.xlsx'];
output_csv = [direct, 'analysis_results.csv'];

% Check if running on macOS (which often has issues with xlswrite)
if ismac
    disp('Running on macOS - will use CSV format by default');
    use_excel = false;
else
    % For other platforms, try to use Excel first
    disp(['Creating Excel file: ', output_excel]);
    use_excel = true;
    
    % Copy template to output
    if exist(template_path, 'file') && use_excel
        try
            copyfile(template_path, output_excel);
            disp('Excel template copied successfully.');
        catch ME
            disp(['Error copying Excel template: ', ME.message]);
            disp('Will use CSV format instead.');
            use_excel = false;
        end
    else
        if ~exist(template_path, 'file')
            disp('Excel template not found.');
        end
        
        % Try to create a basic Excel file to test if Excel writing works
        try
            warning('off', 'MATLAB:xlswrite:AddSheet');
            result = xlswrite(output_excel, {'Analysis Results'}, 'Sheet1', 'A1');
            if ~result
                disp('Excel writing not supported on this system. Will use CSV format instead.');
                use_excel = false;
            end
        catch
            disp('Excel writing failed. Will use CSV format instead.');
            use_excel = false;
        end
    end
end

% Read analysis results from the text file
analysis_file = [direct, 'analysis_results.txt'];
disp(['Reading analysis results from: ', analysis_file]);

if ~exist(analysis_file, 'file')
    disp('Analysis results file not found.');
    
    % Check if individual image analysis folders exist
    img_folders = dir(direct);
    img_folders = img_folders([img_folders.isdir]);
    img_folders = img_folders(~ismember({img_folders.name}, {'.', '..'}));
    
    if ~isempty(img_folders)
        disp('Found individual image analysis folders. Attempting to compile results...');
        compiled_data = {};
        line_num = 1;
        
        % Loop through each image folder and read its analysis_results.txt file
        for i = 1:length(img_folders)
            img_folder = img_folders(i).name;
            img_result_file = [direct, img_folder, delimeter, 'analysis_results.txt'];
            
            if exist(img_result_file, 'file')
                % Add image header to compiled data
                compiled_data{line_num} = ['Image: ', img_folder];
                line_num = line_num + 1;
                
                % Read image analysis file
                img_fid = fopen(img_result_file, 'r');
                while ~feof(img_fid)
                    line = fgetl(img_fid);
                    if ~isempty(line)
                        compiled_data{line_num} = line;
                        line_num = line_num + 1;
                    end
                end
                fclose(img_fid);
                
                % Add a blank line between images
                compiled_data{line_num} = '';
                line_num = line_num + 1;
            end
        end
        
        if line_num > 1
            % Write compiled data to analysis_results.txt
            disp(['Creating compiled analysis file: ', analysis_file]);
            compile_fid = fopen(analysis_file, 'w');
            for i = 1:length(compiled_data)
                fprintf(compile_fid, '%s\n', compiled_data{i});
            end
            fclose(compile_fid);
            
            % Now we have an analysis file, so continue with data processing
            data = compiled_data;
        else
            disp('No analysis results found in image folders.');
            return;
        end
    else
        disp('No analysis data found. Run analysis first.');
        return;
    end
else
    % Read the existing analysis file
    fid = fopen(analysis_file, 'r');
    data = {};
    line_num = 1;
    
    % Read all lines from the file
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line)
            data{line_num} = line;
            line_num = line_num + 1;
        end
    end
    fclose(fid);
end

% Write to Excel or CSV based on platform support
try
    % Create headers
    headers = {'Image', 'Parameter', 'Value'};
    
    if use_excel
        % Try writing to Excel
        warning('off', 'MATLAB:xlswrite:AddSheet'); % Suppress sheet warnings
        excel_success = xlswrite(output_excel, headers, 'Analysis', 'A1');
        if ~excel_success
            disp('Warning: Could not write to Excel. Falling back to CSV format.');
            use_excel = false;
        end
    else
        % Use CSV from the start
        disp(['Creating CSV file: ', output_csv]);
    end
        
        if ~use_excel
            % Write headers to CSV
            csvHeader = cell2table(headers);
            writetable(csvHeader, output_csv, 'WriteMode', 'overwrite');
        end
        
        row = 2;
        current_image = '';
        
        for i = 1:length(data)
            line = data{i};
            
            % Check if this is an image header line
            if contains(line, 'Image:')
                parts = strsplit(line, 'Image: ');
                if length(parts) >= 2
                    current_image = strtrim(parts{2});
                end
                continue;
            end
            
            % Parse parameter and value
            if contains(line, ':')
                parts = strsplit(line, ':');
                if length(parts) >= 2
                    parameter = strtrim(parts{1});
                    value_str = strtrim(parts{2});
                    
                    % Try to convert value to number if possible
                    value = str2double(value_str);
                    if isnan(value)
                        value = value_str;
                    end
                    
                    % Write to Excel or CSV
                    row_data = {current_image, parameter, value};
                    
                    if use_excel
                        % Try writing to Excel without showing warnings
                        try
                            warning('off', 'MATLAB:xlswrite:AddSheet');
                            if isnumeric(value)
                                xlswrite(output_excel, {current_image}, 'Analysis', ['A' num2str(row)]);
                                xlswrite(output_excel, {parameter}, 'Analysis', ['B' num2str(row)]);
                                xlswrite(output_excel, value, 'Analysis', ['C' num2str(row)]);
                            else
                                xlswrite(output_excel, {current_image}, 'Analysis', ['A' num2str(row)]);
                                xlswrite(output_excel, {parameter}, 'Analysis', ['B' num2str(row)]);
                                xlswrite(output_excel, {value}, 'Analysis', ['C' num2str(row)]);
                            end
                        catch ME
                            if row == 2 % Only show warning once on the first row
                                disp(['Error writing to Excel: ', ME.message]);
                                disp('Continuing with CSV only.');
                                use_excel = false;
                            end
                        end
                    end
                    
                    if ~use_excel
                        % Append to CSV
                        csvRow = cell2table(row_data);
                        writetable(csvRow, output_csv, 'WriteMode', 'append');
                    end
                    
                    row = row + 1;
                end
            end
        end
        
        disp(['Wrote ', num2str(row-2), ' data rows.']);
        
        if use_excel && exist(output_excel, 'file')
            disp('Excel data writing completed.');
        else
            disp(['CSV data written to: ', output_csv]);
        end
    catch ME
        disp(['Error writing data: ', ME.message]);
        disp('Creating a basic CSV file instead.');
        
        % Create a simple CSV file
        basic_csv = [direct, 'analysis_results_basic.csv'];
        fid_csv = fopen(basic_csv, 'w');
        fprintf(fid_csv, 'Image,Parameter,Value\n');
        
        for i = 1:length(data)
            fprintf(fid_csv, '%s\n', data{i});
        end
        
        fclose(fid_csv);
        disp(['Basic CSV file created at: ', basic_csv]);
    end

disp('Data writing process completed.');
end
