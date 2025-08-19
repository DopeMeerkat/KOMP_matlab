import matlab.engine
import os
import sys
import argparse
import subprocess

def main():
    """
    Run analysis on stack images and optionally write data to Excel
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run analysis on stack images')
    parser.add_argument('--write-data', action='store_true', 
                        help='Write analysis data to Excel after analysis')
    args = parser.parse_args()
    
    try:
        # Get the absolute path to the src directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        print(f"Current directory: {current_dir}")
        
        # Start MATLAB engine
        print("Starting MATLAB engine...")
        eng = matlab.engine.start_matlab()
        
        # Add the source directory to the MATLAB path
        eng.addpath(current_dir, nargout=0)
        
        # Call the analysis function for the stack data
        print("Running analysis on stack images...")
        eng.c_analysis_stack(nargout=0)
        
        # Check if analysis was successful
        home_dir = os.path.dirname(current_dir)
        output_dir = os.path.join(home_dir, 'CCC_E01_hF_FL1_stack_output', 'Analyzed')
        
        if os.path.exists(output_dir) and len(os.listdir(output_dir)) > 0:
            print(f"Analysis completed successfully. Results in: {output_dir}")
            
            # Check for ROI visualizations
            roi_vis_dir = os.path.join(output_dir, 'ROI_Visualizations')
            roi_report = os.path.join(output_dir, 'roi_analysis_report.html')
            
            if os.path.exists(roi_vis_dir) and len(os.listdir(roi_vis_dir)) > 0:
                print(f"ROI visualizations available in: {roi_vis_dir}")
                
                if os.path.exists(roi_report):
                    print(f"ROI analysis HTML report available: {roi_report}")
                    
                    # On macOS, try to open the report in the default browser
                    if sys.platform == 'darwin':
                        try:
                            subprocess.run(['open', roi_report], check=False)
                            print("Opened ROI analysis report in browser")
                        except:
                            print("Could not automatically open report. Please open it manually.")
        else:
            print("Warning: Analysis output directory empty or not found.")
        
        # If --write-data flag is provided, also run the data writing process
        if args.write_data:
            print("Writing analysis data to Excel...")
            eng.control_data_write_stack(nargout=0)
            
            # Check if the output files exist
            excel_file = os.path.join(output_dir, 'analysis_results.xlsx')
            csv_file = os.path.join(output_dir, 'analysis_results.csv')
            basic_csv_file = os.path.join(output_dir, 'analysis_results_basic.csv')
            
            if os.path.exists(excel_file) and os.path.getsize(excel_file) > 0:
                print(f"Excel file successfully created: {excel_file}")
            elif os.path.exists(csv_file) and os.path.getsize(csv_file) > 0:
                print(f"CSV file created: {csv_file}")
            elif os.path.exists(basic_csv_file) and os.path.getsize(basic_csv_file) > 0:
                print(f"Basic CSV file created: {basic_csv_file}")
            else:
                print("Warning: No output files found. Check MATLAB output for errors.")
            
            print("Data writing completed successfully.")
    
    except Exception as e:
        print(f"Error in runAnalysisStack.py: {e}")
        sys.exit(1)
    
    print("Analysis workflow completed.")
    sys.exit(0)

if __name__ == "__main__":
    main()
