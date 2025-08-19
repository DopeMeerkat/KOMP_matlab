import matlab.engine
import os
import sys

def main():
    """
    Call MATLAB function to write analysis data from stack images to Excel.
    """
    print("Starting writeDataStack.py...")
    
    try:
        # Get the absolute path to the src directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        print(f"Current directory: {current_dir}")
        
        # Start MATLAB engine
        print("Starting MATLAB engine...")
        eng = matlab.engine.start_matlab()
        
        # Add the source directory to the MATLAB path
        eng.addpath(current_dir, nargout=0)
        
        # Call the stack-specific data writing function
        print("Calling control_data_write_stack...")
        eng.control_data_write_stack(nargout=0)
        
        # Check if the output files exist
        home_dir = os.path.dirname(current_dir)
        output_dir = os.path.join(home_dir, 'CCC_E01_hF_FL1_stack_output', 'Analyzed')
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
        print(f"Error in writeDataStack.py: {e}")
        sys.exit(1)
    
    sys.exit(0)

if __name__ == "__main__":
    main()
