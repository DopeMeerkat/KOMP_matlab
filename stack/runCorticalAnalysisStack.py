#!/usr/bin/env python3
"""
Run cortical analysis on stack images using the 2_Dia.png for growth plate detection 
and the NoDAPI overlay from the Registration folder for cortical region detection.
"""

import matlab.engine
import os
import sys
import argparse
import subprocess

def main():
    """
    Run cortical analysis on stack images, using the growth plate data from 2_Dia.png
    and the registered NoDAPI image.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run cortical analysis on stack images')
    parser.add_argument('--visualize', action='store_true', 
                        help='Open visualization images after analysis')
    args = parser.parse_args()
    
    try:
        # Get the absolute path to the src directory and home directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        home_dir = os.path.dirname(current_dir)
        
        print(f"Current directory: {current_dir}")
        
        # Define paths
        stack_dir = os.path.join(home_dir, 'CCC_E01_hF_FL1_stack')
        output_dir = os.path.join(home_dir, 'CCC_E01_hF_FL1_stack_output')
        registration_dir = os.path.join(output_dir, 'Registration')
        analyzed_dir = os.path.join(output_dir, 'Analyzed')
        growth_plate_file = os.path.join(analyzed_dir, 'growth_plate1.mat')
        
        # Check if growth plate data exists
        if not os.path.exists(growth_plate_file):
            print(f"Growth plate data not found at: {growth_plate_file}")
            print("Running growth plate detection first...")
            
            # Check if Dia image exists
            dia_img = os.path.join(stack_dir, '2_Dia.png')
            if not os.path.exists(dia_img):
                print(f"Error: 2_Dia.png not found at {dia_img}")
                sys.exit(1)
                
            # Start MATLAB engine
            print("Starting MATLAB engine...")
            eng = matlab.engine.start_matlab()
            
            # Add the source directory to the MATLAB path
            eng.addpath(current_dir, nargout=0)
            
            # Create find_growth_plate_stack_for_dia function to detect growth plate from 2_Dia.png
            print("Running growth plate detection on 2_Dia.png...")
            eng.eval("find_growth_plate_stack_for_dia('" + dia_img.replace('\\', '/') + "')", nargout=0)
            
            # Check if growth plate data was created
            if not os.path.exists(growth_plate_file):
                print(f"Error: Failed to create growth plate data at {growth_plate_file}")
                sys.exit(1)
        
        # Check if the NoDAPI overlay image exists
        nodapi_img = os.path.join(registration_dir, 'CCC_E01_hF_FL1_shift2_NoDAPI.jpg')
        if not os.path.exists(nodapi_img):
            print(f"Warning: NoDAPI overlay image not found at {nodapi_img}")
            print("Checking for other overlay images in Registration folder...")
            
            # Try to find any overlay image
            if os.path.exists(registration_dir):
                overlay_files = [f for f in os.listdir(registration_dir) if f.endswith('.jpg') or f.endswith('.png')]
                if overlay_files:
                    overlay_img = os.path.join(registration_dir, overlay_files[0])
                    print(f"Found alternative overlay image: {overlay_files[0]}")
                    nodapi_img = overlay_img
                else:
                    print(f"Error: No overlay images found in {registration_dir}")
                    print("Please run the registration visualization script first.")
                    sys.exit(1)
            else:
                print(f"Error: Registration directory not found at {registration_dir}")
                print("Please run the registration visualization script first.")
                sys.exit(1)
        
        # Start MATLAB engine if not already started
        try:
            eng
        except NameError:
            print("Starting MATLAB engine...")
            eng = matlab.engine.start_matlab()
            eng.addpath(current_dir, nargout=0)
        
        # Load the growth plate data to check if it has the left/right boundaries
        print("Loading growth plate data...")
        growth_plate_data = eng.load(growth_plate_file.replace('\\', '/'), nargout=1)
        
        # Check if growth plate data has left and right boundary info
        has_horizontal_boundaries = False
        try:
            # Check if growth_plate1 cell array has left and right boundary indices
            if len(growth_plate_data['growth_plate1'][0][0]) >= 5:
                left_x = growth_plate_data['growth_plate1'][0][0][3][0][0]
                right_x = growth_plate_data['growth_plate1'][0][0][4][0][0]
                has_horizontal_boundaries = True
                print(f"Found horizontal boundaries in growth plate data: left_x={left_x}, right_x={right_x}")
            
            # Check if bt1 struct exists and has left/right fields
            if 'bt1' in growth_plate_data and hasattr(growth_plate_data['bt1'], 'left_x') and hasattr(growth_plate_data['bt1'], 'right_x'):
                left_x = growth_plate_data['bt1'].left_x
                right_x = growth_plate_data['bt1'].right_x
                has_horizontal_boundaries = True
                print(f"Found horizontal boundaries in bt1 struct: left_x={left_x}, right_x={right_x}")
        except Exception as e:
            print(f"Warning: Error accessing horizontal boundaries in growth plate data: {e}")
            
        # Run cortical analysis - the function handles saving the output internally
        print("Running cortical analysis...")
        eng.find_cortical_stack(
            nodapi_img.replace('\\', '/'), 
            growth_plate_data,
            nargout=0  # No need to return anything, files are saved by the MATLAB function
        )
        
        # Create output directory for cortical results if it doesn't exist
        cortical_dir = os.path.join(analyzed_dir, 'Cortical')
        if not os.path.exists(cortical_dir):
            os.makedirs(cortical_dir)
        
        print(f"Cortical analysis completed. Results saved to: {cortical_dir}")
        
        # Open visualization if requested
        viz_file = os.path.join(cortical_dir, 'CCC_E01_hF_FL1_shift2_NoDAPI_cortical_analysis.png')
        if args.visualize and os.path.exists(viz_file) and sys.platform == 'darwin':
            try:
                subprocess.run(['open', viz_file], check=False)
                print(f"Opened visualization: {viz_file}")
            except:
                print(f"Could not open visualization. File is located at: {viz_file}")
    
    except Exception as e:
        print(f"Error in runCorticalAnalysisStack.py: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        # Close MATLAB engine if it was started
        try:
            eng.quit()
            print("MATLAB engine closed.")
        except:
            pass
    
    print("Cortical analysis workflow completed.")
    sys.exit(0)

if __name__ == "__main__":
    main()
