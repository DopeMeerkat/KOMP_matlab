#!/usr/bin/env python3
"""
This script runs the updated find_growth_plate_stack_for_dia.m function
to visualize the ROI using the mineral image for left/right boundaries.
"""

import os
import sys
import subprocess
import platform

# Get the base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(BASE_DIR, 'src')

# Add the src directory to the path
sys.path.append(SRC_DIR)

def get_matlab_command():
    """Get appropriate MATLAB command based on platform"""
    if platform.system() == 'Windows':
        return 'matlab.exe -nodesktop -nosplash -r'
    else:
        return 'matlab -nodesktop -nosplash -r'

def run_growth_plate_detection():
    """Run the growth plate detection and ROI visualization"""
    # Find the 2_Dia.png image
    stack_dir = os.path.join(BASE_DIR, "CCC_E01_hF_FL1_stack")
    
    dia_img_path = os.path.join(stack_dir, "2_Dia.png")
    
    if not os.path.exists(dia_img_path):
        print(f"Error: Could not find {dia_img_path}")
        return False
    
    # Ensure output directory exists
    output_dir = os.path.join(os.path.dirname(stack_dir), "CCC_E01_hF_FL1_stack_output", "Analyzed")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create MATLAB command
    path_for_matlab = dia_img_path.replace('\\', '/')
    matlab_cmd = f"try; find_growth_plate_stack_for_dia('{path_for_matlab}'); catch ME; disp(getReport(ME)); exit(1); end; exit(0);"
    
    # Run MATLAB command
    print(f"Running find_growth_plate_stack_for_dia with {dia_img_path}...")
    cmd = f"{get_matlab_command()} \"{matlab_cmd}\""
    
    process = subprocess.Popen(cmd, shell=True)
    exit_code = process.wait()
    
    if exit_code != 0:
        print("Error: MATLAB execution failed")
        return False
    
    # Check for visualization outputs
    roi_viz_dir = os.path.join(output_dir, "ROI_Visualizations")
    growth_plate_viz = os.path.join(output_dir, "dia_growth_plate_visualization.png")
    roi_viz = os.path.join(roi_viz_dir, "roi_on_nodapi.jpg")
    
    if os.path.exists(growth_plate_viz):
        print(f"Growth plate visualization created: {growth_plate_viz}")
    
    if os.path.exists(roi_viz):
        print(f"ROI visualization created: {roi_viz}")
        
        # Try to open the visualization (platform-specific)
        try:
            if platform.system() == 'Darwin':  # macOS
                subprocess.run(['open', roi_viz], check=False)
                print("Opened ROI visualization")
            elif platform.system() == 'Windows':
                os.startfile(roi_viz)
                print("Opened ROI visualization")
            elif platform.system() == 'Linux':
                subprocess.run(['xdg-open', roi_viz], check=False)
                print("Opened ROI visualization")
        except Exception as e:
            print(f"Could not open visualization: {e}")
    
    growth_plate_data = os.path.join(output_dir, "growth_plate1.mat")
    if os.path.exists(growth_plate_data):
        print(f"Growth plate data saved to: {growth_plate_data}")
        return True
    else:
        print("Warning: Growth plate data file not created")
        return False

if __name__ == "__main__":
    success = run_growth_plate_detection()
    if success:
        print("ROI visualization completed successfully")
    else:
        print("ROI visualization failed")
        sys.exit(1)
