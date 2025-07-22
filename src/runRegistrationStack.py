#!/usr/bin/env python
"""
runRegistrationStack.py

This script runs the visualization functions for the CCC_E01_hF_FL1_stack folder.
It skips the alignment process since the images are already aligned,
but generates the same visualization outputs as the original registration:
- Grid overlays for each image
- Composite visualizations
- RGB composite
- Multi-color overlay
- Edge overlay
- Shift images (_shift.png)
- Shift2_NoDAPI images (_shift2_NoDAPI.png)

Usage:
    python runRegistrationStack.py

Output:
    Creates visualization images in CCC_E01_hF_FL1_stack_output/Registration
"""

import os
import sys
import subprocess
import platform
import time

def run_visualization():
    """
    Run the MATLAB visualization function for the stack images
    """
    # Determine the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Go to parent directory if we're in the src folder
    parent_dir = os.path.dirname(script_dir) if os.path.basename(script_dir) == 'src' else script_dir
    
    # Change to the parent directory to run the MATLAB script
    os.chdir(parent_dir)
    
    # Determine MATLAB command based on platform
    if platform.system() == 'Windows':
        matlab_cmd = 'matlab -nosplash -nodesktop -r "addpath(\'src\'); visualize_stack_registration; exit;"'
    else:
        matlab_cmd = 'matlab -nosplash -nodesktop -r "addpath(\'src\'); visualize_stack_registration; exit;"'
    
    print(f"Running: {matlab_cmd}")
    
    # Run the MATLAB command
    process = subprocess.Popen(
        matlab_cmd, 
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    # Get the output
    stdout, stderr = process.communicate()
    
    # Print the output
    if stdout:
        print("MATLAB Output:")
        print(stdout.decode('utf-8'))
    
    if stderr:
        print("MATLAB Errors:")
        print(stderr.decode('utf-8'))
    
    # Check if visualization was successful
    reg_dir = os.path.join(parent_dir, 'CCC_E01_hF_FL1_stack_output', 'Registration')
    if os.path.exists(reg_dir) and len(os.listdir(reg_dir)) > 0:
        print(f"Registration visualization completed successfully. Results saved in: {reg_dir}")
        
        # List the generated files
        print("\nGenerated visualization files:")
        for f in sorted(os.listdir(reg_dir)):
            print(f"  - {f}")
        
        # Verify we have the expected output files
        shift_files = [f for f in os.listdir(reg_dir) if '_shift.jpg' in f]
        shift2_nodapi_file = 'CCC_E01_hF_FL1_shift2_NoDAPI.jpg'
        
        if not shift_files:
            print("\nWarning: No shift files were generated. Registration process may be incomplete.")
        else:
            print(f"\nGenerated {len(shift_files)} shift files.")
            
        if shift2_nodapi_file not in os.listdir(reg_dir):
            print(f"Warning: {shift2_nodapi_file} was not generated. Registration process may be incomplete.")
        else:
            print(f"Generated NoDAPI overlay: {shift2_nodapi_file}")
    else:
        print("Registration visualization may have failed. Check the MATLAB output.")

    return process.returncode

def check_prerequisites():
    """
    Check if the stack folder exists and contains the required images
    """
    # Determine the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Go to parent directory if we're in the src folder
    parent_dir = os.path.dirname(script_dir) if os.path.basename(script_dir) == 'src' else script_dir
    
    # Check if the stack folder exists
    stack_dir = os.path.join(parent_dir, 'CCC_E01_hF_FL1_stack')
    if not os.path.exists(stack_dir):
        print(f"Error: Stack folder not found: {stack_dir}")
        return False
    
    # Check if there are PNG images in the stack folder
    png_files = [f for f in os.listdir(stack_dir) if f.endswith('.png')]
    if not png_files:
        print(f"Error: No PNG images found in stack folder: {stack_dir}")
        return False
    
    # Check for required image types
    has_epi = any('Epi' in f for f in png_files)
    has_dia = any('Dia' in f for f in png_files)
    
    if not (has_epi and has_dia):
        print("Warning: Could not find Epi and Dia images. Visualization may be incomplete.")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.join(parent_dir, 'CCC_E01_hF_FL1_stack_output')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # Create Registration directory if it doesn't exist
    reg_dir = os.path.join(output_dir, 'Registration')
    if not os.path.exists(reg_dir):
        os.makedirs(reg_dir)
        print(f"Created Registration directory: {reg_dir}")
    
    return True

def main():
    """
    Main function to run the registration visualization
    """
    print("=" * 80)
    print("KOMP MATLAB: Registration Visualization for CCC_E01_hF_FL1_stack")
    print("=" * 80)
    print("This script will generate registration visualizations for the pre-aligned images.")
    print("No image alignment will be performed since the images are already aligned.")
    print("Generating _shift.png and _shift2_NoDAPI.png files for compatibility.")
    print("=" * 80)
    print()
    
    start_time = time.time()
    
    # Check prerequisites
    if not check_prerequisites():
        sys.exit(1)
    
    # Run visualization
    ret_code = run_visualization()
    
    # Report execution time
    elapsed_time = time.time() - start_time
    print(f"\nExecution completed in {elapsed_time:.2f} seconds.")
    
    sys.exit(ret_code)

if __name__ == "__main__":
    main()
