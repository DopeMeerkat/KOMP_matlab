import matlab.engine
import os

# Get the absolute path to the src directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Start MATLAB engine
eng = matlab.engine.start_matlab()

# Add the source directory to the MATLAB path
eng.addpath(current_dir, nargout=0)

# Call the analysis function for the stack data
eng.c_analysis_stack(nargout=0)
