import matlab.engine
import os

# Get the absolute path to the src directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Start MATLAB engine
eng = matlab.engine.start_matlab()

# Add the source directory to the MATLAB path
eng.addpath(current_dir, nargout=0)

# Now call the function
eng.find_growth_plate_stack(nargout=0)
