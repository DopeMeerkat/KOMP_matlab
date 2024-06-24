import matlab.engine

eng = matlab.engine.start_matlab()
eng.find_growth_plate_2(nargout=0)