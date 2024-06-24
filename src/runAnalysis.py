import matlab.engine

eng = matlab.engine.start_matlab()
eng.c_analysis(nargout=0)