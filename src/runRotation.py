import matlab.engine

eng = matlab.engine.start_matlab()
eng.control_rotate(nargout=0)