import matlab.engine

eng = matlab.engine.start_matlab()
eng.control_data_write_2(nargout=0)