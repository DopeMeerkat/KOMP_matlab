import matlab.engine

eng = matlab.engine.start_matlab()
eng.find_registration_errors_4_UCHC_TRAP2(nargout=0)
eng.find_registration_errors_4_UCHC_AP2(nargout=0)
eng.c_thresholding(nargout=0)

#eng.control_rotate.m(nargout=0)