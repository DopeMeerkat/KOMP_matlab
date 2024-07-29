#!/bin/bash
cd "$(dirname "$0")"
cd ..

#path/to/runRegistration.py
python3 runRegistration.py
python3 filterRegistration.py