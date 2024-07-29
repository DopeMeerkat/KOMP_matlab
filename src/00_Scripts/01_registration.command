#!/bin/bash
cd "$(dirname "$0")"
cd ..

python3 runRegistration.py
python3 filterRegistration.py