#!/bin/bash
cd "$(dirname "$0")"
cd ..

python3 repairThreshold.py
python3 runRotation.py
# python3 runGrowthPlate.py