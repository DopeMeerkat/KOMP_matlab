#!/bin/bash
cd "$(dirname "$0")"
cd ..

python3 runAnalysis.py
python3 writeData.py