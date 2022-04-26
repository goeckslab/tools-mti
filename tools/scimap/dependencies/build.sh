#!/bin/bash

# Install sklearn
if [ `uname` == Darwin ]; then
    $PYTHON -m pip install https://files.pythonhosted.org/packages/04/79/306714443069278702eaf767e82c0a619d7e93119fa9d7b075ccc7dd4801/scikit_learn-0.24.2-cp37-cp37m-macosx_10_13_x86_64.whl 
fi

if [ `uname` == Linux ]; then
    $PYTHON -m pip install https://files.pythonhosted.org/packages/6f/6b/10881b09340d69d4a941e5624bfbf1ba853be8cdf2141077e66dda0b088e/scikit_learn-0.24.2-cp37-cp37m-manylinux1_x86_64.whl 
fi

# Install phenograph
$PYTHON -m pip install https://files.pythonhosted.org/packages/fc/37/4aa1d8c2ded0c612031d32ad8606b6222243f9326ca28754122e306680be/PhenoGraph-1.5.7-py3-none-any.whl

# Install scimap
$PYTHON -m pip install . -vv
