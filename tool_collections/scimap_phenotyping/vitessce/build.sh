#!/bin/bash

# Install generate-tiff-offsets
$PYTHON -m pip install https://files.pythonhosted.org/packages/a3/32/0b489c4d19e5b2cd06abbbdcbc0b0a330574d6d50fa024e188928e7a6f85/generate_tiff_offsets-0.1.7-py2.py3-none-any.whl

# Install negspy
$PYTHON -m pip install https://files.pythonhosted.org/packages/75/d1/013c9acf35b77b7eb620f7b93c8d435e2546e1a7427421b5aa6e8acddf44/negspy-0.2.24.tar.gz

# Install Starlette
$PYTHON -m pip install https://files.pythonhosted.org/packages/3b/48/c305e580e6584d8dd0c2c58238dac973f484345d9de4bc1aa5b162c86a54/starlette-0.14.0-py3-none-any.whl

# Install vitessce
$PYTHON -m pip install . -vv
#$PYTHON -m pip install https://files.pythonhosted.org/packages/cd/d5/b3ebf95dded4e10cbc0168226e854de86054a3b9bc1a53ce71ec03130180/vitessce-0.1.0a11-py2.py3-none-any.whl

# Install bftools
cp -r bftools/* $PREFIX/bin/
