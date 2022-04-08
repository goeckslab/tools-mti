#!/bin/bash

set -x

pip install https://files.pythonhosted.org/packages/45/3f/e84aef209001eb8dba16427855b0dc39f2a58905d7de53a7dd7187059bb2/mxnet_mkl-1.6.0-py2.py3-none-manylinux1_x86_64.whl


ln -sf $SP_DIR/mxnet/libgfortran.so.3  $PREFIX/lib/libgfortran.so.3 
ln -sf $SP_DIR/mxnet/libiomp5.so  $PREFIX/lib/libiomp5.so
ln -sf $SP_DIR/mxnet/libmkldnn.so.0  $PREFIX/lib/libmkldnn.so.0
ln -sf $SP_DIR/mxnet/libmkldnn.so.1  $PREFIX/lib/libmkldnn.so.1
ln -sf $SP_DIR/mxnet/libmklml_intel.so  $PREFIX/lib/libmklml_intel.so
ln -sf $SP_DIR/mxnet/libmxnet.so  $PREFIX/lib/libmxnet.so
ln -sf $SP_DIR/mxnet/libquadmath.so.0  $PREFIX/lib/libquadmath.so.0
