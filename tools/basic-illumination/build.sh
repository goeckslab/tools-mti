#!/bin/bash

sharedir=$PREFIX/share
mkdir -p $PREFIX/bin

# Install Fiji
unamestr=`uname`
if [ "$unamestr" == 'Linux' ];
then
	wget https://downloads.imagej.net/fiji/latest/fiji-linux64.zip
	unzip fiji-linux64.zip
	cp -R Fiji.app/* $sharedir/
    ln -s $sharedir/ImageJ-linux64 $PREFIX/bin/ImageJ

elif [ "$unamestr" == 'Darwin' ];
then
	wget https://downloads.imagej.net/fiji/latest/fiji-macosx.zip
	unzip fiji-macosx.zip
	cp -R Fiji.app/* $sharedir/
	ln -s $sharedir/Contents/MacOS/ImageJ-macosx $PREFIX/bin/ImageJ
fi

chmod 0755 "${PREFIX}/bin/ImageJ"

# Install BaSiC
cp BaSiCPlugin/BaSiC_.jar $sharedir/plugins/
cp BaSiCPlugin/Dependent/*.jar $sharedir/jars/
rm $sharedir/jars/jtransforms-2.4.jar