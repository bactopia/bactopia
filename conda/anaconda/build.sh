#!/bin/bash
set -e -x

BACTOPIA_HOME="$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM"
mkdir -p $BACTOPIA_HOME
mkdir -p ${PREFIX}/bin

cp -R ./* ${BACTOPIA_HOME}
chmod 777 ${BACTOPIA_HOME}/bactopia.nf
chmod 777 ${BACTOPIA_HOME}/bin/*
ln -s ${BACTOPIA_HOME}/bactopia.nf ${PREFIX}/bin/bactopia
ln -s ${BACTOPIA_HOME}/bin/setup-datasets.py ${PREFIX}/bin/setup-datasets
ln -s ${BACTOPIA_HOME}/bin/prepare-fofn.py ${PREFIX}/bin/prepare-fofn
