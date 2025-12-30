#!/bin/bash
make DISABLED_MARCH_NATIVE=1
mkdir -p $PREFIX/bin
cp ./halign4 $PREFIX/bin/
