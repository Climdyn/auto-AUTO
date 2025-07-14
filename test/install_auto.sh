#!/bin/bash
# Script to install AUTO on Github runners

cd $1
git clone https://github.com/auto-07p/auto-07p.git
cd auto-07p
./configure
make
