#!/bin/sh

find $1 -name "outfilecf*" -type f | sort > filelist.pipi.in
find $1 -name "outfilekkcf*" -type f | sort > filelist.kk.in
find $1 -name "outfileppcf*" -type f | sort > filelist.pp.in

