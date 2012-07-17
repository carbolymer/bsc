#!/bin/sh

rm -f *.data

IFS_BAK=$IFS
IFS="
"
# Kaons
find $1 -name "outfilekkcf*" -type f | sort | ./merger 1> filelist.kk.in
# FILES=`cat filelist.kk.in`
# for parameter in $FILES
# do
# 	IFS=$IFS_BAK
# 	./fitsh $parameter kk &
# 	./fit1d $parameter kk &
# 	IFS="
# "
# done

# Pions
find $1 -name "outfilecf*" -type f | sort | ./merger 1> filelist.pipi.in
# FILES=`cat filelist.pipi.in`
# for parameter in $FILES
# do
# 	IFS=$IFS_BAK
# 	./fitsh $parameter pipi &
# 	./fit1d $parameter pipi &
# 	IFS="
# "
# done

# Protons
find $1 -name "outfileppcf*" -type f | sort | ./merger 1> filelist.pp.in
# FILES=`cat filelist.pp.in`
# for parameter in $FILES
# do
# 	IFS=$IFS_BAK
# 	./fitsh $parameter pp &
# 	./fit1d $parameter pp &
# 	IFS="
# "
# done

# IFS=$IFS_BAK
# IFS_BAK=