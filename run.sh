#!/bin/sh

CENTRALITY=b5

EVENT_DIR=/mnt/store/mgalazyn/lhyqid3v_LHCPbPb_2760_$CENTRALITY/

echo "Running for centrality $CENTRALITY"

rm -f $CENTRALITY/pp*.out
rm -f $CENTRALITY/kk*.out
rm -f *.log

IFS_BAK=$IFS
IFS="
"

# Kaons
# find $1 -name "outfilekkcf*" -type f | sort | ./merger 1> filelist.kk.in
# find $EVENT_DIR -name "outfilekkcf*" -type f | sort | ./merger 1> filelist.kk.in
FILES=`cat filelist.kk.in`
for parameter in $FILES
do
	IFS=$IFS_BAK
	./fitsh $parameter $CENTRALITY/kk &>> fitsh.log &
	./fit1d $parameter $CENTRALITY/kk &>> fit1d.log &
	IFS="
"
done

# # Pions
# find $1 -name "outfilecf*" -type f | sort | ./merger 1> filelist.pipi.in
# find $EVENT_DIR -name "outfilecf*" -type f | sort | ./merger 1> filelist.pipi.in
# FILES=`cat filelist.pipi.in`
# for parameter in $FILES
# do
# 	IFS=$IFS_BAK
# 	./fitsh $parameter $CENTRALITY/pipi &>> fitsh.log &
# 	./fit1d $parameter $CENTRALITY/pipi &
# 	IFS="
# "
# done

# Protons
# find $1 -name "outfileppcf*" -type f | sort | ./merger 1> filelist.pp.in
# find $EVENT_DIR -name "outfileppcf*" -type f | sort | ./merger 1> filelist.pp.in
FILES=`cat filelist.pp.in`
for parameter in $FILES
do
	IFS=$IFS_BAK
	./fitsh $parameter $CENTRALITY/pp &>> fitsh.log &
	./fit1d $parameter $CENTRALITY/pp &>> fit1d.log &
	IFS="
"
done



IFS=$IFS_BAK
IFS_BAK=

echo "Waiting for all fitting processes..."
while [ `ps aux | grep [f]it | wc -l` != 0 ]; do
	sleep 0.3
done
echo "Plotting..."
./plotter
echo -e "\n\nDONE"
