#!/bin/sh

CENTRALITY=b5

EVENT_DIR=/mnt/store/mgalazyn/lhyqid3v_LHCPbPb_2760_$CENTRALITY/

echo "Running for centrality $CENTRALITY"

rm -f data/$CENTRALITY/kk*.out
rm -f data/$CENTRALITY/pipi*.out
rm -f data/$CENTRALITY/pp*.out
rm -f data/$CENTRALITY
rm -f log/*.log

IFS_BAK=$IFS
IFS="
"

# Kaons
find $1 -name "outfilekkcf*" -type f | sort | ./merger 1> filelist.kk.in
find $EVENT_DIR -name "outfilekkcf*" -type f | sort | ./merger 1> filelist.kk.in
FILES=`cat filelist.kk.in`
for parameter in $FILES
do
	IFS=$IFS_BAK
	./fitsh $parameter data/$CENTRALITY/kk &>> log/kk_fitsh.log &
	./fit1d $parameter data/$CENTRALITY/kk &>> log/kk_fit1d.log &
	IFS="
"
done

# # Pions
find $1 -name "outfilecf*" -type f | sort | ./merger 1> data/filelist.pipi.in
find $EVENT_DIR -name "outfilecf*" -type f | sort | ./merger 1> data/filelist.pipi.in
FILES=`cat filelist.pipi.in`
for parameter in $FILES
do
	IFS=$IFS_BAK
	./fitsh $parameter data/$CENTRALITY/pipi &>> log/pipi_fitsh.log &
	./fit1d $parameter data/$CENTRALITY/pipi &>> log/pipi_fit1d.log
	IFS="
"
done

# Protons
find $1 -name "outfileppcf*" -type f | sort | ./merger 1> data/filelist.pp.in
find $EVENT_DIR -name "outfileppcf*" -type f | sort | ./merger 1> data/filelist.pp.in
FILES=`cat filelist.pp.in`
for parameter in $FILES
do
	IFS=$IFS_BAK
	./fitsh $parameter data/$CENTRALITY/pp &>> log/pp_fitsh.log &
	./fit1d $parameter data/$CENTRALITY/pp &>> log/pp_fit1d.log &
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
