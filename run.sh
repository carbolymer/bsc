#!/bin/sh

CENTRALITY=bb3m6

EVENT_DIR=/home/mgalazyn/workspace/tpi_output/bb3m6

echo "Running for centrality $CENTRALITY"

rm -f data/$CENTRALITY/kk*.out
rm -f data/$CENTRALITY/pipi*.out
rm -f data/$CENTRALITY/pp*.out
rm -f log/*.log

if [ "$1" != "" ]; then
	EVENT_DIR=$1
fi

if [ "$2" != "" ]; then
	CENTRALITY=$2
fi

if [ ! -d data/$CENTRALITY ] ; then
	mkdir data/$CENTRALITY
fi


IFS_BAK=$IFS
IFS="
"

# Kaons
find $EVENT_DIR -name "outfilekkcf*" -type f | sort | ./merger 1> data/$CENTRALITY/filelist.kk.in
FILES=`cat data/$CENTRALITY/filelist.kk.in`
for parameter in $FILES
do
	IFS=$IFS_BAK
	./fitsh $parameter data/$CENTRALITY/kk &>> log/kk_fitsh.log &
	./fit1d $parameter data/$CENTRALITY/kk &>> log/kk_fit1d.log &
	IFS="
"
done

# # Pions
find $EVENT_DIR -name "outfilecf*" -type f | sort | ./merger 1> data/$CENTRALITY/filelist.pipi.in
FILES=`cat data/$CENTRALITY/filelist.pipi.in`
for parameter in $FILES
do
	IFS=$IFS_BAK
	./fitsh $parameter data/$CENTRALITY/pipi &>> log/pipi_fitsh.log &
	./fit1d $parameter data/$CENTRALITY/pipi &>> log/pipi_fit1d.log &
	IFS="
"
done

# Protons
find $EVENT_DIR -name "outfileppcf*" -type f | sort | ./merger 1> data/$CENTRALITY/filelist.pp.in
FILES=`cat data/$CENTRALITY/filelist.pp.in`
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
# make plots
echo -e "\n\nDONE"
