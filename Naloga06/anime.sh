#!/bin/sh

CWD=`pwd`
DIR=`pwd`/$1
SAVE=$2
LENGTH=$3

NUMBER=`ls $DIR | wc -l`
FPS=`echo "$NUMBER/$LENGTH" | bc`

cd $DIR

mencoder -msgcolor "mf://*.jpg" -mf w=600:h=400:fps=${FPS}:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $1.avi

cp $1.avi $CWD
#rm *.jpeg

cd $CWD
#rm -rf $DIR
mplayer $1.avi

if [ "$SAVE" == "0" ]; then
	rm $1.avi
fi

exit 0

