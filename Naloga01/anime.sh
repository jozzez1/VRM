# anime.sh
##########

#!/bin/zsh

CWD=`pwd`
DIR=`pwd`/$1
ANI=$DIR/animation
COUNT=0

NUMBER=`ls $DIR | wc -l`
FPS=`echo "$NUMBER/12" | bc`

echo $FPS && sleep 2

mkdir $ANI
cd $ANI

for FILE in `ls $DIR | grep txt`;
do
	echo "set terminal jpeg size 600,400;
	  set tics nomirror;
	  set y2tics;
	  set grid;
	  set label \"$COUNT\" at screen 0.20,0.87;
	  set ylabel \"|u(x)|^2\";
	  set y2label \"V(x)\";
	  set xlabel \"x\";
	  set yrange [0:1];
	  set output \"$FILE.jpeg\";
	  plot \"$DIR/$FILE\" u 1:2 w l lt 1 title \"|u(x)|^2\" axes x1y1, \
		 \"$DIR/$FILE\" u 1:3 w l lt 3 title \"V\" axes x1y2 " | gnuplot

	  (( COUNT ++ ))
done

mencoder "mf://*.jpeg" -mf w=600:h=400:fps=${FPS}:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $1.avi

rm *.jpeg
cp $1.avi $CWD
rm -r $DIR

exit 0
