# anime.sh
##########

#!/usr/local/bin/zsh

CWD=`pwd`
DIR=`pwd`/$1
ANI=$DIR/animation
SAVE=$2
LENGTH=$3
dT=$4
COUNT="0"

NUMBER=`ls $DIR | wc -l`
FPS=`echo "$NUMBER/$LENGTH" | bc`

echo $FPS && sleep 2

mkdir $ANI
cd $ANI

for FILE in `ls $DIR | grep txt`;
do
	I="$COUNT"
	T=`echo "scale = 4; $dT * $I + 0.5" | bc`

	echo "T = $T"

	echo "set terminal jpeg enhanced font \"/usr/local/lib/X11/fonts/ChromeOS/Arimo-Regular.ttf,12\" size 550,400;
	  unset surf;
	  set view map;
	  set palette rgbformulae 10,10,10 maxcolor 2;
	  set title \"Temperature = $T\";
	  set output \"$FILE.jpg\";
	  set cbrange [-1.1:1.1];
	  plot \"$DIR/$FILE\" matrix with image title '';" | gnuplot

	(( COUNT ++ ))

done

mencoder -msgcolor "mf://*.jpg" -mf w=600:h=400:fps=${FPS}:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $1.avi

cp $1.avi $CWD
rm *.jpeg

cd $CWD
rm -rf $DIR
mplayer $1.avi

if [ "$SAVE" == "0" ]; then

	rm $1.avi

fi

exit 0

