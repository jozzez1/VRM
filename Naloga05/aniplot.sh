#!/usr/local/bin/zsh

# we set up some variables

animate ()
{
	CWD=`pwd`

	base=$1
	length=$2
	count=0
	ani=$CWD/$base/animation

	number=`ls $CWD/$base | wc -l`
	fps=`echo "$number/$length" | bc`

	mkdir -p $ani
	cd $ani

	for file in `ls $CWD/$base | grep txt`
	do
		echo -e "set terminal jpeg size 600,400;
		set output \"$file.jpg\";
		set grid;
		set label \"$count\" at screen 0.20,0.87;
		plot \"$CWD/$base/$file\" u 1:2 w l lt 1 title 'myplot'" | gnuplot

		(( count ++ ))
	done

	mencoder "mf://*.jpg" -mf w=600:h=400:fps=${fps}:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $base.avi

	mv $base.avi $cwd
	rm -rf $ani

	cd $CWD
}

# we animate both T and J
animate $1 $3
animate $2 $3

if [ "$4" != "1" ]; then
	rm -rf $1.avi
	rm -rf $2.avi
fi

exit 0

