# Clean the current directory

rm -f *.dat HK_constructed HK_data  log_file out_file
#
#
# Get user input
echo -n "Which sample do you want to run [enter a number between 1 and 5] ? "
read ns
echo
echo -n "Path for mpirun - Press enter to use `which mpirun` "
read mpr
echo
if [ -z $mpr ]
then
	mpr="mpirun -oversubscribe"
fi
echo
echo "The execution is going to use "`which $mpr`
echo
if [ $ns -le 2 ]
then
	np=1
else
echo -n "How many processors do you want to use? "
read np
fi
#
#
# Copy relevant files for execution
#

hdir=`pwd`

cp -rf $hdir/data/sample_"$ns"/input/* .
if [ $ns -ge 3 ]
then
      gunzip HK_data.gz
fi
#
# Now execute and store the screen output in log_file
#
if [[ ! -f $hdir/MO_IPT.run ]]; then
    echo "The executable doesn't exist! Compile now ..."
fi

date > log_file
$mpr -np "$np" ./MO_IPT.run | tee -a log_file
date  >> log_file
#
#
echo
echo
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "Now a comparison may be made between the output of this run"
echo "and the one saved in data/sample_"$ns"/output."
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
