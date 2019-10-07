bindir=$HOME/bin

mkdir -p $bindir
mkdir -p $bindir/share/kabatman

(cd src; make)
cp src/kabatman $bindir
cp kabattest/kabattest.pl $bindir
if [ -x $bindir/kabattest ]; then
    echo "kabattest already exists in $bindir. Do you wish to remove it?"
    rm -i $bindir/kabattest
fi
(cd $bindir; ln -s kabattest.pl kabattest)
cp data/* $bindir/share/kabatman
export KABATDIR=$bindir/share/kabatman

if grep -q KABATDIR $HOME/.bashrc; then
    echo "Your ~/.bashrc already contains a setting for KABATDIR so not changing it."
else
    echo "export KABATDIR=$KABATDIR" >>$HOME/.bashrc
    echo "The following line has been added to your ~/.bashrc file:"
    echo "   export KABATDIR=$KABATDIR"
fi

