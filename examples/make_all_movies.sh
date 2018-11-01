#
# Small script to run all examples and produce movies
#

CELADRO=../build/celadro
THREADS=4

set -e  # stop on any error
for runcard in $(find */*.dat); do
  directory=${runcard%/*}
  $CELADRO $runcard -fco $directory/output -t$THREADS
  python2 $directory/plot.py $directory/output $directory
done
