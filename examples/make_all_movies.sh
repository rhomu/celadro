#
# Small script to run all examples and produce movies
#

CELADRO=../build/celadro
THREADS=16

set -e  # stop on any error
for runcard in $(find __*/*.dat); do
  directory=${runcard%/*}
  $CELADRO $runcard -fco $directory/output -t$THREADS
  python3 $directory/plot.py $directory/output $directory
done
