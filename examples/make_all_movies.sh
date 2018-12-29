#
# Small script to run all examples and produce movies
#

CELADRO=../build/celadro
THREADS=4

set -e  # stop on any error
for runcard in $(find */config.dat); do
  directory=${runcard%/*}
  $CELADRO $directory -fct$THREADS
  python3 $directory/plot.py $directory/data $directory
done
