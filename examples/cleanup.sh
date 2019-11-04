#
# Cleanup all output data and movies.
#

set -e  # stop on any error
for runcard in $(find */*.dat); do
  directory=${runcard%/*}
  rm -rf $directory/data/ movie_$directory.mp4
done
