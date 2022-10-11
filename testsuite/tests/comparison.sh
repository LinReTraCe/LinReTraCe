#! /bin/bash

./lprint testsuite/tests/Si_output.hdf5 c-intra xx s-intra xx > testsuite/tests/Si_output.txt
python3 testsuite/tests/diff.py
if [ $? -eq 0 ]; then
  echo "-- Transport Output Comparison successful."
else
  echo "-- Transport Output Comparison FAILED."
fi
