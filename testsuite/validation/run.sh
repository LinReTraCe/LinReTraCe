#! /bin/bash

./test

if [ $? -eq 0 ]; then
  echo "-- HDF5 validation test successful."
else
  echo "-- HDF5 validation test FAILED."
  echo "   If the compilation was successful and the"
  echo "   HDF5 library links were not found (e.g. cannot open shared object file ...)"
  echo "   make sure to include the HDF5 library in the LD_LIBRARY_PATH"
  echo "   export LD_LIBRARY_PATH=/path/to/hdf5/lib:\$LD_LIBRARY_PATH"
fi
