#! /bin/bash

# M.Pickem: mini help script to delete all files created via the setup.py command
# sometimes necessary if things need to be cleaned up properly

python3 setup.py install --record files.txt
cat files.txt | xargs rm
rm files.txt
