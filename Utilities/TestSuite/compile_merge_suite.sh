#!/bin/sh

python3.6 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt

source ./initialize_cluster.sh
python3 run_testsuite.py ../Inputfiles/TestSuiteConfiguration/.merge-compile.xml

deactivate
rm -rf venv
