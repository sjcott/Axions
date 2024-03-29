#!/bin/bash

# Testing setting up directories on centaurus and transfering files

# Get directory of script
VAR1="$(dirname "$(readlink -f "$0")")"

# Copy required files to centaurus
scp $VAR1/moore_ic.txt sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/moore_ic.txt

scp $VAR1/mooreEvol.cpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/mooreEvol.cpp

scp $VAR1/header.cpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/header.cpp

scp $VAR1/header.hpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/header.hpp

scp $VAR1/submit_job.sh sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/submit_job.sh



