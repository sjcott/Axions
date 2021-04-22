#!/bin/bash

# Testing setting up directories on centaurus and transfering files

# Get directory of script
VAR1="$(dirname "$(readlink -f "$0")")"

# Copy required files to centaurus
#scp $VAR1/ic.txt sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/ic.txt

#scp $VAR1/evolution.cpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/evolution.cpp

scp $VAR1/Data/SOR_Fields.txt sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/Data/SOR_Fields.txt

scp $VAR1/ic_w_evolution.cpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/ic_w_evolution.cpp

scp $VAR1/array.cpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/array.cpp

scp $VAR1/array.hpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/array.hpp

scp $VAR1/ic_func.cpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/ic_func.cpp

scp $VAR1/ic_func.hpp sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/ic_func.hpp

scp $VAR1/submit_job.sh sjcott@centaurus.jb.man.ac.uk:/home/sjcott/Documents/Axions/submit_job.sh



