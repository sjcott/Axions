#!/bin/bash

ssh -o "ProxyCommand ssh -A sjcott@external.jb.man.ac.uk -W %h:%p" sjcott@centaurus 'cd Documents/Axions/; sbatch submit_job.sh'

