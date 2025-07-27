#!/bin/bash

# copy file at ARG $1 in the remote server at student.ph.ed.ac.uk to the location at ARG $2

scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -r "s2696227@student.ph.ed.ac.uk:~/Documents/Dissertation/1DActiveSolids/$1" "$2"
