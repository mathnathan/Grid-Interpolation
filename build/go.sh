#!/bin/bash
make&&./read_hycom.x&&matlab -nodesktop -nosplash -r "run ../matlab/read_3d_field.m"
