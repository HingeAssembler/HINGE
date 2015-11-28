#!/bin/bash

if [ "$2" == "push" ];
then rsync -ri --delete --exclude '.*' --exclude 'build' src/ $1@shannon.stanford.edu:/data/pacbio_assembly/AwesomeAssembler/src
fi

if [ "$2" == "pull" ];
then rsync -ri --delete --exclude '.*' --exclude 'build' $1@shannon.stanford.edu:/data/pacbio_assembly/AwesomeAssembler/src/ src
fi