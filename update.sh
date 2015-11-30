#!/bin/bash

if [ "$2" == "push" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' src/ $1@shannon.stanford.edu:/data/pacbio_assembly/AwesomeAssembler/src
fi

if [ "$2" == "pull" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' $1@shannon.stanford.edu:/data/pacbio_assembly/AwesomeAssembler/src/ src
fi

if [ "$2" == "update" ];
then ssh -t $1@shannon.stanford.edu "cd /data/pacbio_assembly/AwesomeAssembler && ./build.sh"
fi
