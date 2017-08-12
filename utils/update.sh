#!/bin/bash

if [ "$2" == "scripts" ];
then rsync -rizP --delete --exclude '.*' --exclude '*.pyc' --exclude 'figures' scripts/ $1@shannon.stanford.edu:/home/$1/AwesomeAssembler/scripts
fi

if [ "$2" == "utils" ];
then rsync -rizP --delete --exclude '.*' --exclude '*.pyc' --exclude 'figures' utils/ $1@shannon.stanford.edu:/home/$1/AwesomeAssembler/utils
fi

if [ "$2" == "push" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' src/ $1@shannon.stanford.edu:/home/$1/AwesomeAssembler/src
fi

if [ "$2" == "pull" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' $1@shannon.stanford.edu:/home/$1/AwesomeAssembler/src/ src
fi

if [ "$2" == "update" ];
then ssh -t $1@shannon.stanford.edu "export TEMP=/home/$1/tmp && cd /home/$1/AwesomeAssembler && ./utils/build.sh"
fi

if [ "$2" == "all" ];
then rsync -rizP --delete --exclude '.*' --exclude 'data' --exclude '*.pyc' --exclude 'figures' --exclude 'build' . $1@shannon.stanford.edu:/home/$1/AwesomeAssembler
fi
