DIR=`dirname ${0}`
PPWD=$PWD/$DIR
#echo $PPWD
export PATH="$PATH:$PPWD/DALIGNER:$PPWD/DAZZ_DB:$PPWD/DEXTRACTOR/:$PPWD/blasr/"
export PATH="$PATH:$PPWD/scripts"
export PATH="$PATH:$PPWD/src/build"