#DIR=`dirname ${0}`
PPWD=$PWD
#echo $PWD
#echo $DIR
export PATH="$PATH:$PPWD/thirdparty/DALIGNER:$PPWD/thirdparty/DAZZ_DB:$PPWD/thirdparty/DEXTRACTOR/:$PPWD/thirdparty/DASCRUBBER"
export PATH="$PATH:$PPWD/scripts"
export PATH="$PATH:$PPWD/build/bin/consensus:$PPWD/build/bin/filter:$PPWD/build/bin/layout"