#!sh

DIR=$1
PROJECT=$2

if [ ! -d $DIR ]
then
    echo "$DIR does not exist, exiting"
    exit
fi

if [ -d $DIR/$PROJECT ]
then 
    echo "error, $DIR/$PROJECT already exists, exiting"
fi

mkdir $DIR/$PROJECT

if [ ! -d $DIR/$PROJECT ]
then
    echo "$DIR/$PROJECT could not be made, exiting"
    exit
fi

cp -r skeleton/* $DIR/$PROJECT

cd $DIR/$PROJECT