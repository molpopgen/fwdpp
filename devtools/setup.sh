#!/usr/bin/env bash

usage(){
    >&2 echo "Usage: $0 options"
    >&2 echo "Mandatory options are:"
    >&2 echo " -d/--dir root path to place the project (DIR)"
    >&2 echo " -p/--project Name for project directory (PROJECT)"
    >&2 echo "Together, the -p and -d options will create a skeleton package in DIR/PROJECT."
    >&2 echo "For example: ./setup.sh -d ~/src -p my_fwdpp_project"
    >&2 echo "Optional options are:"
    >&2 echo " -u/--url project url.  This must be single-quoted with special characters escaped, e.g. -u 'https:\/\/github.com\/molpopgen\/fwdpp'"
    exit 1
}
while true; do
    case "$1" in
	-d | --dir ) DIR="$2"; shift 2;;
	-p | --project ) PROJECT="$2"; shift 2;;
	-u | --url ) URL="$2"; shift 2;;
        -- ) shift; break ;;
    * ) break ;;
  esac
done

##VALIDATE THE INPUT PARAMS
if [ -z ${DIR+x} ]; then >&2 echo "Error: no base path name specified"; usage; else echo "Base path set to '$DIR'"; fi
if [ -z ${PROJECT+x} ]; then >&2 echo "Error: project name specified"; usage; else echo "Project name set to '$PROJECT'"; fi

if [ ! -d $DIR ]
then
    echo "$DIR does not exist, exiting"
    exit
fi

if [ -d $DIR/$PROJECT ]
then 
    echo "error, $DIR/$PROJECT already exists, exiting"
    exit
fi

mkdir $DIR/$PROJECT

if [ ! -d $DIR/$PROJECT ]
then
    echo "$DIR/$PROJECT could not be made, exiting"
    exit
fi

cp -r skeleton/* $DIR/$PROJECT

mv $DIR/$PROJECT/src/FWDPPPACKAGE.cc $DIR/$PROJECT/src/$PROJECT.cc
sed -i "s/FWDPPPACKAGE/$PROJECT/" $DIR/$PROJECT/configure.ac
sed -i "s/FWDPPPACKAGE/$PROJECT/g" $DIR/$PROJECT/src/Makefile.am
sed -i "s/FWDPPPROJECTURL/$URL/" $DIR/$PROJECT/configure.ac
