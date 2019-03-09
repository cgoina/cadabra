#!/bin/sh

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

echo "$SCRIPTPATH"

sudo ln -s ${SCRIPTPATH}/bin/libavbin64.dylib /usr/local/lib/libavbin64.dylib
