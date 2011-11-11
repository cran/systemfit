#!/bin/sh
echo "Sweave(\"systemfit.Rnw\")" | LC_ALL="C" R --no-save --no-restore

echo "Stangle(\"systemfit.Rnw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | LC_ALL="C" R --no-save --no-restore
