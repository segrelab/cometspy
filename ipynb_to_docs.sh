#!/bin/bash

for f in *.ipynb
do
    jupyter-nbconvert --to markdown $f --output $HOME/Dropbox/projects/comets-manual/docs/python-module/"${f%%.*}".md
done

