#!/bin/bash

python analysis.py $1 0.89 1.11 -0.002 0.025 1.e20 2.e20 1.e19
ffmpeg -r 30 -i $1/figure/snap%06da_ei.png -vcodec libx264 -pix_fmt yuv420p -r 60 $1/figure/movie.mp4
