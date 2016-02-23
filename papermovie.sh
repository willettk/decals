#!/bin/bash

# Inspiration: https://github.com/brownsys/paper-movies

# Going after a 1280x720 (aka 720p HD) resolution
# 1080p would be 1920x1080, so could have higher-quality video

# Issues: don't have montage or ffmpeg on zrak. Will try on kintana.

filename=analysis/decals_gz_analysis.pdf

git log -- $filename  > pdflog
grep "commit" pdflog | awk '{print $2}' > commits.txt
rm -Rf pdflog

rev=1

while read p; do
    echo "Commit $p"
    git reset --hard $p
    cp $filename movie/$rev.pdf

    echo $rev
    let rev+=1

    cd movie

    montage $rev.pdf -tile 7x2 -background white -geometry 213x275-15+4 $rev.png

done < commits.txt


# Checkout the most recent version again

git pull
git checkout master

# We need to do some renumbering and cropping to make ffmpeg happy

x=0
for f in `ls -1 *.png | grep -v "\-1.png" | sort -n`; do
    n=`printf "%03d" $x`
    mv $f $n.png
    convert $n.png -crop 1280x568+0+0 a$n.png # the formula above results in 1281x568 images
    mv -f a$n.png $n.png
    let x+=1
done

rm *-1.png

ffmpeg -r 5 -i %03d.png -vcodec libx264 -pix_fmt yuv420p -b 8000k movie.mov

cd ..
