#!/bin/bash
#
# Example:
# > ./papermovie.sh analysis/decals_gz_analysis.pdf

# Inspiration: https://github.com/brownsys/paper-movies

# Had to install ImageMagick with homebrew
# Fixed font list via http://stackoverflow.com/questions/24696433/why-font-list-is-empty-for-imagemagick

filename=$1

git pull origin pm2
git checkout pm2

git log -- $filename  > pdflog
grep "commit" pdflog | awk '{print $2}' > commits.txt
mv commits.txt movie/commits.txt
rm -Rf pdflog

rev=1

while read p; do
    echo "Commit $p"
    git reset --hard $p
    cp $filename movie/$rev.pdf

    echo $rev

    montage movie/$rev.pdf -tile 7x2 -background white -geometry 213x275-15+4 movie/$rev.png

    let rev+=1

done < movie/commits.txt

# Renumber, crop, and reverse order of the frames

x=0
for f in `ls -1 movie/*.png | grep -v "movie\-1.png" | sort -nr`; do
    n=`printf "%03d" $x`
    mv $f movie/$n.png
    convert movie/$n.png -crop 1280x568+0+0 movie/a$n.png
    mv -f movie/a$n.png movie/$n.png
    let x+=1
done

arr=$(echo $filename | tr ";" "\n")
IFS='/' read -a fnameonly <<< "$filename"
x=${fnameonly[${#fnameonly[@]}-1]}
IFS='.' read -a stub <<< "$x"

ffmpeg -r 5 -i movie/%03d.png -vcodec mpeg4 -pix_fmt yuv420p -b 8000k movie/$stub.mov

rm -f movie/*.pdf
rm -f movie/*.png
rm -f movie/commits.txt

# Checkout the most recent version again
git pull origin pm2
git checkout pm2

