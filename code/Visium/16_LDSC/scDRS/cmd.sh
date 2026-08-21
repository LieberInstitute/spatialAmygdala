
DEST=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS
mkdir -p "$DEST/logs"

# copy only; never overwrite an existing file, never delete anything
for f in ./*.sh ./*.py ./*.tsv ./*.md; do
    b=$(basename "$f")
    if [ -e "$DEST/$b" ]; then
        echo "EXISTS, left alone : $b"
    else
        cp "$f" "$DEST/$b"
        echo "created            : $b"
    fi
done

chmod u+x "$DEST"/*.sh
echo "---- $DEST ----"
ls -la "$DEST"
