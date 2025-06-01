while IFS= read -r filename; do
    rm -v "electron_ideal_sample/$filename"
done < output.txt

