# Combine all cleaned up raw data into a single file for combined pca

files=(CleanUp_LCMSMS_HILIC_Positive_rawHeight.txt CleanUp_LCMSMS_RP_Negative_rawHeight.txt CleanUp_LCMSMS_RP_Positive_rawHeight.txt)

for i in 0 1 2; do
    wc -l Inputs/${files[$i]}
    if [[ $i -eq 0 ]]; then
        cp Inputs/${files[$i]} Inputs/Combined_rawHeight.txt
    else
        tail -n +3 Inputs/${files[$i]} >> Inputs/Combined_rawHeight.txt
    fi
    wc -l Inputs/Combined_rawHeight.txt
done