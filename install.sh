#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-r', '--reinstall', action='store_true',default=False, help='store a boolean [default %(default)s]')
parser.add_argument('-c', '--clean_only', action='store_true',default=False, help='store a boolean [default %(default)s]')
parser.add_argument('-n', '--native', action='store_true',default=False, help='store a boolean [default %(default)s]')
parser.add_argument('-d', '--db_location', default="NA", type=str,help='Provide location the database to be downloaded [default %(default)s]')

EOF




DIR=`dirname $(readlink -f "$0")`

if [[ "$DB_LOCATION" != "NA" ]]
then
echo the answer: "$DB_LOCATION"
else
DB_LOCATION=$DIR
fi







# Set default options.
CLEAN_ONLY=false
FORCE=false
REINSTALL=false






declare -A DB_ID_HUMAN=(
['viral_vipr']='1fIxhnwNSPj6NL2R44bqYkYu2T8OLfqpk'
['fungi']='1yBeBjnrnHtxZliruu3oC8NjZ3wHg-WQ3'
['BWAindex']='19Uscw8KrPyUiuPcErrXpyOxN0PqUtbOZ'
['protozoa']='1_dPn8kk3I--Icy0gwTorFneV1sor1dU2'
)

declare -A DB_MD5_HUMAN=(
['viral_vipr']='9dce447328dfbc3a62cc7dd5b052242f'
['fungi']='9f2d304fd5c49981682b2bb7a900a30e'
['BWAindex']='4f009e3732d9f513e7b19b58edc41c13'
['protozoa']='23e12115a5e9d526553c901e772731f5'
)


ORGANISM='human'
download_list=$'BWAindex\nviral_vipr\nfungi\nprotozoa'







cd "$DB_LOCATION"
mkdir "db_$ORGANISM"
cd "db_$ORGANISM"





echo '----- Downloading human and mirobial references ----------------------------------------------------'
for download in $download_list; do
echo "Downloading item: $download for $ORGANISM"
success=false
while [ $success = false ]; do
case "$ORGANISM" in
human)
db_id="${DB_ID_HUMAN[$download]}"
db_md5="${DB_MD5_HUMAN[$download]}"
;;
mouse)
db_id="${DB_ID_MOUSE[$download]}"
db_md5="${DB_MD5_MOUSE[$download]}"
;;
*)
echo 'Error: Unknown ORGANISM.' >&2
exit 1
;;
esac
confirm_code=`curl --silent --insecure --cookie-jar cookies.txt "https://docs.google.com/uc?export=download&id=$db_id" | sed -rn 's .*confirm=([0-9A-Za-z_]+).* \1\n p'`
curl --location --insecure --cookie cookies.txt "https://docs.google.com/uc?export=download&confirm=$confirm_code&id=$db_id"  >$download.tar.gz





rm cookies.txt
if [ `md5sum "$download.tar.gz" | sed 's \(.*\)\ .* \1 '` = "$db_md5" ]; then
tar -zxvf "$download.tar.gz"
rm "$download.tar.gz"
success=true
else
echo "Download of $download for $ORGANISM failed (checksum" \
'mismatch. Retrying.'
fi
done
done






# ------------------------------------------------------------------------------
# DOWNLOAD TOOLS
# ------------------------------------------------------------------------------

cd "$DIR/tools"

pwd
ls
echo $REINSTALL

# Skip this section if neither -c nor -r are selected and there is a previous
# installation (as indicated by the presence of the imrep directory).
echo '----- Checking for existing installations --------------------------------------'
if [ $CLEAN_ONLY ] && [ $REINSTALL ] && [ -d 'imrep' ]; then
    echo 'Existing installation found. Skipping tools download. To reinstall,' \
        'please use the -r option.'
else
    echo '----- Removing previous versions -----------------------------------------------'
    rm -fr imrep metaphlan2 MiniConda
    if [ $CLEAN_ONLY ]; then
        echo 'Done: Cleaning complete.'
        exit 0
    fi



    #Download megahit
    echo '----- Downloading Megahit --------------------------------------------------'
    git clone https://github.com/voutcn/megahit.git
    cd megahit
    make
    cd ..



    # Download MiniConda and add shebangs.
    echo '----- Setting up Python environment --------------------------------------------'
    if [ $NATIVE ]; then
        ./install-MiniConda.sh
        cd MiniConda/lib
        ln -s libncursesw.so.5 libtinfow.so.5
        cd ../..
        MiniConda="$PWD/MiniConda/bin/python"
    #    sed -i "1c #!$MiniConda" metaphlan2/metaphlan2.py
    #    sed -i "1c #!$MiniConda" metaphlan2/strainphlan.py
    #    sed -i "1c #!$MiniConda" metaphlan2/utils/read_fastx.py
    #else
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/metaphlan2.py
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/strainphlan.py
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/utils/read_fastx.py
    fi
fi


cd "$DIR"
if [ `readlink -e "$DB_LOCATION"` != "$DIR" ]; then
ln -s "$DB_LOCATION/db_$ORGANISM"
fi
echo "Done: Reference databases are ready"







