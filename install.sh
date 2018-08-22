#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-r', '--reinstall', action='store_true',default=False, help='store a boolean [default %(default)s]')
parser.add_argument('-c', '--clean_only', action='store_true',default=False, help='store a boolean [default %(default)s]')
parser.add_argument('-n', '--native', action='store_true',default=False, help='store a boolean [default %(default)s]')
parser.add_argument('-f', '--force', action='store_true',default=False, help='store a boolean [default %(default)s]')
parser.add_argument('-d', '--db_location', default="NA", type=str,help='Provide location the database to be downloaded [default %(default)s]')
parser.add_argument('-l', '--link', default="NA", type=str,help='Provide location the database to be downloaded [default %(default)s]')


EOF


DIR=`dirname $(readlink -f "$0")`

#------------------------------------- Tools  -------------------------------------
cd ${DIR}/tools/


if [ $CLEAN_ONLY ]
then
echo '----- Removing previous versions -----------------------------------------------'
rm -fr MiniConda megahit
echo 'Done: Cleaning complete.'
exit 0
fi


if [ $REINSTALL ]
then
echo '----- Removing previous versions -----------------------------------------------'
rm -fr  MiniConda megahit
fi


if [ -d 'megahit' ]
then
echo 'Existing installation found. Skipping tools download. To reinstall, please use the -r option.'
else


#Download megahit
echo '----- Downloading Megahit --------------------------------------------------'
git clone https://github.com/voutcn/megahit.git
cd megahit
make
cd ..



# Download MiniConda and add shebangs.
echo '----- Setting up Python environment --------------------------------------------'
if [ $NATIVE ]
then
echo "--native option was choosen. And Minoconda instalation is skipped"
else
./install-MiniConda.sh
cd MiniConda/lib
ln -s libncursesw.so.5 libtinfow.so.5
cd ../bin

./conda install -c bioconda blast



MiniConda="$PWD/MiniConda/bin/python"
fi
fi



#------------------------------------- Databases  -------------------------------------






if [[ "$DB_LOCATION" != "NA" ]]
then
echo the answer: "$DB_LOCATION"
else
DB_LOCATION=$DIR
fi

echo "Database location ", $DB_LOCATION


cd $DIR
ORGANISM='human'


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


# if
if [ $FORCE ]
then
echo "FORSE option was choosen"
fi




download_list=$'BWAindex\nviral_vipr\nfungi\nprotozoa'

echo '----- Downloading human and mirobial references ----------------------------------------------------'


echo '----- Checking for existing databases ------------------------------------------'
if [ -h "db_$ORGANISM" ] || [ -d "db_$ORGANISM" ]; then
if [ $FORCE ]
then
echo 'Unlinking existing database.'
if [ -h "db_$ORGANISM" ]; then
rm "db_$ORGANISM"
else
rm -r "db_$ORGANISM"
fi
else
echo 'Existing database found. Skipping database download. To unlink the current database, please use the -f option.'
exit 0
fi
fi



if [ "$LINK" != 'NA' ]; then
echo '----- Linking database -----------------------------------------------------'
if [ -d "$LINK" ]; then
ln -s "$LINK"
echo 'Done: Database linked.'
exit 0
else
echo "Error: Link target doesn't exist." >&2
exit 1
fi
fi





cd "$DB_LOCATION"
mkdir "db_$ORGANISM"
cd "db_$ORGANISM"

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



cd "$DIR"
if [ `readlink -e "$DB_LOCATION"` != "$DIR" ]; then
ln -s "$DB_LOCATION/db_$ORGANISM"
fi
echo "Done: Reference databases are ready"




