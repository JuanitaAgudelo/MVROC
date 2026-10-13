# Download observations 
if [ ! -f data/obs_sbn_X05_full.parquet ]; then
    # wget -P data/ https://storage.googleapis.com/asteroid-institute-public/production/rubin/mpc/obs_sbn/full/csv/obs_sbn_X05_full.csv.gz
    wget -P data/ https://storage.googleapis.com/asteroid-institute-public/production/rubin/mpc/obs_sbn/full/parquet/obs_sbn_X05_full.parquet
else
    echo "Data already exists: data/obs_sbn_X05_full.parquet"
fi 

if [ ! -f ]
wget -P https://minorplanetcenter.net/Extended_Files/mpcorb_extended.dat.gz 