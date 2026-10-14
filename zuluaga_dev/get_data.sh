# Download full observations 
if [ ! -f data/obs_sbn_X05_full.parquet ]; then
    # wget -P data/ https://storage.googleapis.com/asteroid-institute-public/production/rubin/mpc/obs_sbn/full/csv/obs_sbn_X05_full.csv.gz
    wget -P data/ https://storage.googleapis.com/asteroid-institute-public/production/rubin/mpc/obs_sbn/full/parquet/obs_sbn_X05_full.parquet
else
    echo "Data already exists: data/obs_sbn_X05_full.parquet"
fi 

# Download extended MPC observations 
if [ ! -f data/mpcorb_extended.dat.gz ]
then 
    wget -P data/ https://minorplanetcenter.net/Extended_Files/mpcorb_extended.dat.gz 
else
    echo "Data already exists: data/mpcorb_extended.dat.gz"
fi

# Download asteroids with rotation from Google Drive
if [ ! -f data/asteroids-with-rotation.csv ]
then 
    gdown 1afQHJYa6bAyIhbP1_dwBoyY0YYzMI7y5
else
    echo "Data already exists: data/asteroids-with-rotation.csv"
fi

# Download X05 simplified
if [ ! -f data/X05_simplified.csv ]
then 
    gdown 1bre3OmMF1xvr3puWUfhhUq_OwRXFltS5
else
    echo "Data already exists: data/X05_simplified.csv"
fi