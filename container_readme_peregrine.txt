# get the folder to store the image
img_folder=/data/${USER}/singularity/single_cell_container/
# get the folder to store the simulated disk
sim_disk=/data/${USER}/singularity/single_cell_container//simulated_disk/

# the paths accessible to the container, separated by a comma
paths=,/data/${USER}/

# create the simulated disk location
mkdir -p ${sim_disk}

# move to the location to store the image
cd ${img_folder}
# download the image
#wget https://drive.google.com/uc?export=download&id=1zzvl8eBsYyt9MHyTggBpjc3jCjMA5dOk
# download the image with a nifty little tool
git clone https://github.com/circulosmeos/gdown.pl.git
./gdown.pl/gdown.pl https://drive.google.com/file/d/1Y1tyXK6ojNZW4sj5O8hkeNFFLM0SIacT/view?usp=sharing singlecell_container.simg

# add an alias to start the container
echo 'alias Rsing="singularity exec --bind '${sim_disk}':/home/'${USER}${paths}' '${img_folder}'singlecell_container.simg R"' >> /home/${USER}/.bashrc
echo 'alias Rsingscript="singularity exec --bind '${sim_disk}':/home/'${USER}${paths}' '${img_folder}'singlecell_container.simg Rscript"' >> /home/${USER}/.bashrc
echo 'alias sccontainer="singularity shell --bind '${sim_disk}':/home/'${USER}${paths}' '${img_folder}'singlecell_container.simg"' >> /home/${USER}/.bashrc

# in case you've never used singularity before
#mkdir -p ~/singularity # if you do have this directory already, move everything to /data/${USER}/singularity/
mkdir -p ~/.singularity
# move the directories
mv ~/singularity /data/${USER}/
mv ~/.singularity /data/${USER}/
# symlink the directories
ln -s /data/${USER}/.singularity ~/.singularity
ln -s /data/${USER}/singularity ~/singularity
