# Temporarily change directory to $HOME to install software
pushd .
cd $HOME
# Make sure some level of pip is installed
python -m ensurepip
sudo apt-get update

# Install Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda

# Configure miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda

# Restore original directory
popd
