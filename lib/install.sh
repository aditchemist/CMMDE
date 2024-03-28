# Install the requirements
sudo apt install python3-pip
pip install ase pymatgen gmxscript prettytable acpype MDAnalysis statsmodels
## Aktifkan dua baris di bawah jika anda ingin menginstall cmmde-gui.  
#pip install -Iv panel-chemistry==0.0.11
#pip install -Iv jinja2==3.0.1
#pip install py3Dmol
sudo apt install -y cmake 
sudo apt install -y libblas-dev
sudo apt install -y libarpack2-dev
sudo apt install -y openbabel  
# Make the bashrc
## Merging the Python libraries of CMMDE
echo "export PYTHONPATH=`pwd`/lib:$PYTHONPATH" >> ~/.bashrc
## Merging the CMMDE binaries to your PATH
echo "export PATH=`pwd`/bin:$PATH" >> ~/.bashrc
echo "export ORCA_COMMAND=/opt/orca/orca" >> ~/.bashrc
echo "export XTB_COMMAND=/usr/local/bin/xtb" >> ~/.bashrc
## PATH to DCDFTBMD program, please refer the README.md to register the program at Nakai Lab of Waseda University
## Aktifkan baris di bawah jika anda sudah menginstall DCDFTBMD
#echo "export DCDFTB_COMMAND=/path/to/DCDFTBMD/bin/binary" >> ~/.bashrc
# If you have GROMACS
#echo "export GROMACS_COMMAND=`pwd`/lib/cmmde_solution.py" >> ~/.bashrc
#echo "export GROMACS_DIR=/path/to/gromacs/build/bin" >> ~/.bashrc
#echo "alias gmx='path/to/gromacs/bin/gmx'" >> ~/.bashrc
#echo "export GROMACS_FF=/path/to/gromacs/source/share/top"
#echo "export AMBER_FF=/path/to/miniconda/dat/leap/cmd"
#echo "export DFTB_COMMAND=/path/to/dftb+/bin/dftb+"
#echo "export DPTOOLS_DIR=/path/to/dftbplus/tools/dptools/bin"
# If you have DOCK6
#echo "export DMS_COMMAND=/usr/local/bin/dms" >> ~/.bashrc
#echo "export DOCK_DIR=/path/to/dock6/bin" >> ~/.bashrc
#echo "export DOCK_PARM=/path/to/dock6/parameters" >> ~/.bashrc

## Pada kasus-kasus tertentu, source ~/.bashrc akan perlu dilakukan secara manual jika cmmde.py tidak terdeteksi di $PATH
source ~/.bashrc
