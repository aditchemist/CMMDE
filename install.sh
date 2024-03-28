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
# Make the zshrc
## Merging the Python libraries of CMMDE
echo "export PYTHONPATH=`pwd`/lib:$PYTHONPATH" >> ~/.zshrc
## Merging the CMMDE binaries to your PATH
echo "export PATH=`pwd`/bin:$PATH" >> ~/.zshrc
echo "export ORCA_COMMAND=/opt/orca/orca" >> ~/.zshrc
echo "export XTB_COMMAND=/usr/local/bin/xtb" >> ~/.zshrc
## PATH to DCDFTBMD program, please refer the README.md to register the program at Nakai Lab of Waseda University
## Aktifkan baris di bawah jika anda sudah menginstall DCDFTBMD
#echo "export DCDFTB_COMMAND=/path/to/DCDFTBMD/bin/binary" >> ~/.zshrc
# If you have GROMACS
#echo "export GROMACS_COMMAND=`pwd`/lib/cmmde_solution.py" >> ~/.zshrc
#echo "export GROMACS_DIR=/path/to/gromacs/build/bin" >> ~/.zshrc
#echo "alias gmx='path/to/gromacs/bin/gmx'" >> ~/.zshrc
#echo "export GROMACS_FF=/path/to/gromacs/source/share/top"
#echo "export AMBER_FF=/path/to/miniconda/dat/leap/cmd"
#echo "export DFTB_COMMAND=/path/to/dftb+/bin/dftb+"
#echo "export DPTOOLS_DIR=/path/to/dftbplus/tools/dptools/bin"
# If you have DOCK6
#echo "export DMS_COMMAND=/usr/local/bin/dms" >> ~/.zshrc
#echo "export DOCK_DIR=/path/to/dock6/bin" >> ~/.zshrc
#echo "export DOCK_PARM=/path/to/dock6/parameters" >> ~/.zshrc

## Pada kasus-kasus tertentu, source ~/.zshrc akan perlu dilakukan secara manual jika cmmde.py tidak terdeteksi di $PATH
source ~/.zshrc
