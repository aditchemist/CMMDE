#!/usr/bin/env python3
import sys
import os
import argparse
from cmmde_orca import orca
from cmmde_dcdftb import dcdftb
from cmmde_dock import readpdb, splitpdb, addH, addcharge, sphgen, showsphere, gridgen, rigiddock,flexdock, translig, sdf2xyz, multiopt, checkopt, multiflexdock
from cmmde_mdpro import proprep, ligprep
from cmmde_dftb import poscar2gen, vasp2gen, dftb
from cmmde_qe import qe
from cmmde_nw import nwchem
from cmmde_xtb import xtb
import time


parser = argparse.ArgumentParser(description='Computational Molecular & Material Design Interface')
# Input Geometri, tanpa ini perhitungan tidak akan berjalan.
parser.add_argument('-i','--input',type=str,default='geom.xyz',help='Input geometri dalam koordinat Cartesian')
parser.add_argument('-cons','--constraints',type=str,default='None',help='Membekukan ikatan, sudut ikatan, dan sebagainya selama proses optimasi geometri. Periksa manual Orca 5.0 lebih lanjut mengenai ini.')
# Informasi muatan dan multiplisitas spin molekul.
parser.add_argument('-c','--charge',type=int,default=0, help='Muatan molekul')
parser.add_argument('-mult','--mult',type=int,default=1, help='Multiplisitas molekul')
# Input yang berkaitan dengan software, metode, dan jenis pekerjaan
parser.add_argument('-s','--software',type=str,default='orca',help='Software yang digunakan. Pilihan: orca, dcdftb, gromacs')
parser.add_argument('-j','--job',type=str,default='sp',help='Jenis perhitungan yang dilakukan. Pilihan: sp, opt, freq, md, mtd, td, dock.')
parser.add_argument('-m','--method',type=str,default='XTB2',help='Metode yang digunakan dalam perhitungan. Pilihan: XTB, XTB2, DFTB, DFTB2, DFTB3-diag, dan sebagainya.')
# Input yang berkaitan dengan perhitungan frekuensi menggunakan software Orca
parser.add_argument('-sfreq','--scalefreq',type=float,default=1, help='Faktor skala frekuensi yang digunakan.')
# Input yang berkaitan dengan perhitungan menggunakan software DCDFTBMD
parser.add_argument('-disp','--dispersion',type=str,default='None', help='Model koreksi dispersi jika menggunakan software DCDFTB. Pilihan: None, D3, D3BJ, D3H5.')
parser.add_argument('-para','--parapath',type=str,default='/Users/adit/opt/dftbplus/external/slakos/origin/3ob-3-1',help='Lokasi folder berisikan himpunan parameter DFTB yang akan digunakan.')
parser.add_argument('-iter','--iter',type=int,default=9999, help='Jumlah iterasi dalam optimasi geometri dan jenis perhitungan lainnya')
parser.add_argument('-ens','--ensembel',type=str,default='NVE', help='Ensembel yang digunakan dalam simulasi dinamika molekul. Pilihan: NVE, NVT, dan NPT.')
parser.add_argument('-tstat','--thermostat',type=str,default='andersen', help='Termostat yang digunakan dalam simulasi NVT. Pilihan: andersen,berendsen,dan nose')
parser.add_argument('-t','--temp', type=str, default='298.15',help='Suhu yang digunakan dalam perhitungan frekuensi maupun simulasi MD dan perhitungan frekuensi dalam satuan Kelvin.')
parser.add_argument('-dt','--deltat',type=float, default=1.0, help='Selang waktu integrasi dalam simulasi dinamika molekul dalam satuan femtodetik.')
parser.add_argument('-mdprint','--mdprint',type=int,default=10, help='Berapa step sekali struktur dicetak. Default=10.')
parser.add_argument('-rest','--restart',type=str,default='false',help='Apakah dilakukan restart simulasi MD?')
parser.add_argument('-nvtdir','--nvtdir',type=str,default='../NVT')
parser.add_argument('-ns','--nstep', type=int, default=50000,help='Banyaknya step dalam simulasi MD.')
parser.add_argument('-np','--nproc',type=int,default=1, help='Jumlah CPU yang digunakan')
# Opsi kondisi PBC
parser.add_argument('-a1','--a1',type=float, default=0, help='panjang vektor a1')

parser.add_argument('-a2','--a2',type=float, default=0, help='panjang vektor a2')
parser.add_argument('-a3','--a3',type=float, default=0, help='panjang vektor a3')
parser.add_argument('-b1','--b1',type=float, default=0, help='panjang vektor b1')
parser.add_argument('-b2','--b2',type=float, default=0, help='panjang vektor b2')
parser.add_argument('-b3','--b3',type=float, default=0, help='panjang vektor b3')
parser.add_argument('-c1','--c1',type=float, default=0, help='panjang vektor c1')
parser.add_argument('-c2','--c2',type=float, default=0, help='panjang vektor c2')
parser.add_argument('-c3','--c3',type=float, default=0, help='panjang vektor c3')
# Input spesifik program DFTB+
parser.add_argument('-hcorr','--hcorr',type=str, default='hdamp',help='Koreksi ikatan hidrogen. Default: hdamp. Pilihan: H5.')
# Input yang berkaitan dengan penyiapan sistem larutan
parser.add_argument('-mt', '--terlarut', type=str, help='Nama file molekul terlarut dalam format .xyz (ditulis tanpa ekstensi).')
parser.add_argument('-ct', '--c_terlarut', type=str, help='Muatan bersih molekul terlarut.',default=0)
parser.add_argument('-pt', '--persen_terlarut', type=str, default=50, help='Persen massa terlarut. Jika lebih dari satu terlarut, pisahkan dengan koma.')
parser.add_argument('-pp', '--persen_pelarut', type=str, default=50, help='Persen massa pelarut. Default:50')
parser.add_argument('-l', '--lapang', type=float, default=5, help='Panjang tambahan (Angstrom) pada rusuk kubus untuk menghindari bad contact. Default: 10.0.')
parser.add_argument('-rc','--rcol',type=float,default=0, help='Additional Coulombic radius for PME calculation')
parser.add_argument('-gen','--generate_dftbinp',type=str,default='false',help='Apakah ingin mengkonversi ke dalam format koordinat xyz?')
parser.add_argument('-mp', '--pelarut', type=str, help='Nama file molekul pelarut dalam format .xyz (ditulis tanpa ekstensi).')
parser.add_argument('-nequil','--nequil',type=int, default=50000, help='Jumlah step pada saat ekuilibrasi.')
parser.add_argument('-nnpt','--nnpt',type=int, default=50000, help='Jumlah step pada saat ekuilibrasi NPT.')
parser.add_argument('-comp','--compress',type=float, default=4.5e-6, help='Nilai kompresibilitas. Default=4.5e-6.')
parser.add_argument('-nprod','--nprod',type=int, default=400000, help='Jumlah step pada saat production.')
parser.add_argument('-prod', '--production', type=str, default='None', help='Ensembel yang digunakan dalam production run. Pilihan:nve,npt,nvt.')
parser.add_argument('-cp', '--c_pelarut', type=int, default=0, help='Muatan bersih molekul pelarut')
parser.add_argument('-ctype','--charge_type',type=str, default='gas', help='Tipe muatan yang digunakan untuk parameterisasi muatan.')
parser.add_argument('-Nump','--NumPelarut', type=int, default=100, help='Jumlah molekul pelarut maksimum dalam sistem larutan. Default = 100.')
parser.add_argument('-cat','--cation',type=str, default='none',help='Kation yang digunakan untuk menetralkan muatan sistem larutan.')
parser.add_argument('-p', '--pressure', default=1.0, type=float, help='Tekanan dalam satuan bar')
parser.add_argument('--packmol', type=str, default='/Users/adit/opt/packmol/packmol')
# Opsi untuk melakukan restart simulasi MD
parser.add_argument('-traj','--traject',type=str, default='../NVT/traject',help='File trayektori dari simulasi sebelumnya')
parser.add_argument('-vel','--velocity',type=str, default='../NVT/velocity',help='File kecepatan atom dari simulasi sebelumnya')
parser.add_argument('-inp','--dftbinp',type=str, default='../NVT/dftb.inp',help='File dftb.inp dari simulasi sebelumnya.')

# Opsi untuk dinding potensial virtual
parser.add_argument('-soft','--soft',type=str,default='false',help='Apakah dinding potensial akan digunakan?. Pilihan: true, false.')
parser.add_argument('-softtype','--softtype',type=str,default='SPHERE',help='Jenis dinding apakah yang akan anda pakai? Baca lebih lanjut manual DCDFTBMD.')
parser.add_argument('-softrange','--softrange',type=float,default=10,help='Ukuran dinding potensial yang anda gunakan. Default = 10 angstrom')
parser.add_argument('-softcenter','--softcenter',type=str,default='COM',help='Jenis pusat koordinat dinding potensial. Default = COM')
parser.add_argument('-metarest','--metarest',type=str,default='false',help='Apakah dilakukan restart simulasi metadinamika?. Pilihan: true, false.')
parser.add_argument('-metafreq','--metafreq',type=int,default=100, help='Dalam berapa step sekali potensial Gaussian ditambahkan?')
parser.add_argument('-metaheight','--metaheight',type=float,default=3.0e-3,help='Ketinggian potensial Gaussian dalam satuan Hartree. Default = 3.0e-3')
parser.add_argument('-cv','--cvtype',type=str,default='coordnum',help='Pilihan collective variable (CV) yang digunakan. Pilihan: coordnum, distance, angle, dihedral, distancediff, distanceadd, meandistance, pointplanedistance.')
parser.add_argument('-metawidth','--metawidth',type=float,default=0.1,help='Lebar potensial Gaussian yang digunakan (dalam satuan yang sama dengan satuan cv). Default = 0.1')
parser.add_argument('-pow1','--pow1',type=float,default=6,help='Nilai pangkat pertama pada definisi bilangan koordinasi rasional. Default = 6')
parser.add_argument('-pow2','--pow2',type=float,default=12,help='Nilai pangkat kedua pada definisi bilangan koordinasi rasional. Default = 12')
parser.add_argument('-rcut','--rcut',type=float,default=1.6,help='Jarak cutoff pada definisi bilangan koordinasi dalam satuan angstrom. Default = 1.6')
parser.add_argument('-fesstart','--fesstart',type=float,default=0, help='Titik minimum CV. Default = 0.')
parser.add_argument('-fesend','--fesend',type=float,default=1, help='Titik maksimum CV. Defalt = 1.')
parser.add_argument('-fesbin','--fesbin',type=float,default=0.01,help='Selang CV. Default = 0.01')
parser.add_argument('-ag1','--ag1',type=str,help='Grup atom pertama dalam sebuah CV.')
parser.add_argument('-ag2','--ag2',type=str,help='Grup atom kedua dalam sebuah CV.')
parser.add_argument('-ag3','--ag3',type=str,help='Grup atom ketiga dalam sebuah CV.')
parser.add_argument('-ag4','--ag4',type=str,help='Grup atom keempat dalam sebuah CV.')
# Input terkait DC setup
parser.add_argument('-bufrad','--bufrad',type=float,default=3.0,help='Radius buffer yang digunakan dalam sebuah subsistem. Default: 3.0')
parser.add_argument('-delta','--delta',type=float,default=3.0,help='Panjang rusuk kubus virual yang dibuat untuk membelah-belah subsistem. Default: 3.0')
parser.add_argument('-opttype','--opttype',type=str,default='bfgs', help='Jenis algoritma optimasi geometri yang digunakan dalam program DFTBUP maupun DCDFTBMD. Default: bfgs. Pilihan: sd, cg, qm, fire.')
parser.add_argument('-freqtype','--freqtype',type=int,default=1, help='Jenis perhitungan frekuensi vibrasi yang dilakukan, analitik (1) atau numerik (2). Default: 1, pilihan: 2.')
parser.add_argument('-econv','--econv',type=float,default=1e-6, help='Batas konvergensi perhitungan energi.')
parser.add_argument('-dconv','--dconv',type=float,default=1e-6, help='Batas konvergensi perhitungan gradien.')
# Input yang berkaitan dengan persiapan docking
parser.add_argument('-ligname','--ligname',type=str,help='Nama ligan yang ingin diekstrak strukturnya')
parser.add_argument('-ligand','--ligand',type=str,help='File ligand dalam format .mol2. Format harus .xyz jika ligand bukan merupakan native ligand.')
parser.add_argument('-dockrange','--dockrange',type=float,default=10,help='Radius yang diperhitungkan sebagai sisi aktif di sekitar ligan asli.')
parser.add_argument('-chargetype','--chargetype',type=str,default='qtpie',help='Tipe muatan parsial untuk ligan. Pilihan: qtpie, gasteiger, eem, eem2015ba, eem2015bm, eem2015bn, eem2015ha, eem2015hm, eem2015hn, eqeq, mmff94, dan qeq.')
parser.add_argument('-calcrmsd','--calcrmsd',type=str,default='no',help='Apakah ingin menghitung RMSD saat docking?')
parser.add_argument('-nlig','--nligands',type=int,help='Jumlah ligan yang akan didocking.')
# Model Pelarut Implisit
parser.add_argument('-solvent','--solvent',type=str,default='none',help='Pelarut yang digunakan dalam perhitungan.')
# Input untuk simulasi MD sistem protein
parser.add_argument('-protein','--protein',type=str,help='Struktur protein dalam format pdb')
# Input untuk perkiraan jalur reaksi menggunakan xTB
parser.add_argument('-nrun','--nrun',type=int,default=1,help='Jumlah banyaknya sampling path yang dilakukan. Default: 1.')
parser.add_argument('-npoint','--npoint',type=int,default=25,help='Jumlah banyaknya titik untuk interpolasi jalur reaksi. Default: 25.')
parser.add_argument('-anopt','--anopt',type=int,default=10,help='Jumlah maksimum step optimasi geometri yang dilakukan. Default: 10.')
parser.add_argument('-kpush','--kpush',type=float,default=0.003,help='Faktor skala untuk variasi nilai RMSD. Default: 0.003 atomic unit.')
parser.add_argument('-kpull','--kpull',type=float,default=-0.015,help='Tetapan pegas untuk menarik posisi atom. Default: -0.015 atomic unit.')
parser.add_argument('-ppull','--ppull',type=float,default=0.05,help='Kekuatan tarikan pegas teroptimasi. Default: 0.05 atomic unit.')
parser.add_argument('-alp','--alp',type=float,default=1.2,help='Lebar potensial Gaussian dalam satuan Angstrom. Default: 1.2 Angstrom.')
# Scanning geometri menggunkaan xTB standalone
parser.add_argument('-dist','--distance',type=str,default='None',help='Jarak antara dua buah atom dalam format: Serial1,Serial2,Jarak. Contoh: 4,1,1.5.')
parser.add_argument('-ang','--angle',type=str,default='None',help='Sudut antara tiga buah atom dalam format: Serial1,Serial2,Serial3,Sudut. Contoh: 4,1,3,180.')
parser.add_argument('-dih','--dihedral',type=str, default='None',help='Sudut dihedral antara empat buah atom dalam format: Serial1,Serial2,Serial4,Dihedral. Contoh: 4,1,3,2,120.')
parser.add_argument('-scanmode','--scanmode',type=str,default='None',help='Mode scan geometri yang diinginkan. Pilihan: concerted, sequential.')
parser.add_argument('-scan','--scan',type=str,default='None',help='Rentang Scanning dalam format: Titik1,Titik2,JumlahTitik. Contoh:2.5,1.0,100.')
# Konstrain dan Fixing Atom pada Program XTB Standalone
parser.add_argument('-fixatm','--fixedatoms',type=str,default='None',help='Label atom-atom yang dibuat fix.')
parser.add_argument('-fixele','--fixedelements',type=str,default='None',help='Simbol unsur atom yang dibuat fix.')
# Input untuk perhitungan NEB
parser.add_argument('-produk','--produk',type=str,help='Struktur produk yang digunakan dalam interpolasi NEB dalam format koordinat Cartesian.')
parser.add_argument('-trans','--transitionstate',default='None',type=str,help='Struktur perkiraan keadaan transisi yang digunakan dalam format koordinat Cartesian.')
# Input untuk IRC
parser.add_argument('-irciter','--irciter',type=int,default=20,help='Iterasi maksimum IRC Orca.')
parser.add_argument('-printlevel','--printlevel',type=str,default='1',help='Print level dalam IRC Orca. Default=1')
parser.add_argument('-inithess','--inithess',type=str,default='Read',help='Cara menginisiasi Hessian. Default: Read.')
parser.add_argument('-grid','--grid',type=int,default=2,help='Ukuran grid dalam integrasi Orca.')
parser.add_argument('-finalgrid','--finalgrid',type=int,default=4,help='Ukuran final grid dalam integrasi Orca.')


# Input untuk perhitungan bergantung waktu
parser.add_argument('-nr','--nroots',type=int,default=5,help='Jumlah orbital aktif yang diperhitungkan dalam perhitungan bergantung waktu.')
parser.add_argument('-tda','--tda',type=str,default='false',help='Apakah akan digunakan pendekatan Tamm-Dancoff. Default = false.')
# Input untuk perhitungan TD-DFTB
parser.add_argument('-ts','--targetstate',type=int,default=0, help='Target orbital transisi yang dituju.')
parser.add_argument('-multtrans','--multtrans',type=int,default=1, help='Multiplisitas spin transisi elektronik.')
parser.add_argument('-ocsstr','--ocsstr',type=str,default='true', help='Menghitung osscilator strength dan momen dipol transisi elektronik.')
parser.add_argument('-wt','--writetrans',type=str,default='true', help='Menulis informasi detail mengenai transisi elektronik yang terjadi.')
parser.add_argument('-lc','--longrange',type=str,default='false', help='Apakah akan dilakukan koreksi interaksi jarak jauh.')
# Input untuk metode multiskala
parser.add_argument('-qmatoms','--qmatoms',type=str,default='None',help='Indeks dari atom-atom di lapisan QM. Indeks dimulai dari nol.')
parser.add_argument('-activeatoms','--activeatoms',type=str,default='1:-1',help='Indeks atom aktif yang digunakan untuk perhitungan.')
parser.add_argument('-hessfile','--hessfile',type=str,default='None',help='File Hessian hasil optimasi geometri atau berbagai perhitungan sebelumnya.')
parser.add_argument('-tc','--totalcharge',type=int,default=0,help='Muatan total molekul.')
parser.add_argument('-tmult','--totalmult',type=int,default=1,help='Multiplisitas spin total molekul.')
parser.add_argument('-qm2method','--qm2method',type=str,default='None',help='Custom method untuk lapisan QM2')
parser.add_argument('-qm2basis','--qm2basis',type=str,default='None',help='Custom himpunan basis untuk lapisan QM2')
# Informasi input padatan
parser.add_argument('-kpts','--kpts',type=str,default='1x1x1',help='Informasi K-points yang digunakan. Ditulis sebagai: 1x1x1, 2x2x2, dsb.')
# Input terkait perhitungan Atom in Molecules
parser.add_argument('-aim','--aim',type=str, default='false',help='Apakah akan dilakukan perhitungan AIM untuk interaksi non-kovalen?')
# Input terkait penggunaan Quantum-Espresso
parser.add_argument('-mode','--mode',type=str,default='from_scratch',help='Mode perhiungan dalam software Quantum-Espresso. Default: from_scratch. Pilihan: from_scratch, restart. ')
parser.add_argument('-pseudo','--pseudo',type=str,default='/home/adit/opt/qe-7.0/pseudo',help='Folder tempat menyimpan file pseudo potential.')
parser.add_argument('-outdir','--outdir',type=str,default='./out',help='Folder tempat menyimpan output.')
parser.add_argument('-unit','--unit',type=str,default='angstrom',help='Satuan yang digunakan dalam mendefinisikan parameter sel dan posisi atom. Default: angstrom. Pilihan: angstrom, bohr')
parser.add_argument('-bravais','--bravais',type=int,default=0,help='Indeks Bravais yang diperlukan untuk mengidentifikasi kristal pada software Quantum-Espresso.')
parser.add_argument('-ecutwfc','--ecutwfc',type=float,default=60.0,help='Cutoff energi kinetik untuk fungsi gelombang dalam satuan Ry. Default: 60.0 Ry.')
parser.add_argument('-ecutrho','--ecutrho',type=float,default=720.0,help='Cutoff energi kinetikr untuk rapat muatan dan potensial. Default: 720.0 Ry.')
parser.add_argument('-mix','--mixing_beta',type=float,default=0.7,help='Koefisien mixing untuk self-consistency.')
parser.add_argument('-conv_thr','--conv_thr',type=float,default=1e-8,help='Batas nilai konvergensi.')
parser.add_argument('-dftfunc','--dftfunc',type=str,default='pbe',help='Fungsional DFT yang digunakan. Default: pbe')
parser.add_argument('-extpseudo','--extpseudo',type=str,default='UPF',help='Ekstensi dari file pseudopotensial. Default: UPF.')
parser.add_argument('-optalgo','--optalgo',type=str,default='bfgs',help='Algoritma untuk optimasi geometri. Default: bfgs. Pilihan: damp, fire, verlet, langevin, langevin-smc, beeman.')
parser.add_argument('-cpress','--cellpress',type=float,default=0.0,help='Tekanan sel dalam satuan kilobar yang digunakan dalam Quantum-Espresso. Default:0.0.')
parser.add_argument('-press_conv','--press_conv_thr',type=float,default=0.5,help='Kriteria konvergensi untuk tekanan sel. Default=0.5 kbar.')
parser.add_argument('-nband','--nband',type=int,default=8,help='Jumlah tingkat energi yang ingin diplot dalam DOS.')
parser.add_argument('-occ','--occ',type=str,default='tetrahedra',help='Metode smearing yang digunakan. Default: tetrahedra. Pilihan: smearing, tetrahedra_lin, tetrahedra_opt, fixed, from_input.')
parser.add_argument("-q","--queue",type=str,default="true",help="Apakah akan menggunakan sistem antrian atau tidak. Default: true")

opt=parser.parse_args(sys.argv[1:])

# Menggunakan molekul dengan format smiles
if ('.xyz' or '.pdb' or 'POSCAR' or '.vasp' or '.poscar') in opt.input:
    geom = opt.input
elif '.mol2' in opt.input:
    if opt.queue == "true":
        with open('run_babel.sh','w') as fout:
            print("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:0:0
export OMP_NUM_THREADS={}
cd $PWD
obabel {} -O geom.xyz""".format(opt.nproc,opt.input),file=fout)
        os.system('sbatch run_babel.sh')
    else:
        with open('run_babel.sh','w') as fout:
            print("""#!/bin/bash
export OMP_NUM_THREADS={}
cd $PWD
obabel {} -O geom.xyz""".format(opt.nproc,opt.input),file=fout)
        os.system("chmod +x run_babel.sh")
        os.system('./run_babel.sh')
    geom = 'geom.xyz'
elif ('.mol2' not in opt.input and '.xyz' not in opt.input and '.pdb' not in opt.input and '.vasp' not in opt.input and '.poscar' not in opt.input and 'POSCAR' not in opt.input and '.gen' not in opt.input):
    os.system("echo '{}' > geom.smi".format(opt.input))
    if opt.queue == "true":
        with open('run_babel.sh','w') as fout:
            print("""#!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --time=168:0:0
    export OMP_NUM_THREADS={}
    cd $PWD
    obabel geom.smi -O geom.xyz --gen3d""".format(opt.nproc),file=fout)
        os.system('sbatch run_babel.sh')
    else:
        with open('run_babel.sh','w') as fout:
            print("""#!/bin/bash
    export OMP_NUM_THREADS={}
    cd $PWD
    obabel geom.smi -O geom.xyz --gen3d""".format(opt.nproc),file=fout)
        os.system("chmod +x run_babel.sh")
        os.system('./run_babel.sh')
    while not os.path.exists('geom.xyz'):
        time.sleep(3)
    geom = 'geom.xyz'
### ORCA
if opt.software == 'orca':
    orca(opt.job,opt.method,opt.nproc,geom,opt.charge,opt.mult,opt.scalefreq,opt.temp,opt.pressure,opt.nroots,opt.tda,opt.solvent,opt.constraints,opt.qmatoms,opt.totalcharge,opt.totalmult,opt.qm2method,opt.qm2basis,opt.activeatoms,opt.hessfile,opt.dispersion,opt.aim,opt.produk,opt.transitionstate,opt.irciter,opt.printlevel, opt.inithess,opt.grid, opt.finalgrid,opt.iter)

### NWChem
if opt.software == 'nwchem':
    nwchem(opt.job,opt.method,opt.nproc,geom,opt.charge,opt.mult,opt.scalefreq,opt.restart,opt.conv_thr,opt.iter)

### DCDFTBMD
if opt.software == 'dcdftb':
    dcdftb(opt.job,opt.method,geom,opt.charge,opt.mult,opt.dispersion,opt.parapath,opt.temp,opt.pressure,opt.ensembel,opt.thermostat,opt.deltat,opt.nstep,opt.mdprint,opt.a1,opt.a2,opt.a3,opt.b1,opt.b2,opt.b3,opt.c1,opt.c2,opt.c3,opt.restart,opt.traject,opt.velocity,opt.dftbinp,opt.soft,opt.softtype,opt.softrange,opt.softcenter,opt.metarest,opt.metafreq,opt.metaheight,opt.cvtype,opt.metawidth,opt.pow1,opt.pow2,opt.rcut,opt.fesstart,opt.fesend,opt.fesbin,opt.ag1,opt.ag2,opt.ag3,opt.ag4,opt.solvent,opt.nroots,opt.targetstate,opt.multtrans,opt.ocsstr,opt.writetrans,opt.longrange,opt.bufrad,opt.delta,opt.opttype,opt.freqtype,opt.econv, opt.dconv)

if opt.software == 'dock':
    if opt.job == 'readpdb':
        readpdb(opt.input)
    if opt.job == 'splitpdb':
        splitpdb(opt.input,opt.ligname)
    if opt.job == 'addH':
        addH(opt.input)
    if opt.job == 'addcharge':
        addcharge(opt.input,opt.chargetype)
    if opt.job == 'sphgen':
        sphgen(opt.input)
    if opt.queue == 'true':
        with open('run_sphgen.sh','w') as fout:
            print("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:0:0
export OMP_NUM_THREADS={}
cd $PWD
$DOCK_DIR/sphgen
$DOCK_DIR/sphere_selector protein.sph {} {}""".format(opt.nproc,opt.ligand,opt.dockrange),file=fout)
        os.system("sbatch run_sphgen.sh")
    else:
        with open('run_sphgen.sh','w') as fout:
            print("""#!/bin/bash
export OMP_NUM_THREADS={}
cd $PWD
$DOCK_DIR/sphgen
$DOCK_DIR/sphere_selector protein.sph {} {}""".format(opt.nproc,opt.ligand,opt.dockrange),file=fout)
        os.system("chmod +x run_sphgen.sh")
        os.system("./run_sphgen.sh")
    if opt.job == 'showsphere':
        showsphere()
    if opt.job == 'gridgen':
        gridgen(opt.input)
        if opt.queue == "true":
            with open('run_grid.sh','w') as fout:
                print("""#!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --time=168:0:0
    export OMP_NUM_THREADS={}
    cd $PWD
    $DOCK_DIR/grid -i grid.in""".format(opt.nproc),file=fout)
            os.system("sbatch run_grid.sh")
        else:
            with open('run_grid.sh','w') as fout:
                print("""#!/bin/bash
    export OMP_NUM_THREADS={}
    cd $PWD
    $DOCK_DIR/grid -i grid.in""".format(opt.nproc),file=fout)
            os.system("chmod +x run_grid.sh")
            os.system("./run_grid.sh")
    
    if opt.job == 'rigiddock':
        rigiddock(opt.ligand,opt.calcrmsd)
        if opt.queue == "true":
            with open('run_rigiddock.sh','w') as fout:
                print("""#!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --time=168:0:0
    export OMP_NUM_THREADS={}
    cd $PWD
    $DOCK_DIR/dock6 -i rigid.in -o rigid.out""".format(opt.nproc),file=fout)
            os.system("sbatch run_rigiddock.sh")
        else:
            with open('run_rigiddock.sh','w') as fout:
                print("""#!/bin/bash
    export OMP_NUM_THREADS={}
    cd $PWD
    $DOCK_DIR/dock6 -i rigid.in -o rigid.out""".format(opt.nproc),file=fout)
            os.system("chmod +x run_rigiddock.sh")
            os.system("./run_rigiddock.sh")
    if opt.job == 'flexdock':
        flexdock(opt.ligand,opt.calcrmsd)
        if opt.queue == "true":
            with open('run_flexdock.sh','w') as fout:
                print("""#!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --time=168:0:0
    export OMP_NUM_THREADS={}
    cd $PWD
    $DOCK_DIR/dock6 -i flex.in -o flex.out""".format(opt.nproc),file=fout)
            os.system("sbatch run_flexdock.sh")
        else:
            with open('run_flexdock.sh','w') as fout:
                print("""#!/bin/bash
    export OMP_NUM_THREADS={}
    cd $PWD
    $DOCK_DIR/dock6 -i flex.in -o flex.out""".format(opt.nproc),file=fout)
            os.system("chmod +x run_flexdock.sh")
            os.system("./run_flexdock.sh")
    if opt.job == 'translig':
        translig(opt.ligand)

if opt.job=='sdf2xyz' and opt.software == 'dock':
    sdf2xyz(opt.ligand)

if opt.job == 'multiflexdock' and opt.software == 'dock':
    multiflexdock(opt.nligands,opt.chargetype)

if opt.job == 'multiopt' and opt.software == 'dock':
    multiopt(opt.nligands)
if opt.job == 'checkopt' and opt.software == 'dock':
    checkopt(opt.nligands)

if opt.software == 'dftb':
    if '.xyz' in opt.input or '.gen' in opt.input:
        geom = opt.input
    elif 'POSCAR' in opt.input or '.poscar' in opt.input:
        poscar2gen(opt.input)
        geom = 'in.gen'
    elif 'vasp' in opt.input:
        vasp2gen(opt.input)
        geom = 'in.gen'

    dftb(geom,opt.job,opt.activeatoms,opt.method,opt.parapath,opt.dispersion,opt.kpts,opt.hcorr)
    if opt.queue == "true":
        with open('run.sh','w') as fout:
            print("""#!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --time=168:0:0
    export OMP_NUM_THREADS={}
    cd $PWD
    cp cmmd.in dftb_in.hsd
    $DFTB_COMMAND cmmd.in > cmmd.out""".format(opt.nproc),file=fout)
        os.system('sbatch run.sh')
    else:
        with open('run.sh','w') as fout:
            print("""#!/bin/bash
    export OMP_NUM_THREADS={}
    cd $PWD
    cp cmmd.in dftb_in.hsd
    $DFTB_COMMAND cmmd.in > cmmd.out""".format(opt.nproc),file=fout)
        os.system("chmod +x run.sh")
        os.system('./run.sh')


if opt.software == 'qe':
    qe(opt.input,opt.job,opt.mode,opt.pseudo,opt.outdir,opt.bravais,opt.unit,opt.ecutwfc,opt.ecutrho,opt.mixing_beta,opt.conv_thr,opt.dftfunc,opt.extpseudo,opt.kpts,opt.optalgo,opt.cellpress,opt.press_conv_thr,opt.nband,opt.occ)

# RUNNING SCRIPT
if opt.software == 'orca':
    if opt.queue=="true":
        with open('run.sh','w') as fout:
            print("""#!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task={}
    #SBATCH --time=168:0:0
    export LD_LIBRARY_PATH=/home/adit/opt/openmpi411/lib:$LD_LIBRARY_PATH
    export PATH=/home/adit/opt/openmpi411/bin:$PATH
    export OMP_NUM_THREADS=1
    cd $PWD
    $ORCA_COMMAND cmmd.in > cmmd.out --oversubscribe""".format(opt.nproc),file=fout)
        os.system('sbatch run.sh')
    else:
        with open('run.sh','w') as fout:
            print("""#!/bin/bash
    cd $PWD
    $ORCA_COMMAND cmmd.in > cmmd.out --oversubscribe""",file=fout)
        os.system("chmod +x run.sh")
        os.system('./run.sh')
if opt.software == 'nwchem':
    if opt.queue=="true":
        with open('run.sh','w') as fout:
            print("""#!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task={}
    #SBATCH --time=168:0:0
    export LD_LIBRARY_PATH=/home/adit/opt/openmpi411/lib:$LD_LIBRARY_PATH
    export PATH=/home/adit/opt/openmpi411/bin:$PATH
    export OMP_NUM_THREADS=1
    cd $PWD
    $NWCHEM_COMMAND cmmd.in > cmmd.out""".format(opt.nproc),file=fout)
        os.system('sbatch run.sh')
    else:
        with open('run.sh','w') as fout:
            print("""#!/bin/bash
    cd $PWD
    $NWCHEM_COMMAND cmmd.in > cmmd.out""",file=fout)
        os.system("chmod +x run.sh")
        os.system('./run.sh')
if opt.software == 'qe':
    with open('run.sh','w') as fout:
        print("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={}
#SBATCH --time=168:0:0
export OMP_NUM_THREADS=1
cd $PWD
mpirun -np {} $QE_COMMAND < cmmd.in > cmmd.out""".format(opt.nproc,opt.nproc),file=fout)
    os.system('sbatch run.sh')

if opt.software == 'dcdftb':
    if opt.queue == 'true':
        os.system("cp cmmd.in dftb.inp")
        with open('run.sh','w') as fout:
            print("""#!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --time=168:0:0
    export OMP_NUM_THREADS={}
    cd $PWD
    $DCDFTB_COMMAND
    mv dftb.out cmmd.out
    mv traject traject.xyz""".format(opt.nproc),file=fout)
        os.system('sbatch run.sh')
    else:
        os.system("cp cmmd.in dftb.inp")
        with open('run.sh','w') as fout:
            print("""#!/bin/bash
cd $PWD
$DCDFTB_COMMAND
mv dftb.out cmmd.out
if [ -f traject.xyz ]; then
mv traject traject.xyz
fi""",file=fout)
        os.system("chmod +x run.sh")
        os.system('./run.sh')

if opt.software == 'gromacs':
	if opt.queue == 'true':
		with open('run.sh','w') as fout:
			print("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:0:0
export OMP_NUM_THREADS=1
cd $PWD
$GROMACS_COMMAND -mt {} -mp {} -ct {} -pt {} -pp {} -l {} -cat {} -gen {} -Nump {} -prod {} -nprod {} -nequil {} -dt {} -ctype  {} -nnpt {} -comp {} -np {} -rc {} -t {}""".format(opt.terlarut,opt.pelarut,opt.c_terlarut,opt.persen_terlarut,opt.persen_pelarut,opt.lapang,opt.cation,opt.generate_dftbinp,opt.NumPelarut,opt.production,opt.nprod,opt.nequil,opt.deltat,opt.charge_type,opt.nnpt,opt.compress,opt.nproc,opt.rcol,opt.temp),file=fout)
		os.system('sbatch run.sh')
	else:
		with open('run.sh','w') as fout:
			print("""#!/bin/bash
export OMP_NUM_THREADS=1
cd $PWD
$GROMACS_COMMAND -mt {} -mp {} -ct {} -pt {} -pp {} -l {} -cat {} -gen {} -Nump {} -prod {} -nprod {} -nequil {} -dt {} -ctype  {} -nnpt {} -comp {} -np {} -rc {} -t {}""".format(opt.terlarut,opt.pelarut,opt.c_terlarut,opt.persen_terlarut,opt.persen_pelarut,opt.lapang,opt.cation,opt.generate_dftbinp,opt.NumPelarut,opt.production,opt.nprod,opt.nequil,opt.deltat,opt.charge_type,opt.nnpt,opt.compress,opt.nproc,opt.rcol,opt.temp),file=fout)
		os.system('chmod +x run.sh')
		os.system('./run.sh')		


if opt.job == 'proprep' and opt.software == 'gmx':
    proprep(opt.protein)
if opt.job == 'ligprep' and opt.software == 'gmx':
    ligprep(opt.ligand,opt.charge)

# Program XTB Standalone
if opt.software == 'xtb':
    xtb(opt.job,geom,opt.nproc,opt.produk,opt.temp,opt.nrun,opt.npoint,opt.anopt,opt.kpush,opt.kpull,opt.ppull,opt.alp,opt.distance,opt.angle,opt.dihedral,opt.scanmode,opt.iter,opt.scan,opt.solvent,opt.charge,opt.mult,opt.method,opt.fixedatoms,opt.fixedelements,opt.queue)
