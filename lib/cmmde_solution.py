#!/usr/bin/env python3

import subprocess
import os
import shutil
import gmxscript
from scipy.constants import N_A
import argparse
import sys
import time
from cmmde_mass import mass
parser = argparse.ArgumentParser(description='SPE Solution Program')
parser.add_argument('-mt', '--terlarut', type=str, help='Nama file molekul terlarut dalam format .xyz (ditulis tanpa ekstensi).')
parser.add_argument('-ctype','--charge_type',type=str, default='gas', help='Tipe muatan yang digunakan untuk parameterisasi muatan. Pilihan: gas, bcc, dftb')
parser.add_argument('-ct', '--c_terlarut', type=str, help='Muatan bersih molekul terlarut. Jika terlarut lebih dari satu buah, pisahkan dengan koma')
parser.add_argument('-pt','--persen_terlarut',type=str, help='Fraksi massa terlarut dalam persen. Jika pelarut lebih dari satu buah, pisahkan menggunakan koma.')
parser.add_argument('-pp','--persen_pelarut',type=float, help='Fraksi massa pelarut dalam persen.')
parser.add_argument('-mp', '--pelarut', type=str, help='Nama file molekul pelarut dalam format .xyz (ditulis tanpa ekstensi).')
parser.add_argument('-cp', '--c_pelarut', type=int, help='Muatan bersih molekul pelarut',default=0)
parser.add_argument('-Nump','--NumPelarut', type=int, default=100, help='Jumlah molekul pelarut maksimum dalam sistem larutan. Default = 100.')
parser.add_argument('-cat','--cation',type=str, default='none',help='Kation yang digunakan dalam sistem elektrolit baterai.')
parser.add_argument('-p', '--pressure', default=1.0, type=float, help='Tekanan dalam satuan bar')
parser.add_argument('-t', '--temperature', default=298.15, type=float, help='Suhu yang digunakan dalam simulasi.')
parser.add_argument('-nequil','--nequil',type=int, default=50000, help='Jumlah step pada saat ekuilibrasi.')
parser.add_argument('-nnpt','--nnpt',type=int, default=50000, help='Jumlah step pada saat ekuilibrasi NPT.')
parser.add_argument('-nprod','--nprod',type=int, default=400000, help='Jumlah step pada saat production.')
parser.add_argument('-dt','--dt',type=float, default=1, help='Step simulasi dinamika molekul yang dilakukan.')
parser.add_argument('-comp','--compress',type=float, default=4.5e-6, help='Nilai kompresibilitas. Default=4.5e-6.')
parser.add_argument('--packmol', type=str, default='/Users/adit/opt/packmol/packmol')
parser.add_argument('-l','--lapang',type=float, default=5, help='Ruang lapang untuk atom-atom bergerak di dalam kotak larutan (angstrom).')
parser.add_argument('-gen','--generate_dftbinp',type=str,default='false',help='Apakah ingin mengkonversi ke dalam format koordinat xyz?')
parser.add_argument('-prod','--production',type=str,default='None',help='Jenis production run yang ingin dilakukan. Pilihan: NPT (recommended), NVE, dan NVT.')
parser.add_argument('-np','--nproc',type=int,default=1, help='Jumlah CPU yang digunakan')
parser.add_argument('-rc','--rcol',type=float,default=0, help='Additional Coulombic radius for PME calculation')

opt=parser.parse_args(sys.argv[1:])



# Input file minimisasi energi dalam GROMACS
minim_mdp = """
integrator  = cg         ; Algorithm (steep or cg)
emtol       = 1        ; Stop minimization when the maximum force < 1 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 5         ; Frequency to update the neighbor list and long range forces
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = {}       ; Short-range electrostatic cut-off
rvdw            = {}       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
;define = -DFLEXIBLE
"""
nve_mdp = """
;define                  = -DFLEXIBLE
integrator              = md        
nsteps                  = {}    ; 1 * 500000 = 50 ps
dt                      = 0.001     ; 0.5 fs
nstxout                 = 100       ; save coordinates every 1.0 ps
nstenergy               = 100       ; save energies every 1.0 ps
nstlog                  = 100       ; update log file every 1.0 ps

; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = {}       ; short-range electrostatic cutoff (in nm)
rvdw                    = {}       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is off
tcoupl                  = no             ; no temperature coupling in NVE
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVE
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; Velocity generation is off
"""

# Input file simulasi menggunakan ensembel kanonik (NVT) dalam GROMACS
nvt_mdp="""
;define                  = -DFLEXIBLE
integrator              = md        ; leap-frog integrator
nsteps                  = {}    
dt                      = 0.001     ; 1 fs
nstxout                 = 100       ; save coordinates every 1.0 ps
;nstenergy               = 100       ; save energies every 1.0 ps
;nstlog                  = 100       ; update log file every 1.0 ps

; Bond parameters
continuation            = yes        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = {}       ; short-range electrostatic cutoff (in nm)
rvdw                    = {}       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = System   ; two coupling groups - more accurate
tau_t                   = 0.01               ; time constant, in ps
ref_t                   = {}              ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no       ; assign velocities from Maxwell distribution
"""

# Input file simulasi isoterm-isobarik (NPT) dalam GROMACS. Harapannya, setelah dilakukan NPT, akan diperoleh sistem dengan massa jenis yang optimal.
npt_mdp = """
title   = NPT equilibration
; Run parameters
integrator  = md    ; leap-frog integrator
nsteps    = {}     
dt        = 0.001   ; 1 fs
; Output control
nstxout   = 500   ; save coordinates every 1.0 ps
nstvout   = 500   ; save velocities every 1.0 ps
nstenergy = 500   ; save energies every 1.0 ps
nstlog    = 500   ; update log file every 1.0 ps
; Bond parameters
continuation            = yes        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 6         ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type       = grid    ; search neighboring grid cells
nstlist       = 10      ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb      = {}    ; short-range electrostatic cutoff (in nm)
rvdw        = {}    ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype     = PME   ; Particle Mesh Ewald for long-range electrostatics
pme_order     = 4       ; cubic interpolation
fourierspacing  = 0.16    ; grid spacing for FFT
; Temperature coupling is on
tcoupl    = V-rescale             ; modified Berendsen thermostat
tc-grps     = System   ; two coupling groups - more accurate
tau_t   = 0.1             ; time constant, in ps
ref_t   = {}              ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl            = Berendsen     ; Pressure coupling on in NPT
pcoupltype          = isotropic             ; uniform scaling of box vectors
tau_p           = 2               ; time constant, in ps
ref_p           = {}                ; reference pressure, in bar
compressibility     = {}              ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc   = xyz   ; 3-D PBC
; Dispersion correction
DispCorr  = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel   = no    ; Velocity generation is off
"""

nptp_mdp = """
title       = NPT Production
; Run parameters
integrator  = md        ; leap-frog integrator
nsteps      = {}        ;
dt          = {}     ; in ps
; Output control
nstxout     = 500       ; save coordinates every 1.0 ps
nstvout     = 500       ; save velocities every 1.0 ps
nstenergy   = 500       ; save energies every 1.0 ps
nstlog      = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes        ; if no => first dynamics run
; Neighborsearching
cutoff-scheme   = Verlet
ns_type         = grid      ; search neighboring grid cells
nstlist         = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb        = {}        ; short-range electrostatic cutoff (in nm)
rvdw            = {}        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.16      ; grid spacing for FFT
; Temperature coupling is off
tcoupl    = V-rescale             ; modified Berendsen thermostat
tc-grps     = System   ; two coupling groups - more accurate
tau_t   = 0.1             ; time constant, in ps
ref_t   = {}              ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl            = Berendsen     ; Pressure coupling on in NPT
pcoupltype        = isotropic             ; uniform scaling of box vectors
tau_p           = 10               ; time constant, in ps
compressibility   = {}              ; isothermal compressibility of water, bar^-1
ref_p           = {}                ; reference pressure, in bar
refcoord_scaling    = com
; Periodic boundary conditions
pbc     = xyz       ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel     = no        ; Velocity generation is off
"""

nvtp_mdp = """
title       = NVT Production
; Run parameters
integrator  = md        ; leap-frog integrator
nsteps      = {}        ;
dt          = {}     ; in ps
; Output control
nstxout     = 500       ; save coordinates every 1.0 ps
nstvout     = 500       ; save velocities every 1.0 ps
nstenergy   = 500       ; save energies every 1.0 ps
nstlog      = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes        ; if no => first dynamics run
; Neighborsearching
cutoff-scheme   = Verlet
ns_type         = grid      ; search neighboring grid cells
nstlist         = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb        = {}        ; short-range electrostatic cutoff (in nm)
rvdw            = {}        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.16      ; grid spacing for FFT
; Temperature coupling is off
tcoupl    = V-rescale             ; modified Berendsen thermostat
tc-grps     = System   ; two coupling groups - more accurate
tau_t   = 0.1             ; time constant, in ps
ref_t   = {}              ; reference temperature, one for each group, in K
; Pressure coupling is off in NVT simulation
pcoupl            = no     ; Pressure coupling on in NPT

; Periodic boundary conditions
pbc     = xyz       ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel     = no        ; Velocity generation is off
"""

nvep_mdp = """
title       = NVE Production
; Run parameters
integrator  = md        ; leap-frog integrator
nsteps      = {}        ;
dt          = {}     ; in ps
; Output control
nstxout     = 500       ; save coordinates every 1.0 ps
nstvout     = 500       ; save velocities every 1.0 ps
nstenergy   = 500       ; save energies every 1.0 ps
nstlog      = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes        ; if no => first dynamics run
; Neighborsearching
cutoff-scheme   = Verlet
ns_type         = grid      ; search neighboring grid cells
nstlist         = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb        = {}        ; short-range electrostatic cutoff (in nm)
rvdw            = {}        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.16      ; grid spacing for FFT
; Temperature coupling is off
tcoupl    = no             ; modified Berendsen thermostat

; Pressure coupling is off in NVE simulation
pcoupl            = no     ; Pressure coupling on in NPT


; Periodic boundary conditions
pbc     = xyz       ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel     = no        ; Velocity generation is off
"""

def chargeparm(mol2file,dftbdat):
    header = []
    footer = []
    serial = []
    symbol = []
    x = []
    y = []
    z = []
    type = []
    resid = []
    resname = []
    charge = []
    with open(mol2file, 'r') as f:
        for i in range(1,9):
            header.append(next(f).strip())
        for line in range(1,int(header[2].split()[0])+1):
            arr = next(f).split()
            serial.append(arr[0])
            symbol.append(arr[1])
            x.append(arr[2])
            y.append(arr[3])
            z.append(arr[4])
            type.append(arr[5])
            resid.append(arr[6])
            resname.append(arr[7])
        for line in f:
            footer.append(line.strip())
    with open(dftbdat,'r') as f:
        for line in f:
            if 'MULLIKEN' in line:
                next(f)
                for i in range(1,int(header[2].split()[0])+1):
                    arr = next(f).split()
                    charge.append(arr[2])

    fname = mol2file.split(".")[0]
    with open("{}_mod.mol2".format(fname),'w') as f:
        for line in header[0:4]:
            print(line,file=f)
        print("DFTBCharge",file=f)
        print("""

    @<TRIPOS>ATOM""",file=f)
        for serial, symbol, x, y, z, type, resid, resname, charge in zip(serial, symbol, x, y, z, type, resid, resname, charge):
            print("     {} {}         {}     {}    {} {}         {} {}      {}".format(serial, symbol, x, y, z, type, resid, resname, charge), file=f)
        for line in footer:
            print(line, file=f)
def make_solution(terlarut,c_terlarut,pelarut,c_pelarut,N_pelarut,persen_terlarut,persen_pelarut,cation, temperature,pressure,packmol,charge_type,lapang,generate_dftbinp,prod,nproc):

  if temperature:
    print (f'Suhu: {temperature:12.4f} K')

# Penyiapan box molekul
  terlarut = terlarut.split(",")
  c_terlarut = c_terlarut.split(",")
  c_terlarut = dict(zip(terlarut,c_terlarut))
  for i in terlarut:
    os.chdir('{}'.format(i))
    if charge_type == 'dftb':
      os.system("obabel cmmd.xyz -O {}.pdb".format(i))
    else:
      os.system("rm geom.xyz")
      os.system("rm cmmd_trj.xyz")
      os.system("obabel *.xyz -O {}.pdb".format(i))
    with open('{}.pdb'.format(i),'r') as f:
      fdata = f.read()
      fdata = fdata.replace('UNL','{}'.format((i[0]+i[1]+i[-1]).upper()))
      fdata = fdata.replace('CONECT', '#CONECT')
      with open('{}.pdb'.format(i),'w') as f:
        f.write(fdata)


    if opt.charge_type == 'dftb':
        charge_type='gas'
        os.system("antechamber -j 5 -at gaff -dr no -i {}.pdb -fi pdb -o {}.mol2 -fo mol2 -c {} -s 2 -nc {}".format(i,i,charge_type,float(c_terlarut[i])))
        chargeparm("{}.mol2".format(i), "dftb.dat")
    else:
        os.system("antechamber -j 5 -at gaff -dr no -i {}.pdb -fi pdb -o {}.mol2 -fo mol2 -c {} -s 2 -nc {}".format(i,i,charge_type,float(c_terlarut[i])))
        os.system("cp {}.mol2 {}_mod.mol2".format(i,i))
    os.system("parmchk2 -i {}_mod.mol2 -f mol2 -o {}.frcmod".format(i,i))
    with open("tleap.in", 'w') as fout:
        print("""
    source leaprc.gaff
        {} = loadmol2 {}_mod.mol2
        loadamberparams {}.frcmod
        check {}
        saveoff {} {}.lib
        saveamberparm {} {}.parmtop {}.inpcrd
        savepdb {} {}_gaff.pdb
        quit
    """.format((i[0]+i[1]+i[-1]).upper(), i, i, (i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper()), file=fout)
    os.system("tleap -f tleap.in")
    os.chdir('../')
  os.chdir('{}'.format(pelarut))

# Penyiapan molekul pelarut menggunakan GAFF
# Memanggil antechamber untuk membuatkan file mol2 yang mengandung informasi muatan bersih atom

  if charge_type == 'dftb':
    os.system("obabel cmmd.xyz -O {}.pdb".format(pelarut))
  else:
    os.system("rm geom.xyz")
    os.system("rm cmmd_trj.xyz")
    os.system("obabel *.xyz -O {}.pdb".format(pelarut))
  with open('{}.pdb'.format(pelarut),'r') as f:
    fdata = f.read()
  fdata = fdata.replace('UNL','SLV')
  fdata = fdata.replace('HOH','SLV')
  fdata = fdata.replace('CONECT', '#CONECT')
  with open('{}.pdb'.format(pelarut),'w') as f:
    f.write(fdata)

  if opt.charge_type == 'dftb':
    charge_type='gas'
    os.system("antechamber -j 5 -at gaff -dr no -i {}.pdb -fi pdb -o {}.mol2 -fo mol2 -c {} -s 2 -nc {}".format(pelarut,pelarut,charge_type,c_pelarut))
    chargeparm("{}.mol2".format(pelarut), "dftb.dat")
# Mengidentifikasi parameter yang tidak terdefinisikan dalam GAFF
  else:
    os.system("antechamber -j 5 -at gaff -dr no -i {}.pdb -fi pdb -o {}.mol2 -fo mol2 -c {} -s 2 -nc {}".format(pelarut,pelarut,charge_type,c_pelarut))
    os.system("cp {}.mol2 {}_mod.mol2".format(pelarut,pelarut))

  os.system("parmchk2 -i {}_mod.mol2 -f mol2 -o {}.frcmod".format(pelarut,pelarut))

  with open("tleap.in", 'w') as fout:
    print("""
    source leaprc.gaff
        SLV = loadmol2 {}_mod.mol2
        loadamberparams {}.frcmod
        check SLV
        saveoff SLV {}.lib
        saveamberparm SLV {}.parmtop {}.inpcrd
        savepdb SLV {}_gaff.pdb
        quit
  """.format(pelarut,pelarut,pelarut,pelarut,pelarut,pelarut), file=fout)
  os.system("tleap -f tleap.in")

  if cation != 'none':
      os.system("mkdir ../{}".format(cation))
      os.chdir('../{}'.format(cation))
      with open('{}.mol2'.format(cation),'w') as f:
          print("""@<TRIPOS>MOLECULE
ION
 1 0 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 {}+          0.0000    0.0000    0.0000 {}+      1  ion        1.0000
@<TRIPOS>BOND
    """.format(cation,cation),file=f)

# Penyiapan box
  os.chdir('../')

  with open ("{}/cmmd.xyz".format(pelarut), 'r') as f:
      elements = []
      Mr_pelarut = 0
      next(f)
      next(f)
      for line in f:
          arr = line.split()
          elements.append(arr[0])
      for ele in elements:
          Mr_pelarut += float(mass(ele))
  Mr_terlarut = []
  for i in terlarut:
      ele_terlarut = []
      Mr_t = 0
      with open("{}/cmmd.xyz".format(i), 'r') as f:
          next(f)
          next(f)
          for line in f:
              arr = line.split()
              ele_terlarut.append(arr[0])
          for ele in ele_terlarut:
              Mr_t += float(mass(ele))
      Mr_terlarut.append(Mr_t)

  Mr_terlarut = dict(zip(terlarut,Mr_terlarut))

  m_pelarut = N_pelarut*(Mr_pelarut/N_A) # dalam g

  m_larutan = m_pelarut/persen_pelarut*100
  persen_terlarut = persen_terlarut.split(",")
  persen_terlarut = dict(zip(terlarut, persen_terlarut))
  N_terlarut = []
  for i in terlarut:
      N = round((float(persen_terlarut[i])/100)*m_larutan/Mr_terlarut[i]*N_A)
      N_terlarut.append(N)
  N_terlarut = dict(zip(terlarut,N_terlarut))
  V_box = m_larutan # Dalam satuan cm3
# Panjang rusuk kubus dalam satuan cm:
  box_size = V_box**(1/3)
# Panjang rusuk kubus dalam satuan Angstroem:
  cmToAng = 1e8
  box_size = box_size * cmToAng + lapang
# Hentikan perhitungan jika ukuran kubus terlalu kecil
  half_box = min(1.0, box_size/20.0-0.1)
  print(half_box)

  def cmmde_solution():
      cwd = os.path.abspath(os.curdir)
      print('Membuat sistem larutan menggunakan PACKMOL')
      os.makedirs('SistemLarutan', exist_ok=True)
      os.chdir('SistemLarutan')


      with open('packmol.inp', 'w') as fout:
          print("""tolerance 2.5
filetype pdb
output system_init.pdb

structure ../{}/{}_gaff.pdb
    number {}
    inside cube -{} -{} -{} {}
    resnumbers 3
  end structure""".format(pelarut, pelarut, N_pelarut, box_size/2.0, box_size/2.0, box_size/2.0, box_size), file=fout)
          for i in terlarut:
            print("""structure ../{}/{}_gaff.pdb
    number {}
    inside cube -{} -{} -{} {}
    resnumbers 3
  end structure""".format(i,(i[0]+i[1]+i[-1]).upper(),N_terlarut[i],box_size/2.0, box_size/2.0, box_size/2.0, box_size), file=fout)
      print('='*80)

      with open ('tleap.in', 'w') as fout:
          print("""source leaprc.gaff
loadamberparams /Users/adit/opt/amber22/dat/leap/parm/frcmod.ionslm_iod_opc3""", file=fout)
          for i in terlarut:
            print("""loadoff ../{}/{}.lib
loadamberparams ../{}/{}.frcmod""".format(i,(i[0]+i[1]+i[-1]).upper(),i,i),file=fout)
          print("""loadoff ../{}/{}.lib
loadamberparams ../{}/{}.frcmod
SYSTEM = loadpdb system_init.pdb""".format(pelarut,pelarut,pelarut,pelarut),file=fout)
          if cation != 'none':
            print("""ION = loadmol2 ../{}/{}.mol2
addions SYSTEM ION 0""".format(cation,cation),file=fout)
          print("""list
saveamberparm SYSTEM system.prmtop system.inpcrd
savepdb SYSTEM system.pdb
quit""",file=fout)


      with open('packmol.inp', 'r') as f:
          print(f.read())
      print('='*80)

      os.system("{} < packmol.inp".format(packmol))


      os.system("tleap -f tleap.in")
      os.system("acpype -p system.prmtop -x system.inpcrd")

      gmxscript.editconf(
          f = "system.pdb",
          box = [box_size/10.0, box_size/10.0, box_size/10.0],
          bt = "cubic",
          o = "system_box.pdb"
      )

  # Menuliskan semua input file yang dibutuhkan, meliputi simulasi NVE, NVT, dan NPT
      nequil = opt.nequil
      nnpt = opt.nnpt
      nprod = opt.nprod
      compress = opt.compress
      rcoulomb = half_box+opt.rcol
      rvdw = half_box
      with open('minim.mdp', 'w') as f:
          print(minim_mdp.format(rcoulomb, rvdw), file=f)
      with open('nve.mdp', 'w') as f:
          print(nve_mdp.format(nequil,rcoulomb,rvdw), file=f)
      with open('nvt.mdp', 'w') as f:
          print(nvt_mdp.format(nequil,rcoulomb, rvdw,temperature), file=f)
      with open('npt.mdp', 'w') as f:
        print(npt_mdp.format(nnpt,rcoulomb, rvdw,temperature,pressure,compress), file=f)
      os.system("cp system.amb2gmx/system_GMX.* .")
      gmxscript.grompp(f="minim.mdp", c="system_box.pdb", o="mini.tpr", p="system_GMX.top", maxwarn=2)
      gmxscript.mdrun(deffnm="mini",ntmpi='{}'.format(opt.nproc))
      gmxscript.select(f="mini.trr", s="mini.tpr", on="mini", select="System")
      gmxscript.grompp(f="nve.mdp", c="mini.gro", p="system_GMX.top", o="nve.tpr",maxwarn=2)
      gmxscript.mdrun(deffnm="nve", v=True, ntmpi='{}'.format(opt.nproc))
      gmxscript.grompp(f="nvt.mdp", c="nve.gro", p="system_GMX.top", o="nvt.tpr",maxwarn=2)
      gmxscript.mdrun(deffnm="nvt", v=True, ntmpi='{}'.format(opt.nproc))
      gmxscript.grompp(f="npt.mdp", c="nvt.gro", p="system_GMX.top", t="nvt.cpt",o="npt.tpr",maxwarn=2)
      gmxscript.mdrun(deffnm="npt", v=True, ntmpi='{}'.format(opt.nproc))

      if prod == 'npt':
        with open('nptp.mdp', 'w') as f:
          print(nptp_mdp.format(nprod,opt.dt/2000,rcoulomb,rvdw,temperature,compress,pressure), file=f)
        gmxscript.grompp(f="nptp.mdp", c="npt.gro", p="system_GMX.top", o="nptp.tpr",maxwarn=2)
        gmxscript.mdrun(deffnm="nptp", v=True, ntmpi='{}'.format(opt.nproc))
        gmxscript.trjconv(f="nptp.trr", s="nptp.tpr", n="mini.ndx", o="finalsystem.pdb")
        print("Melakukan perhitungan MSD untuk pelarut")
        os.system("""gmx select -f nptp.gro -s nptp.tpr -on SLV.ndx -select "(resname SLV)" """)
        os.system("gmx msd -f nptp.trr -s nptp.tpr -n SLV.ndx -o msd_SLV.xvg")
        print("Melakukan perhitungan MSD untuk terlarut, masukkan masing-masing terlarut")
        for i in terlarut:
          os.system("""gmx select -f nptp.gro -s nptp.tpr -on {}.ndx -select "(resname {})" """.format((i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper()))
          os.system("gmx msd -f nptp.trr -s nptp.tpr -n {}.ndx -o msd_{}.xvg".format((i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper()))

      if prod == 'nvt':
        with open('nvtp.mdp', 'w') as f:
          print(nvtp_mdp.format(nprod,opt.dt/2000,rcoulomb,rvdw,temperature), file=f)
        gmxscript.grompp(f="nvtp.mdp", c="npt.gro", p="system_GMX.top", o="nvtp.tpr",maxwarn=2)
        gmxscript.mdrun(deffnm="nvtp", v=True, ntmpi='{}'.format(opt.nproc))
        gmxscript.trjconv(f="nvtp.trr", s="nvtp.tpr", n="mini.ndx", o="finalsystem.pdb")
        print("Melakukan perhitungan MSD untuk pelarut")
        os.system("""gmx select -f nvtp.gro -s nvtp.tpr -on SLV.ndx -select "(resname SLV)" """)
        os.system("gmx msd -f nvtp.trr -s nvtp.tpr -n SLV.ndx -o msd_SLV.xvg")
        print("Melakukan perhitungan MSD untuk terlarut, masukkan masing-masing terlarut")
        for i in terlarut:
          os.system("""gmx select -f nvtp.gro -s nvtp.tpr -on {}.ndx -select "(resname {})" """.format((i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper()))
          os.system("gmx msd -f nvtp.trr -s nvtp.tpr -n {}.ndx -o msd_{}.xvg".format((i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper()))
      if prod == 'nve':
        with open('nvep.mdp', 'w') as f:
          print(nvep_mdp.format(nprod,opt.dt/2000,rcoulomb,rvdw), file=f)
        gmxscript.grompp(f="nvep.mdp", c="npt.gro", p="system_GMX.top", o="nvep.tpr",maxwarn=2)
        gmxscript.mdrun(deffnm="nvep", v=True, ntmpi='{}'.format(opt.nproc))
        gmxscript.trjconv(f="nvep.trr", s="nvep.tpr", n="mini.ndx", o="finalsystem.pdb")
        print("Melakukan perhitungan MSD untuk pelarut")
        os.system("""gmx select -f nvep.gro -s nvep.tpr -on SLV.ndx -select "(resname SLV)" """)
        os.system("gmx msd -f nvep.trr -s nvep.tpr -n SLV.ndx -o msd_SLV.xvg")
        print("Melakukan perhitungan MSD untuk terlarut, masukkan masing-masing terlarut")
        for i in terlarut:
          os.system("""gmx select -f nvep.gro -s nvep.tpr -on {}.ndx -select "(resname {})" """.format((i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper()))
          os.system("gmx msd -f nvep.trr -s nvep.tpr -n {}.ndx -o msd_{}.xvg".format((i[0]+i[1]+i[-1]).upper(),(i[0]+i[1]+i[-1]).upper()))

      if generate_dftbinp != 'false':
          os.system("obabel finalsystem.pdb -O finalsystem.xyz")
          with open('finalsystem.xyz', 'r') as f:
              fdata = f.read()
          fdata = fdata.replace('*','{}'.format(opt.cation))
          with open('finalsystem.pdb', 'r') as f:
              tv1 = 0
              tv2 = 0
              tv3 = 0
              for line in f:
                  if 'CRYST' in line:
                      arr = line.split()
                      tv1 += float(arr[1])
                      tv2 += float(arr[2])
                      tv3 += float(arr[3])
                      break
          with open('finalsystem.xyz'.format(terlarut),'w') as fout:
              fout.write(fdata)
              print("TV {} 0.00 0.00".format(tv1),file=fout)
              print("TV 0.00 {} 0.00".format(tv2), file=fout)
              print("TV 0.00 0.00 {}".format(tv3), file=fout)

      os.chdir(cwd)

  shutil.rmtree('SistemLarutan', ignore_errors=True)
  cmmde_solution()
# Eksekusi simulasi
make_solution(opt.terlarut,opt.c_terlarut,opt.pelarut,opt.c_pelarut,opt.NumPelarut,opt.persen_terlarut,opt.persen_pelarut,opt.cation,opt.temperature,opt.pressure,opt.packmol,opt.charge_type,opt.lapang,opt.generate_dftbinp,opt.production,opt.nproc)
