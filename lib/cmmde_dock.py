import os
def readpdb(geom):
    others = []
    ligands = []
    water = []
    with open(geom, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                others.append(line.strip())
    for i in others:
        if 'HOH' in i:
            water.append(i)
        else: 
            ligands.append(i)
    ligandname = []
    for i in ligands:
        ligandname.append(i.split()[3])
    ligandname = set(ligandname)

    print("File pdb ini mengandung {} ligan, yaitu:".format(len(ligandname)))
    for i in ligandname:
        print(i)

def splitpdb(geom,ligname):
    ligand = []
    with open(geom, 'r') as f:
        for line in f:
            if line.startswith('ATOM' or 'TER'):
                protein.append(line.strip())
            if 'HETATM' in line:
                if "{}".format(ligname) in line:
                    ligand.append(line.strip())           
    
    with open('lig_{}.pdb'.format(ligname),'w') as f:
        for i in ligand:
            print(i,file=f)
## Jangan lakukan addH dan addcharge untuk molekul protein, ini khusus untuk molekul-molekul organik kecil. Untuk protein, gunakan Chimera!
def addH(geom):
    arr = geom.split(".")
    filename = arr[0]
    os.system("obabel {} -p 7.4 -O {}_withH.pdb".format(geom,filename))
    
def addcharge(geom,chargetype):
    arr = geom.split(".")
    filename = arr[0]
    os.system("obabel {} --partialcharge {} -O {}.mol2".format(geom,chargetype,filename))

def sphgen(geom):
    os.system("$DMS_COMMAND {} -n -w 1.4 -v -o protein.ms".format(geom))
    with open('INSPH','w') as f:
        print("""protein.ms
R
X
0.0
4
1.4
protein.sph""",file=f)
    if os.path.exists('OUTSPH'):
        os.remove('OUTSPH')
    if os.path.exists('protein.sph'):
        os.remove('protein.sph')
    
    # os.system("$DOCK_DIR/sphgen")
    # os.system("$DOCK_DIR/sphere_selector protein.sph {} 10".format(ligname))
        
def showsphere():
    with open('selected_spheres.in','w') as f:
        print("""selected_spheres.sph
1
N
selected_spheres.pdb
""",file=f)
    os.system('$DOCK_DIR/showsphere < selected_spheres.in')

def gridgen(geom):
    os.system('cp $DOCK_PARM/vdw_AMBER_parm99.defn .')
    with open('box.in','w') as f:
        print("""Y
5
selected_spheres.sph
1
protein_box.pdb
""",file=f)
    os.system('$DOCK_DIR/showbox < box.in')
    with open('grid.in','w') as f:
        print("""compute_grids                  yes
grid_spacing                   0.3
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  {}
box_file                       protein_box.pdb
vdw_definition_file            vdw_AMBER_parm99.defn
score_grid_prefix              grid
""".format(geom),file=f)
    # os.system('$DOCK_DIR/grid -i grid.in')

def rigiddock(ligname,rmsd):
    arr = ligname.split(".")
    filename = arr[-2]
    if os.path.exists('RigidDock_{}'.format(filename)):
        os.system('rm -rf RigidDock_{}'.format(filename))
    os.mkdir('RigidDock_{}'.format(filename))
    os.chdir('RigidDock_{}'.format(filename))
    os.system('cp $DOCK_PARM/flex.defn .')
    os.system('cp $DOCK_PARM/flex_drive.tbl .')
    with open('rigid.in','w') as f:
        print("""conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ../{}
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               {}
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ../selected_spheres.sph
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ../grid
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                ../vdw_AMBER_parm99.defn
flex_defn_file                                               flex.defn
flex_drive_file                                              flex_drive.tbl
ligand_outfile_prefix                                        rigid
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
""".format(ligname,rmsd),file=f)
    # os.system('$DOCK_DIR/dock6 -i rigid.in -o rigid.out')

def flexdock(ligname,rmsd):
    arr = ligname.split(".")
    filename = arr[-2]
    if os.path.exists('FlexDock_{}'.format(filename)):
        os.system('rm -rf FlexDock_{}'.format(filename))
    os.mkdir('FlexDock_{}'.format(filename))
    os.chdir('FlexDock_{}'.format(filename))
    os.system('cp $DOCK_PARM/flex.defn .')
    os.system('cp $DOCK_PARM/flex_drive.tbl .')
    with open('flex.in','w') as f:
        print("""conformer_search_type                                        flex
write_fragment_libraries                                     no
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              40
pruning_use_clustering                                       yes
pruning_max_orients                                          100
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               25.0
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            no
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ../{}
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               {}
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ../selected_spheres.sph
max_orientations                                             500
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ../grid
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                1000
simplex_grow_max_iterations                                  1000
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                ../vdw_AMBER_parm99.defn
flex_defn_file                                               flex.defn
flex_drive_file                                              flex_drive.tbl
ligand_outfile_prefix                                        flex
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no

""".format(ligname,rmsd),file=f)
    # os.system('$DOCK_DIR/dock6 -i flex.in -o flex.out')

def translig(ligand):
    arr = ligand.split(".")
    filename = arr[-2]
    os.system("obabel {} -O {}_center.xyz -c".format(ligand,filename))
    center_x = 0
    center_y = 0
    center_z = 0
    with open('protein_box.pdb','r') as f:
        for line in f:
            if 'CENTER' in line:
                arr = line.split()
                center_x+=float(arr[5])
                center_y+=float(arr[6])
                center_z+=float(arr[7])           
    os.system("$DOCK_DIR/transform.py {}_center.xyz -tx {} -ty {} -tz {} > {}_trans.xyz".format(filename,center_x,center_y,center_z,filename))
    os.system('obabel {}_trans.xyz -O {}.pdb'.format(filename,filename))
    os.system('rm {}_center.xyz'.format(filename))
    os.system('rm {}_trans.xyz'.format(filename))

def sdf2xyz(ligand):
    with open('run_babel.sh','w') as fout:
        with open('run_babel.sh','w') as fout:
            print("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:0:0
cd $PWD
obabel {} -O lig.xyz -m""".format(ligand),file=fout)    
    os.system('sbatch run_babel.sh')

def multiopt(nligands):
    for i in range(1,nligands+1):
        os.system('mkdir lig{}'.format(i))
        os.system('cp lig{}.xyz lig{}'.format(i,i))
        os.chdir('lig{}'.format(i))
        os.system('cmmde.py -s dcdftb -j opt -i lig{}.xyz -m DFTB3 -disp D3BJ -np 2 -solvent water'.format(i))
        os.chdir('../')

def checkopt(nligands):
    for i in range(1,nligands+1):
        os.chdir('lig{}'.format(i))
        with open('cmmd.out','r') as f:
            for line in f:
                if "Stationary point found" in line:
                    print("Optimasi geometri lig{} sukses!!!".format(i))
                    
        os.chdir('../')
def multiflexdock(nligands,chargetype):
    for i in range(1,nligands+1):
        os.chdir('lig{}'.format(i))
        os.system('cmmdepost.py -s dcdftb -j opt -i lig{}.xyz'.format(i))
        os.system('cp ../protein_box.pdb .')
        os.system('cmmde.py -s dock -j translig -ligand cmmd.xyz')
        os.system('rm protein_box.pdb')
        os.system('cmmde.py -s dock -j addcharge -i cmmd.pdb -chargetype {}'.format(chargetype))
        os.system('cp cmmd.mol2 ../lig{}.mol2'.format(i))
        os.chdir('../')
        os.system('cmmde.py -s dock -j flexdock -ligand lig{}.mol2'.format(i))

        