import os
def proprep(protein):
    filename = protein.split('.pdb')[0]
    os.system('$GROMACS_DIR/gmx pdb2gmx -f {} -o {}.gro'.format(protein,filename))

def ligprep(ligand,charge):
    if '.mol2' in ligand:
        filename = ligand.split('.mol2')[0]
        os.system('obabel {} -O temp_{}_h.pdb'.format(ligand,filename))
        
    elif '.xyz' in ligand:
        filename = ligand.split('.xyz')[0]
        os.system('obabel {} -O temp_{}_h.pdb'.format(ligand,filename))

    elif '.pdb' in ligand:
        filename = ligand.split('.pdb')[0]
        os.system('reduce {}.pdb > temp_{}_h.pdb'.format(filename,filename))
    
    delete_entries = ['CONECT']
    with open('temp_{}_h.pdb'.format(filename),'r') as oldpdb, open('{}_h.pdb'.format(filename),'w') as newpdb:
        for line in oldpdb:
            if not any(del_entry in line for del_entry in delete_entries):
                newpdb.write(line)
    os.system('rm temp_{}_h.pdb'.format(filename))
    
    os.system('antechamber -j 5 -at gaff -dr no -i {}_h.pdb -fi pdb -o {}_h.mol2 -fo mol2 -c gas -s 2 -nc {}'.format(filename,filename,charge))
    
    
    os.system('parmchk2 -i {}_h.mol2 -f mol2 -o {}_h.frcmod'.format(filename,filename))
    with open('tleap.in','w') as f:
        print("""source oldff/leaprc.ff99SB
source leaprc.gaff
LIG = loadmol2 {}_h.mol2 
loadamberparams {}_h.frcmod
check LIG 
saveoff LIG {}_h.lib 
saveamberparm LIG {}_h.prmtop {}_h.inpcrd 
quit""".format(filename,filename,filename,filename,filename),file=f)
    os.system('tleap -f tleap.in')
    os.system('acpype -p {}_h.prmtop -x {}_h.inpcrd'.format(filename,filename))