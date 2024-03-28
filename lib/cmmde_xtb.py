# CMMDE function for xTB standalone program
import os
def xtb(job,geom,nproc,product,temperature,nrun,npoint,anopt,kpush,kpull,ppull,alp,distance,angle,dihedral,scanmode,maxiter,scan,solvent,charge,mult,method,fixedatoms,fixedelements,slurm):
	with open('run.sh','w') as fout:
		if slurm == "true":
			print("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:0:0
export OMP_NUM_THREADS={}
cd $PWD""".format(nproc), file=fout)
		else:
			print("""#!/bin/bash
export OMP_NUM_THREADS={}""".format(nproc),file=fout)
		meth = ''
		if 'GFN-FF' in method or 'gfnff' in method:
			meth += '--gfnff'
		if 'XTB2' in method or 'xtb2' in method:
			meth += '--gfn2'
		if 'XTB1' in method or 'xtb1' in method:
			meth += '--gfn1'
		if 'XTB0' in method or 'xtb0' in method:
			meth += '--gfn0'
		if 'opt' in job and 'freq' not in job:
			if solvent != 'none':
				print("$XTB_COMMAND {} --opt -P {} --alpb {} --chrg {} --uhf {} {} > cmmd.out".format(geom,nproc,solvent,charge,mult-1,meth),file=fout)
			else:
				print("$XTB_COMMAND {} --opt -P {} --chrg {} --uhf {} {} > cmmd.out".format(geom,nproc,charge,mult-1,meth),file=fout)
		if 'fix' in job:
			if solvent != 'none':
				print("$XTB_COMMAND {} --opt -P {} --alpb {} --chrg {} --uhf {} {} --input cmmd.in > cmmd.out".format(geom,nproc,solvent,charge,mult-1,meth),file=fout)
			else:
				print("$XTB_COMMAND {} --opt -P {} --chrg {} --uhf {} {} --input cmmd.in > cmmd.out".format(geom,nproc,charge,mult-1,meth),file=fout)
		if 'freq' in job and 'opt' not in job:
			if solvent != 'none':
				print("$XTB_COMMAND {} --hess -P {} --input cmmd.in --alpb {} --chrg {} --uhf {} {} > cmmd.out".format(geom,nproc,solvent,charge,mult-1,meth), file=fout)
			else:
				print("$XTB_COMMAND {} --hess -P {} --input cmmd.in --chrg {} --uhf {} {} > cmmd.out".format(geom,nproc,charge,mult-1,meth), file=fout)
		if 'opt,freq' in job or 'opt, freq' in job:
			if solvent != 'none':
				print("$XTB_COMMAND {} --ohess -P {} --input cmmd.in --alpb {} --chrg {} --uhf {} {} > cmmd.out".format(geom,nproc,solvent,charge, mult-1,meth),file=fout)
			else:	
				print("$XTB_COMMAND {} --ohess -P {} --input cmmd.in --chrg {} --uhf {} {} > cmmd.out".format(geom,nproc,charge,mult-1,meth),file=fout)
		if 'path' in job:
			if solvent != 'none':
				print("$XTB_COMMAND {} --path {} --input cmmd.in -P {} --alpb {} --chrg {} --uhf {} ,method > cmmd.out".format(geom,product,nproc,solvent,charge,mult-1,meth),file=fout)
			else:	
				print("$XTB_COMMAND {} --path {} --input cmmd.in -P {} --chrg {} --uhf {} {} > cmmd.out".format(geom,product,nproc,charge,mult-1,meth),file=fout)
		if 'scan' in job:
			if solvent != 'none':
				print("$XTB_COMMAND {} --opt -P {} --input cmmd.in --alpb {} --chrg {} --uhf {} {} > cmmd.out".format(geom,nproc,solvent,charge,mult-1,meth), file=fout)
			else:	
				print("$XTB_COMMAND {} --opt -P {} --input cmmd.in --chrg {} --uhf {} {} > cmmd.out".format(geom,nproc,charge,mult-1,meth), file=fout)
	if 'freq' in job:
		with open('cmmd.in', 'w') as f:
			print("""$thermo
temp={}""".format(temperature),file=f)

	if 'path' in job:
		with open("cmmd.in", 'w') as f:
			print("""$path
nrun={}
   npoint={}
   anopt={}
   kpush={}
   kpull={}
   ppull={}
   alp={}
$end""".format(nrun,npoint,anopt,kpush,kpull,ppull,alp),file=f)
	if 'fix' in job:
		with open("cmmd.in",'w') as f:
			print("$fix",file=f)
			if distance != 'None':
				distance = distance.split(";")
				for index, i in enumerate(distance):
					print("distance: {}".format(i),file=f)
			if angle != 'None':
				angle = angle.split(";")
				for index, i in enumerate(angle):
					print("angle: {}".format(i),file=f)
			if dihedral != 'None':
				dihedral = dihedral.split(";")
				for index, i in enumerate(dihedral):
					print("dihedral: {}".format(i),file=f)
			if fixedatoms != 'None':
				print("atoms: {}".format(fixedatoms),file=f)
			if fixedelements != 'None':
				print("elements: {}".format(fixedelements),file=f)

	if 'scan' in job:
		with open("cmmd.in",'w') as f:
			print("$constrain",file=f)
			if distance != 'None':
				distance = distance.split(";")
				for index, i in enumerate(distance):
					print("distance: {}".format(i),file=f)
			if angle != 'None':
				angle = angle.split(";")
				for index, i in enumerate(angle):
					print("angle: {}".format(i),file=f)
			if dihedral != 'None':
				dihedral = dihedral.split(";")
				for index, i in enumerate(dihedral):
					print("dihedral: {}".format(i),file=f)
			print("$scan",file=f)
			if scanmode != 'None':
				print("mode = {}".format(scanmode),file=f)
			scan = scan.split(";")
			for index, i in enumerate(scan):
				print("{}: {}".format(index+1,i),file=f)
			print("""$opt
  maxcycle = {}""".format(maxiter),file=f)

	if slurm == "true":
		os.system('sbatch run.sh')
	else:
		os.system('chmod +x run.sh')
		os.system('./run.sh')
	

	

	
