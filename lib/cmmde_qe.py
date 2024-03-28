from cmmde_mass import mass
def qe(geom,job,mode,pseudo,outdir,bravais,unit,ecutwfc,ecutrho,mixing,conv_thr,dftfunc,extpseudo,kpts,optalgo,press,press_conv_thr,nband,occ):
	# Mapping istilah job pada CMMDE dengan istilah job di Quantum Espresso
	jobmap = {
		'sp':'scf',
		'opt':'relax',
		'optcell':'vc-relax',
		'md':'md',
		'band':'bands',
		'mdnpt':'vc-md'
		}
	x = []
	y = []
	z = []
	symb = []
	v1 = []
	v2 = []
	v3 = []

	k_points = kpts.split("x")
	shift = []
	for k in k_points:
		if int(k)%2 == 0:
			shift.append(1)
		else:
			shift.append(0)

	if 'POSCAR' in geom or 'CONTCAR' in geom:
		with open(geom,'r') as f:
			next(f)
			next(f)
			v1.append(next(f).split())
			v2.append(next(f).split())
			v3.append(next(f).split())
			next(f)
			next(f)
			next(f)
			for line in f:
				arr = line.split()
				x.append(arr[0])
				y.append(arr[1])
				z.append(arr[2])
				symb.append(arr[3])

	if '.xyz' in geom:
		with open(geom,'r') as f:
			N_Atom = next(f)
			for line in f:
				arr = line.split()
				symb.append(arr[0])
				x.append(arr[1])
				y.append(arr[2])
				z.append(arr[3])
				

	with open("cmmd.in",'w') as f:
		# CONTROL card in input file
		print("""&CONTROL
 calculation = {},
 restart_mode = {},
 prefix = 'cmmd',
 pseudo_dir = {},
 outdir = {},
 /
&SYSTEM 
 ibrav = {},
 nat = {},
 ntyp = {},	
 ecutwfc = {},
 ecutrho = {},
 occupations = {}, 
 /
 &ELECTRONS
  mixing_beta = {},
  conv_thr = {},
 /
 """.format(jobmap[job],mode,pseudo,outdir,bravais,len(symb),len(set(symb)),ecutwfc,ecutrho,occ,mixing,conv_thr),file=f)
		if 'opt' in job:
			print("""&IONS
 ion_dynamics = {},
 /""".format(optalgo),file=f)
		if 'optcell' in job:
			print("""&CELL
 cell_dynamics = {},
 press = {},
 press_conv_thr = {},
 /""".format(optalgo,press,press_conv_thr),file=f)
		print("ATOMIC_SPECIES",file=f)
		for i in set(symb):
			print("{} {} {}.{}.{}".format(i,mass(i),i,dftfunc,extpseudo),file=f)
		if bravais == 0:
			print("""CELL_PARAMETERS ({})
 {} {} {}
 {} {} {}
 {} {} {}""".format(unit,v1[0][0],v1[0][1],v1[0][2] ,v2[0][0],v2[0][1],v2[0][2] ,v3[0][0],v3[0][1],v3[0][2] ),file=f)
		print("ATOMIC_POSITIONS ({})".format(unit),file=f)
		for i,x,y,z in zip(symb,x,y,z):
			print("{} {} {} {}".format(i,x,y,z),file=f)
		print("K_POINTS automatic",file=f)
		print("{} {} {} {} {} {}".format(k_points[0],k_points[1],k_points[2],shift[0],shift[1],shift[2]),file=f)

	if 'charge' in job:
		with open("cmmd.in",'w') as f:
			print("""&INPUTPP
 outdir = {},
 prefix = 'cmmd',
 plot_num = 0,
/
&PLOT
 iflag = 3,
 output_format = 6,
 fileout = 'cmmd_rho.cube',
 nx = 64, ny = 64, nz = 64,
 /""".format(outdir),file=f)

	if 'nscf' in job:
		with open('cmmd.in', 'w') as f:
			print("""&CONTROL
 calculation = 'nscf',
 restart_mode = {},
 prefix = 'cmmd',
 pseudo_dir = {},
 outdir = {},
/
&SYSTEM
 ibrav = {},
 nat = {},
 ntyp = {},
 ecutwfc = {},
 ecutrho = {},
 nbnd = {},
 occupations = {},
/
&ELECTRONS
 mixing_beta = {},
 conv_thr = {},
/""".format(mode,pseudo,outdir,bravais,len(symb),len(set(symb)),ecutwfc,ecutrho,nband,occ,mixing,conv_thr),file=f)
			print("ATOMIC_SPECIES",file=f)
			for i in set(symb):
				print("{} {} {}.{}.{}".format(i,mass(i),i,dftfunc,extpseudo),file=f)
			if bravais == 0:
				print("""CELL_PARAMETERS ({})
 {} {} {}
 {} {} {}
 {} {} {}""".format(unit,v1[0][0],v1[0][1],v1[0][2] ,v2[0][0],v2[0][1],v2[0][2] ,v3[0][0],v3[0][1],v3[0][2] ),file=f)
		print("ATOMIC_POSITIONS crystal",file=f)
		for i,x,y,z in zip(symb,x,y,z):
			print("{} {} {} {}".format(i,x,y,z),file=f)
		print("K_POINTS automatic",file=f)
		print("{} {} {} {} {} {}".format(k_points[0],k_points[1],k_points[2],shift[0],shift[1],shift[2]),file=f)