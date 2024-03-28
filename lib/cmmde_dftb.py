from cmmde_hubbard import hubbard, azimuth
def xyz2gen(geom,a1,a2,a3,b1,b2,b3,c1,c2,c3):
    sym = []
    x = []
    y = []
    z = []
    elements = []
    with open(geom,'r') as f:
        Natom = int(next(f))
        next(f)
        for line in f:
            arr = line.split()
            sym.append(arr[0])
            x.append(float(arr[1]))
            y.append(float(arr[2]))
            z.append(float(arr[3]))
    for index,element in enumerate(sym):
        if sym[index] != sym[index-1]:
            elements.append(element)

    with open("in.gen", 'w') as f:
        types = ""
        if a1 != 0:
            types += "S"
        else:
            types += "C"
        print("{} {}".format(Natom,types),file=f)
        for i in elements:
            print(i, end=' ',file=f)
        print("",file=f)
        indx = 1
        sindx = 0
        sym_indx = []
        for i,symbol in enumerate(sym):
            if sym[i] != sym[i-1]:
                sindx+=1
                sym_indx.append(sindx)
            else:
                sindx+=0
                sym_indx.append(sindx)
        for sym_indx,i,j,k in zip(sym_indx,x,y,z):
            print(indx,sym_indx,i,j,k,file=f)
            indx+=1
        if 'S' in types:
            print("0 0 0",file=f)
            print("{} {} {}".format(a1,a2,a3),file=f)
            print("{} {} {}".format(b1,b2,b3),file=f)
            print("{} {} {}".format(c1,c2,c3),file=f)

def poscar2gen(geom):
    sym = []
    a = []
    b = []
    c = []
    x = []
    y = []
    z = []
    coord_type = []
    with open(geom, 'r') as f:
        next(f)
        next(f)
        acell = next(f).split()
        bcell = next(f).split()
        ccell = next(f).split()
        a.append(acell)
        b.append(bcell)
        c.append(ccell)
        next(f)
        next(f)
        if 'direct' in next(f):
            coord_type.append("F")
        else:
            coord_type.append("S")
        for line in f:
            arr = line.split()
            x.append(arr[0])
            y.append(arr[1])
            z.append(arr[2])
            sym.append(arr[3])
    with open('in.gen','w') as f:
        elements = []
        print("{} {}".format(len(sym),coord_type[-1]),file=f)
        for index,element in enumerate(sym):
            if sym[index] != sym[index-1]:
                elements.append(element)
        for i in elements:
            print(i, end=' ',file=f)
        print("",file=f)
        indx = 1
        sindx = 0
        sym_indx = []
        for i,symbol in enumerate(sym):
            if sym[i] != sym[i-1]:
                sindx+=1
                sym_indx.append(sindx)
            else:
                sindx+=0
                sym_indx.append(sindx)
        for sym_indx,i,j,k in zip(sym_indx,x,y,z):
            print(indx,sym_indx,i,j,k,file=f)
            indx+=1
        if 'S' or 'F' in coord_type:
            print("0 0 0",file=f)
            print("{} {} {}".format(a[0][0],a[0][1],a[0][2]),file=f)
            print("{} {} {}".format(b[0][0],b[0][1],b[0][2]),file=f)
            print("{} {} {}".format(c[0][0],c[0][1],c[0][2]),file=f)

def vasp2gen(geom):
    sym = []
    a = []
    b = []
    c = []
    x = []
    y = []
    z = []
    coord_type = []
    with open(geom, 'r') as f:
        next(f)
        next(f)
        acell = next(f).split()
        bcell = next(f).split()
        ccell = next(f).split()
        a.append(acell)
        b.append(bcell)
        c.append(ccell)
        next(f)
        next(f)
        # next(f)
        if 'direct' in next(f):
            coord_type.append("F")
        else:
            coord_type.append("S")
        for line in f:
            arr = line.split()
            x.append(arr[0])
            y.append(arr[1])
            z.append(arr[2])
            sym.append(arr[3])
            
    with open('in.gen','w') as f:
        elements = []
        print("{} {}".format(len(sym),coord_type[-1]),file=f)
        if len(list(set(sym))) == 1:
            for i in list(set(sym)):
                elements.append(i) 
        else: 
            for index,element in enumerate(sym):
                sym[index] != sym[index-1]
                elements.append(element)
        
        for i in elements:
            print(i, end=' ',file=f)
        print("",file=f)
        indx = 1
        sindx = 0
        sym_indx = []
        for i,symbol in enumerate(sym):
            if len(set(sym)) == 1:
                 sindx=1
                 sym_indx.append(sindx)
            else:
                if sym[i] != sym[i-1]:
                    sindx+=1
                    sym_indx.append(sindx)
                else:
                    sindx+=0
                    sym_indx.append(sindx)
        for sym_indx,i,j,k in zip(sym_indx,x,y,z):
            print(indx,sym_indx,i,j,k,file=f)
            indx+=1
        if 'S' or 'F' in coord_type:
            print("0 0 0",file=f)
            print("{} {} {}".format(a[0][0],a[0][1],a[0][2]),file=f)
            print("{} {} {}".format(b[0][0],b[0][1],b[0][2]),file=f)
            print("{} {} {}".format(c[0][0],c[0][1],c[0][2]),file=f)

def dftb(geom,job,activeatoms,method,parapath,dispersion,kpts,hcorr):
    elements = ""
    coord_type = ""
    if '.gen' in geom:
        with open(geom,'r') as f:
            arr = next(f).split()
            coord_type+=arr[1]
            elements+=next(f)
        elements = elements.split()
    if '.xyz' in geom:
        sym = []
        elements = []
        with open(geom,'r') as f:
            next(f)
            next(f)
            for line in f:
                arr = line.split()
                sym.append(arr[0])
        for i in set(sym):
            elements.append(i)
    with open ("cmmd.in",'w') as f:
        if '.gen' in geom:
            print("""Geometry = GenFormat {{
    <<< {}
}}""".format(geom),file=f)
        if '.xyz' in geom:
            print("""Geometry = xyzFormat {{
    <<< {}
    }}""".format(geom),file=f)
        if job == 'opt':
            print("""Driver = ConjugateGradient {{
    MovedAtoms = {}
    MaxForceComponent = 1e-4
    MaxSteps = 1000
    OutputPrefix = "cmmd"
                  """.format(activeatoms),file=f)
            print('}',file=f)
        if job == 'optcell':
            print("""Driver = ConjugateGradient {{
    MovedAtoms = {}
    MaxForceComponent = 1e-4
    MaxSteps = 1000
    LatticeOpt = Yes
    OutputPrefix = "cmmd"
                  """.format(activeatoms),file=f)
            print('}',file=f)
        if method == 'XTB1':
        	kpts = kpts.split('x')
        	shift = 0
        	if int(kpts[0])%2 == 0:
        		shift+=0.5
        	else:
        		shift +=0 
        	print("""Hamiltonian = xTB {{
   Method = "GFN1-xTB"
   kPointsAndWeights = SuperCellFolding{{
      {} 0 0
      0 {} 0
      0 0 {}
      {} {} {}
   }}
        		}}""".format(kpts[0],kpts[1],kpts[2],shift,shift,shift),file=f)
        else:
        	print("Hamiltonian = DFTB {",file=f)
        	if method == 'DFTB2':
        		print("scc = Yes",file=f)
        		print("MaxSCCIterations = 1000",file=f)
        		if method == 'DFTB3diag':
        			print("""scc = Yes
    ThirdOrder = Yes
    MaxSCCIterations = 1000
""",file=f)
        	if method == 'DFTB3':
        		print("""scc = Yes
    ThirdOrderFull = Yes
              """,file=f)
        ## Koreksi ikatan hidrogen
        hdamp = {
           'DFTB3': '4.0',
           'DFTB3diag': '4.95',
           'DFTB2': '4.5'
        }
        if hcorr == 'hdamp':
        	print("""HCorrection = Damping {{
 Exponent = {}
}}""".format(hdamp[method]),file=f)
        	if hcorr == 'H5':
        		print("""HCorrection = H5{ }""",file=f)
        	if dispersion == 'D3':
        		print("""Dispersion = DftD3 {
    Damping = ZeroDamping {
        sr6 = 0.7461
        alpha6 = 14.0
    }
    s6 = 1.0
    s8 = 3.209
}""",file=f)
        	if dispersion == 'D3BJ':
        		print("""Dispersion = DftD3 {
    Damping = BeckeJohnson {
        a1 = 0.5719
        a2 = 3.6017
    }
    s6 = 1.0
    s8 = 0.5883
        }""",file=f)
        	if dispersion == 'D3H5':
        		print("""Dispersion = DftD3{
    Damping = ZeroDamping{
     sr6 = 1.25
     alpha6 = 29.61
     }
     s6 = 1.0
     s8 = 0.49
     HHRepulsion = Yes
}""",file=f)
        	print("""SlaterKosterFiles = Type2FileNames {{
        Prefix = {}/
        Separator = "-"
        Suffix = ".skf" """.format(parapath),file=f)
        	print('}',file=f)
    # Mapping bilangan kuantum azimuth ke penamaan orbital
        	azi2orb = {'1':'s','2':'p','3':'d','4':'f'}

        	print("MaxAngularMomentum {",file=f)
        	for element in elements:
        		print("""{} = "{}" """.format(element,azi2orb[azimuth(element)]),file=f)
        	print("}",file=f)

        	if method == 'DFTB3':
        		print("HubbardDerivs {",file=f)
        		for element in elements:
        			print("""{} = {}""".format(element,hubbard(element)),file=f)
        		print("}",file=f)
    # Informasi K-points
        	kpts = kpts.split('x')
        	shift = 0
        	if int(kpts[0])%2 == 0:
        		shift+=0.5
        	else:
        		shift +=0
        # if 'F' or 'S' in coord_type:
        	if '.gen' in geom:
        		print("""KPointsAndWeights = SuperCellFolding {{
                	{} 0 0
                	0 {} 0
                	0 0 {}
                	{} {} {}
            	}}""".format(kpts[0],kpts[1],kpts[2],shift,shift,shift),file=f)

        	print("}", file=f)
        if job == 'dos':
            print("Analysis {",file=f)
            print("   ProjectStates {",file=f)
            for element in elements:
                print("""Region {{
        Atoms = {}
        ShellResolved = Yes
        Label = "dos_{}"
                    }}""".format(element,element),file=f)
            print("    }",file=f)
            print("}",file=f)

        print("""ParserOptions {
  ParserVersion = 7
}""",file=f)
