import os
from cmmde_hubbard import hubbard, azimuth
def dcdftb(job,method,geom,charge,mult,dispersion,para_path,temp,pressure,ensembel,thermostat,deltat,step,mdprint,a1,a2,a3,b1,b2,b3,c1,c2,c3,restart,traject,velocity,dftbinp,softpot,softtype,softrange,softcenter,metarest,metafreq,metaheight,cvtype,metawidth,pow1,pow2,rcut,fesstart,fesend,fesbin,ag1,ag2,ag3,ag4,solvent,nroots,target,transmult,ocsstr,writetrans,lc,bufrad,delta,opttype,freqtype,econv,dconv):
    symb = []
    x = []
    y = []
    z = []
    with open(geom, 'r') as f:
        Natom = int(next(f).strip())
        next(f)
        for line in f:
            arr = line.split()
            symb.append(arr[0])
            x.append(arr[1])
            y.append(arr[2])
            z.append(arr[3])
    # Menghitung jumlah atom
    
    elements = list(set(symb))
    Numelements = len(elements)
    if 'TV' in elements:
        elements.remove('TV')
        Numelements = len(elements)
    with open("cmmd.in",'w') as f:
        if method == "DFTB2_gammah":
            print("DC=FALSE",file=f)
            print("SCC=(LC={} DAMPXH=TRUE DAMPXHZETA=4.5 ECONV={} DCONV={})".format(lc,econv,dconv), file=f)
        elif method == "DFTB2":
            print("DC=FALSE",file=f)
            print("SCC=(LC={} ECONV={} DCONV={})".format(lc, econv, dconv), file=f)
        elif method == "DFTB3":
            print("DC=FALSE",file=f)
            if dispersion != 'D3H5':
                print("SCC=(THIRDFULL=TRUE DAMPXH=TRUE DAMPXHZETA=4.0 LC={} ECONV={} DCONV={})".format(lc, econv, dconv),file=f)
            else:
                print("SCC=(THIRDFULL=TRUE DAMPXH=FALSE LC={} ECONV={} DCONV={})".format(lc, econv, dconv),file=f)
            
        elif method == "DFTB3-diag":
            print("DC=FALSE",file=f)
            print("SCC=(THIRDDIAG=TRUE DAMPXH=TRUE DAMPXHZETA=4.95 LC={} ECONV={} DCONV={})".format(lc, econv, dconv),file=f)
        elif method == "DCDFTB2_gammah":
            print("DC=(BUFRAD={} DELTARXYZ={})".format(bufrad,delta),file=f)
            print("SCC=(LC={} DAMPXH=TRUE DAMPXHZETA=4.5 ECONV={} DCONV={})".format(lc,econv,dconv), file=f)
        elif method == "DCDFTB2":
            print("DC=(BUFRAD={} DELTARXYZ={})".format(bufrad,delta),file=f)
            print("SCC=(LC={} ECONV={} DCONV={})".foramt(lc, econv, dconv),file=f)
        elif method == "DCDFTB3":
            print("DC=(BUFRAD={} DELTARXYZ={})".format(bufrad,delta),file=f)
            print("SCC=(THIRDFULL=TRUE DAMPXH=TRUE DAMPXHZETA=4.0 LC={} ECONV={} DCONV={})".format(lc, econv, dconv),file=f)
        elif method == "DCDFTB3-diag":
            print("DC=(BUFRAD={} DELTARXYZ={}))".format(bufrad,delta),file=f)
            print("SCC=(THIRDDIAG=TRUE DAMPXH=TRUE DAMPXHZETA=4.95 LC={} ECONV={} DCONV={})".format(lc, econv, dconv),file=f)
        # Mapping pelarut implisit
        gbsolventtype = { 'water': '1','acetonitril': '2','dimethylsulfoxide':'3','methanol':'4','chloroform':'5','dichloromethane':'6','benzene':'7','toluene':'8','acetone':'9','tetrahydrofuran':'10','diethylether':'11','carbondisulfide':'12','dimethylformamide':'13','hexane':'14','nitromethane':'15','furan':'16','dioxane':'17','ethylacetate':'18','phenol':'19','aniline':'20','benzaldehyde':'21','octanol':'23','dimethylacetamide':'26','ethanol':'25','formamide':'27','methylethylketone':'28'}
        if solvent != 'none':
            if job == 'td':
                if solvent == 'ethanol': gbsolventtype[solvent] = 4
                print("SCC=(ALPB=TRUE GBSOLVENTTYPE={} GBPARAMTYPE=7)".format(gbsolventtype[solvent]),file=f)
                print("SCC=(SA=TRUE GBSAHBOND=TRUE)",file=f) 
            else:   
                print("SCC=(ALPB=TRUE GBSOLVENTTYPE={} GBPARAMTYPE=4)".format(gbsolventtype[solvent]),file=f)
                print("SCC=(SA=TRUE GBSAHBOND=TRUE)",file=f)
        if job == 'td':
            print("TD = (NSTATE={} TARGETSTATE={} MULT={} OSCSTR={} WRITETRANSITION = {})".format(nroots,target,transmult,ocsstr,writetrans),file=f)
        if dispersion == "D3":
            print("DISPERSION=(DISPTYPE=4)",file=f)
        elif dispersion == "D2":
            print("DISPERSION=(DISPTYPE=1)",file=f)
        elif dispersion == "SK":
            print("DISPERSION=(DISPTYPE=2)",file=f)
        elif dispersion == "LJ":
            print("DISPERSION=(DISPTYPE=3)",file=f)
        elif dispersion == "D3":
            print("DISPERSION=(DISPTYPE=4)",file=f)
        elif dispersion == "D3BJ":
            print("DISPERSION=(DISPTYPE=5)",file=f)
        elif dispersion == "D3H4":
            print("DISPERSION=(DISPTYPE=6)",file=f)
        elif dispersion == "DFTulg":
            print("DISPERSION=(DISPTYPE=7)",file=f)
        elif dispersion == "dDMC":
            print("DISPERSION=(DISPTYPE=8)",file=f)    
        elif dispersion == "D3H5":
            print("DISPERSION=(DISPTYPE=9)",file=f)
        elif dispersion == "kubilius":
            print("DISPERSION=(DISPTYPE=10)",file=f)
        elif dispersion == "goursot":
            print("DISPERSION=(DISPTYPE=11)",file=f)
        elif dispersion == "vdwts":
            print("DISPERSION=(DISPTYPE=12)",file=f)
        elif dispersion == "manybody":
            print("DISPERSION=(DISPTYPE=13)",file=f)
        elif dispersion == "D4":
            print("DISPERSION=(DISPTYPE=14)",file=f)

        if "opt" in job:
            opt_type={'bfgs':'1','sd':'2','cg':'3','qm':'4','fire':'5'}
            print("OPT=(MAXITER=1000 OPTTYPE={})".format(opt_type[opttype]),file=f)
        if "cell" in job:
            opt_type={'bfgs':'1','sd':'2','cg':'3','qm':'4','fire':'5'}
            print("OPT=(MAXITER=1000 OPTTYPE={} LATTICEOPT=TRUE)".format(opt_type[opttype]),file=f)
        if "freq" in job:
            print("FREQ=(THERMOTEMP={} THERMOPRES={} FREQTYPE={})".format(temp,pressure,freqtype),file=f)
        if ensembel == 'NVE':
            nvt = 'FALSE'
        elif ensembel == 'NVT':
            nvt = 'TRUE'
        if thermostat == 'nose':
            thermo = 3
        elif thermostat == 'berendsen':
            thermo = 4
        elif thermostat == 'andersen':
            thermo = 5
        if job == "md":
            if ensembel == 'NVT':
                nvt = 'TRUE'
            else:
                nvt = 'FALSE'
            print('MD=(NVT={} NVTTYPE={} BATHTEMP={} ERRORTEMP=100)'.format(nvt,thermo,temp),file=f)
            print('MD=(DELTAT={} NSTEP={} PRINT={})'.format(deltat*1e-15,step,mdprint),file=f)
            if softpot == 'true':
                print('MD=(SOFT=TRUE SOFTSHAPETYPE={} SOFTRANGE={} SOFTCENTERTYPE={})'.format(softtype,softrange,softcenter),file=f)
            if restart == 'true':
                print('MD=(READVELOCITY=TRUE)',file=f)
        if job == "mtd":
            if ensembel == 'NVT':
                nvt = 'TRUE'
            else:
                nvt = 'FALSE'
            if metarest == 'false':
                metarest = 'FALSE'
            else:
                metarest = 'TRUE'                
            print('MD=(NVT={} NVTTYPE={} BATHTEMP={} ERRORTEMP=100)'.format(nvt,thermo,temp),file=f)
            print('MD=(DELTAT={} NSTEP={} PRINT={})'.format(deltat*1e-15,step,mdprint),file=f)
            print('MD=(METADYNAMICS=TRUE METARESTART={} METAPRINTFES=TRUE METAFREQ={} METAHEIGHT={})'.format(metarest,metafreq,metaheight),file=f)
            if softpot == 'true':
                print('MD=(SOFT=TRUE SOFTSHAPETYPE={} SOFTRANGE={} SOFTCENTERTYPE={})'.format(softtype,softrange,softcenter),file=f)
            if restart == 'true':
                print('MD=(READVELOCITY=TRUE)',file=f)    
            
            with open('metacv.dat','w') as fcv:
                atomgroup1 = ag1.split()
                atomgroup2 = ag2.split()
                nag1 = len(atomgroup1)
                nag2 = len(atomgroup2)
                if cvtype == 'coordnum':
                    print("""RATIONALCOORDINATIONNUMBER {} {} {} {} {} {} LABELS AVERAGE {} {} {} 
{}
{}
                """.format(metawidth,nag1,nag2,pow1,pow2,rcut,fesstart,fesend,fesbin,ag1,ag2),file=fcv)
                if cvtype == 'distance':
                    print("BONDDISTANCE {} {} {} {} {} {}".format(metawidth,ag1,ag2,fesstart,fesend,fesbin),file=fcv)
                if cvtype == 'angle':
                    print("BONDANGLE {} {} {} {} {} {} {}".format(metawidth,ag1,ag2,ag3,fesstart,fesend,fesbin),file=fcv)
                if cvtype == 'dihedral':
                    print("BONDDIHEDRAL {} {} {} {} {} {} {} {}".format(metawidth,ag1,ag2,ag3,ag4,fesstart,fesend,fesbin), file=fcv)
                if cvtype == 'distancediff':
                    print("BONDDISTANCEDIFFERENCE {} {} {} {} {} {} {} {}".format(metawidth,ag1,ag2,ag3,ag4,fesstart,fesend,fesbin),file=fcv)
                if cvtype == 'distanceadd':
                    print("BONDDISTANCEADDITION {} {} {} {} {} {} {} {}".format(metawidth,ag1,ag2,ag3,ag4,fesstart,fesend,fesbin),file=fcv)
                if cvtype == 'meandistance':
                    print("""MEANDISTANCE {} {} {} {} LABELS {} {} {}
{}
{}
                    """.format(metawidth,nag1,nag2,pow1,fesstart,fesend,fesbin,ag1,ag2),file=fcv)
                if cvtype == 'pointplanedistance':
                    print("ATOMPOINTPLANCEDISTANCE {} {} {} {} {}".format(metawidth,ag1,ag2,ag3,ag4),file=fcv)
        
        print("",file=f)
        print("DCDFTB input generated from CMMDE code",file=f)
        print("",file=f)
        print(Numelements,file=f)
        # Pemetaan turunan Hubbard untuk atom-atom tertentu:
        for element in elements:
            if method == 'DFTB2' or method == 'DCDFTB2' or method == 'DFTB3-diag' or method == 'DCDFTB3-diag' or method == 'DFTB2_gammah' or method == 'DCDFTB2_gammah':
                print(element,azimuth(element),file=f)
            elif method == 'DFTB3' or method == 'DCDFTB3' or method == 'DCDFTB3-diag':
                print(element,azimuth(element),hubbard(element),file=f)
            
            ind = 0
            for element2 in elements:
                os.system("cp {}/{}-{}.skf .".format(para_path,element,element2))
                if ind < len(elements):
                    print(element+'-'+element2+'.skf',end=' ',file=f)
                    ind+=1
                if ind == len(elements):
                    print(' ',file=f)
        print("",file=f)
        print("{} {} {}".format(Natom,charge,mult),file=f)
        if restart == 'false':
            for symb, x, y, z in zip(symb,x,y,z):
                print("{} {} {} {}".format(symb,x,y,z), file=f)
        else:
            with open(traject, 'r') as ftraj:
                lines = ftraj.readlines()
                coords = []
                for i in range(-Natom,0):
                    coords.append(lines[i].strip('\n'))
            with open(velocity,'r') as fvel:
                lines = fvel.readlines()
                velocities = []
                for i in range(-Natom,0):
                    velocities.append(lines[i].strip('\n'))
            with open('veloc.dat','w') as fveloc:
                for i in velocities:
                    print(i[3:], file=fveloc)

            with open(dftbinp, 'r') as inp:
                latt_info = []
                for line in inp:
                    if 'TV' in line:
                        arr = line.strip('\n')
                        latt_info.append(arr)

            for i in coords:
                print(i, file=f)
            for i in latt_info:
                print(i, file=f)
                            
        if a1 != 0:
            print("""TV {} {} {}
TV {} {} {}
TV {} {} {}""".format(a1,a2,a3,b1,b2,b3,c1,c2,c3),file=f)
        print(" ",file=f)

            

