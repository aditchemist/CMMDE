def orca(job,method,nproc,geom,charge,mult,scalfreq,temperature,pressure,nroots,tda,solvent,constraints,qmatoms,totalcharge,totalmult,qm2method,qm2basis,activeatoms,hessfile,disp,aim,product,ts,irciter, printlevel, inithess,grid, finalgrid,maxiter):
    if 'opt' in job or 'Opt' in job or 'OPT' in job and not 'ts' in job: 
        method += ' opt'
    if 'freq' in job and ('XTB2' in method or 'XTB' in method):
        method += ' Numfreq'
    if 'irc' in job:
        method += ' IRC'
    if 'ts' in job and 'neb' not in job:
        method += ' OptTS'
    if 'neb-ts' in job:
        method += ' NEB-TS'
    if 'neb' in job and 'ts' not in job:
        method += ' NEB'
    if 'neb-ci' in job:
        method += ' NEB-CI'
    if 'freq' in job and ('XTB2' not in method and 'XTB' not in method and 'XTBFF' not in method):
        method += ' freq'
    if disp != 'None':
        method += ' {}'.format(disp)
    if aim == 'true' or aim == 'True' or aim == 'TRUE':
        method += ' AIM'
    if solvent != 'none' and ('XTB2' in method or 'XTB' in method or 'XTBFF' in method):
        method += ' ALPB({})'.format(solvent)
    if solvent != 'none' and ('XTB2' not in method and 'XTB' not in method and 'XTBFF' not in method):
        method += ' CPCM({})'.format(solvent)
    
    with open('cmmd.in','w') as f:
        print("#CMMDE generated Orca input file", file=f)
        print("!{}".format(method), file=f)
        print("""%pal 
 nprocs {} 
end""".format(nproc),file=f)
        if 'opt' in job:
        	print("""%geom
  maxiter {}
end""".format(maxiter),file=f)
        if qmatoms != 'None':
            print("""%qmmm
   QMAtoms {{ {} }} end
   Charge_Total {}
   Mult_Total {}""".format(qmatoms,totalcharge,totalmult),file=f)
            if qm2method != 'None':
                print("""   QM2CUSTOMMETHOD "{}"
   QM2CUSTOMBASIS "{}" """.format(qm2method,qm2basis),file=f)
            if activeatoms != 'None':
                print("ActiveAtoms {} end".format(activeatoms))
            print("end",file=f)
        if constraints != 'None':
            print("%geom", file=f)
            print("Constraints",file=f)
            for i in constraints.split(','):
                print("""{{ {} }}""".format(i),file=f)
            print("end",file=f)
            print("end",file=f)
        if hessfile != 'None':
            print("""%geom
    MaxIter {}
    InHess Read
    InHessName "{}"
end""".format(maxiter,hessfile),file=f)
        # Baris koordinat atom
        print("""     
*xyzfile {} {} {}
""".format(charge,mult,geom)
        , file=f)
        if 'freq' in job:
            print("""%freq 
 scalfreq {} 
 Temp {}
 Pressure {}
end""".format(scalfreq,temperature,pressure),file=f)
        if 'td' in job:
            print("""%tddft 
  nroots {}
  tda {}
end""".format(nroots,tda),file=f)
        if 'neb' in job:
        	print("""%geom 
  maxiter {}
end""".format(maxiter),file=f)
        	print("""%NEB
  neb_end_xyzfile "{}" """.format(product),file=f)
        	if ts != 'None':
        		print("""  neb_ts_xyzfile "{}"
end""".format(ts),file=f)
        	else:
        		print("end",file=f)     
        if 'irc' in job:
            print("""%irc
  MaxIter {}
  PrintLevel {}
  InitHess {}
  Hess_filename "{}"
end

%method
  Grid {}
  FinalGrid {}
end
""".format(irciter, printlevel, inithess, hessfile, grid, finalgrid),file=f)
    return

