def nwchem(job,method,nproc,geom,charge,mult,scalfreq,restart,conv_thr,iter):
    # Extracting geometry information
    symb = []
    x = []
    y = []
    z = []
    with open(geom, 'r') as f:
        nAtom = next(f)
        next(f)
        for line in f:
            arr = line.split()
            symb.append(arr[0])
            x.append(arr[1])
            y.append(arr[2])
            z.append(arr[3])
    with open('cmmd.in','w') as f:
        if restart == 'true':
            print("restart cmmde_nw",file=f)
        else:
            print("start cmmde_nw",file=f)
        print("""title "NWChem calculation"

geometry units au""",file=f)
        for i,x,y,z in zip(symb,x,y,z):
            print("{} {} {} {}".format(i,x,y,z),file=f)
        print("end",file=f)
        basis_set = method.split()[1]
        print("""basis
* library {}
end""".format(basis_set),file=f)
        print("""scf
thresh {}
maxiter {}
end""".format(conv_thr,iter),file=f)
        if job == 'opt':
            print("task scf optimize",file=f)
        else:
            print("task scf",file=f)
