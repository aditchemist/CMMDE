def gen2poscar(input):
    fname = input.split(".")
    type = ""
    Natom = 0
    elem_type = []
    coord = []
    index = []
    cell = []
    N_elem = []
    with open(input, 'r') as f:
        header = next(f).split()
        Natom += int(header[0])
        type += header[1]
        header2 = next(f).split()
        for i in header2:
            elem_type.append(i)
        for i in range(Natom):
            arr = next(f).split()
            coord.append(arr[2:])
            index.append(arr[1])
        next(f)
        for line in f:
            arr = line.split()
            cell.append(arr[0:])
    for i in set(index):
        Nele = index.count(i)
        N_elem.append(Nele)
    sym = []
    for k in range(len(N_elem)):
        for j in elem_type:
            for i in range(int(N_elem[k])):
                sym.append(j)
    with open('{}.vasp'.format(fname[0]),'w') as f:
        print("CMMDE poscar file",file=f)
        print("1.0",file=f)
        if cell[0] !=0:
            for i in cell:
                print("{} {} {}".format(i[0],i[1],i[2]),file=f)
        for i in elem_type:
            print("{}".format(i),end=" ", file=f)
        print("",file=f)
        for i in N_elem:
            print("{}".format(i),end=" ", file=f)
        print("",file=f)
        print("Cart",file=f)
        for coord, sym in zip(coord,sym):
            print("{:16.12f} {:16.12f} {:16.12f} {} ".format(float(coord[0]),float(coord[1]),float(coord[2]),sym),file=f)
