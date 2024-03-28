def poscar2gen(input):
    fname = input.split(".")
    vx = []
    vy = []
    vz = []
    with open(input, 'r') as f:
        next(f)
        next(f)
        for i in range(3):
            arr = next(f).split()
            vx.append(arr[0])
            vy.append(arr[1])
            vz.append(arr[2])
        elem_type = next(f).split()
        elem_num = next(f).split()
