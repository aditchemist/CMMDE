#!/usr/bin/env python3

import sys
import numpy as np
from pymatgen.core import lattice, structure 
import statsmodels.api as sm
import matplotlib
matplotlib.use('PS')
#matplotlib.rc('text', usetex=True)
#matplotlib.use('Qt5agg')
import matplotlib.pyplot as plt

def msd_com(groups,traject,start,latt,dt):
    index_groups = []

    entries = groups.split(':')
    for entry in entries:
        arr = entry.split('*')
        if (len(arr) != 3):
            print('Salah memasukkan format definisi molekul, pastikan dengan format     Nama*JumlahAtom*JumlahMolekul')
            parser.print_help()
            sys.exit(1)
        index_groups.append( (arr[0], int(arr[1]), int(arr[2])) )
    # print(index_groups)


    traj_file = traject
    print('Titik acuan: ', start)

    symbols = []
    symbol_list = []
    nat_list = []
    cell = []

    with open(latt, 'r') as f:
        for line in f:
            if 'TV' in line:
                cell.append(list(map(float, line.split()[1:4])))

    box = lattice.Lattice(cell)

    #membaca geometri awal dari file trayektori dinamika molekul
    with open(traj_file, 'r') as f:
        nat = int(f.readline())
        info = f.readline()
        arr = info.split()
        # t1 = float(arr[3])
        for i in range(0, nat):
            line = f.readline()
            arr = line.split()
            symbols.append(arr[0])
            if (len(symbol_list) == 0):
                symbol_list.append(arr[0])
                nat_list.append(1)
            else:
                if (arr[0] == symbol_list[-1]):
                    nat_list[-1] += 1
                else:
                    symbol_list.append(arr[0])
                    nat_list.append(1)
        nat = int(f.readline())
        info = f.readline()
        arr = info.split()
        # t2 = float(arr[3])
        # dt = t2-t1

    groups_nat = sum([k[1]*k[2] for k in index_groups])
    if (groups_nat != nat):
        raise ValueError('Salah dalam mendefinisikan molekul! Jumlah atom tidak konsisten.')

    print('Kelompok:')
    for k, nat, nmole in index_groups:
        print ('    Nama: {}, Jumlah atom dalam molekul: {}, Jumlah molekul dalam sistem: {}'.  format(k, nat, nmole))
    print("")


    first = True
    first_coords = []


    def load_xyz(lattice, f):
        nat = int(next(f))
        title = next(f)
        species = []
        coords = []
        for i in range(nat):
            line = next(f)
            arr = line.split()
            species.append(arr[0])
            coords.append(list(map(float, arr[1:4])))
        mole = structure.Structure(lattice=lattice, species=species, coords=coords,     coords_are_cartesian=True)

        return mole


    fout = open('msd.out', 'w')
    format_str = '{:>15s}(A^2)'*len(index_groups)

    print ('#Time(ps) ', format_str.format(*[k for k, nat, nmole in index_groups]), file=fout)

    first_coms = None
    last_coms = None
    image_lasts = None

    times = []
    msds = []

    with open(traj_file, 'r') as f:
        step = 0
        while(True):
            try:
                step += 1

                struct = load_xyz(box, f)
                atom_index = 0
                group_index = 0
                com_step = []
                for group in index_groups:
                    # print('Menghitung pusat massa untuk molekul ke-{}'.format(group))
                    for mole_index in range(group[2]):
                        # print('Menghitung pusat massa untuk molekul ke-{}'.format (mole_index))
                        local_indexs = [atom_index + x for x in range(group[1])]
                        # print(local_indexs)
                        if (len(local_indexs) == 1):
                            com_step.append(struct.sites[local_indexs[0]].frac_coords)
                        else:
                            atom1 = struct[local_indexs[0]]
                            com_local = atom1.specie.data['Atomic mass']*atom1.frac_coords
                            total_mass = atom1.specie.data['Atomic mass']
                            for at2 in local_indexs[1:]:
                                atom2 = struct[at2]
                                dis_, image_ = atom1.distance_and_image(atom2)
                                com_local += atom2.specie.data['Atomic mass']*(atom2.   frac_coords+image_)
                                total_mass += atom2.specie.data['Atomic mass']
                            com_local /= total_mass
                            com_step.append(com_local)
                        atom_index = local_indexs[-1]+1

                if (step == 1):
                    first_coms = com_step.copy()
                    last_coms = com_step.copy()
                    image_lasts = [(0,0,0)]*len(last_coms)
                else:
                    msd = []
                    image_cur_list = []

                    for first_com, last_com, cur_com, image_last in zip(first_coms,     last_coms, com_step, image_lasts):
                        image_current = struct.lattice.get_distance_and_image(last_com,     cur_com)[1]
                        image_final = image_last + image_current
                        final_dis = struct.lattice.get_distance_and_image(first_com, cur_com,   jimage=image_final)[0]
                        msd.append(final_dis*final_dis)
                        image_cur_list.append(image_final)

                    image_lasts = image_cur_list
                    last_coms = com_step.copy()

                    msd_group = []
                    for i in range(len(index_groups)):
                        items = msd[0:index_groups[i][2]]
                        msd = msd[index_groups[i][2]:]
                        msd_group.append(sum(items)/len(items))
                    format_str = '{:<10.5f} ' + '{:20.12f}'*len(index_groups)
                    print (format_str.format(step*dt/1000.0, *msd_group), file=fout)

                print('Tahap {} selesai'.format(step))

            except StopIteration as e:
                break
def msd_fit(msd,noheader,start):    
    data_file = msd
    x = []
    ys = []
    ncol = -1
    lineno = 0
    headers = []
    with open(data_file, 'r') as f:
        for line in f:
            lineno +=1
            if (line.startswith('#')):
                if noheader:
                    continue
                if (lineno == 1):
                    headers = line.split()[1:]
                continue
            arr = line.split()
            x.append(float(arr[0]))
            ys.append(list(map(float, arr[1:])))

            if (ncol == -1):
                ncol = len(arr) -1
            else:
                if (ncol != (len(arr) -1)):
                    raise ValueError('The number of columns of data file is not consistent!')


    print('Terdapat {} kolom dalam data MSD yang diberikan.'.format(ncol))
    print()

    print(headers)
    if (len(headers) != ncol):
        headers = []
        for col in range(ncol):
            headers.append('Kolom {}'.format(col+1))

    for col in range(ncol):
        print ('Kolom {}:'.format(col+1))

        all_data_x = np.array(x, dtype=np.double).reshape(-1, 1)
        data_x = all_data_x[start:-1]

        all_data_y = [d[col] for d in ys]
        data_y = np.array(all_data_y[start:-1], dtype=np.double).reshape(-1, 1)
        X = sm.add_constant(data_x)
        ols = sm.OLS(data_y, X)
        ols_result = ols.fit()

        summ = ols_result.summary()
        # print(summ)
        print('*'*80)
        diff_coeff = ols_result.params[1]/6.0
        print(' Koefisien Difusi = {:.2f} A^2/ps = {:.2e} m^2/s'.format  (diff_coeff, diff_coeff*1.0e-8))
        # diff_coeff = ols_result.params[1]/2.0
        # print(' Diffusion coefficient (1/2): {:16.10F} A^2/ps, {:17.12E} m^2/s'.format  (diff_coeff, diff_coeff*1.0e-8))
        # diff_coeff = ols_result.params[1]/4.0
        # print(' Diffusion coefficient (1/4): {:16.10F} A^2/ps, {:17.12E} m^2/s'.format  (diff_coeff, diff_coeff*1.0e-8))
        print('*'*80)
        print()

        predicted_y = ols_result.predict(X)
        
            
        plt.plot(x, all_data_y, label=headers[col].strip('(A^2)'))
        plt.plot(data_x, predicted_y, 'k--')
    
    plt.tick_params(direction='in')
    plt.legend(loc=2,frameon=False)
    plt.grid(linestyle=':',linewidth=1)
    plt.xlabel('Time [ps]')
    plt.ylabel('MSD [$\AA^2$]')

    plt.savefig('msd.pdf', dpi=1200, format='pdf')
    # plt.show()
