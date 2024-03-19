#!/usr/bin/env python3
import os
import sys
import argparse
import re
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from prettytable import PrettyTable
matplotlib.use('PS')
#matplotlib.rc('text', usetex=True)
#matplotlib.use('Qt5agg')
from cmmde_msd_com import msd_com, msd_fit
from cmmde_rdf import rdf
from cmmde_dock import checkopt

parser = argparse.ArgumentParser(description='Program Analisis Hasil Perhitungan Dalam MOWS CMMD 2021')
parser.add_argument('-j','--job',type=str,default='None',help='Jenis perhitungan yang dilakukan. Pilihan: sp, opt, thermo, ir, reax, msd_com, msd, msd_fit rdf, td.')
parser.add_argument('-m', '--method', type=str, default='XTB2', help='Metode yang digunakan dalam perhitungan.')
parser.add_argument('-i','--input',type=str,default='geom.xyz',help='Koordinat awal dalam format .xyz (hanya diperlukan untuk program DCDFTBMD) saat optimasi geometri.')
parser.add_argument('-s', '--software', type=str, default='orca',help='Jenis software yang digunakan. Pilihan: orca, dcdftb, dftb+.')
parser.add_argument('-irx','--inputreaction',type=str, default='None', help='Input reaksi kimia')
parser.add_argument('-t','--temp',type=float,default=298.15, help='Suhu yang digunakan dalam perhitungan frekuensi dan termokimia')
parser.add_argument('-traj','--traject',type=str,default='traject.xyz',help='File trayektori dinamika molekul dalam format .xyz')
parser.add_argument('-l', '--latt', type=str, default='cmmd.in', help='Ukuran rusuk sel satuan.')
parser.add_argument('-start', '--start', type=int, default=0, help='Titik acuan dalam bilangan bulat. Default = 0.')
parser.add_argument('-end','--end', type=int, default=-1, help='Titik ahir pembacaan trayektori')
parser.add_argument('-msd','--msd',type=str, default='msd.out', help='File berisikan MSD keluaran cmmdepost.py -j msdcom.')
parser.add_argument('-dt','--dt',type=float,default=0.5, help='Time step yang didefinisikan selama simulasi berlangsung.')
parser.add_argument('-n', '--noheader', action="store_true",help='Mendefinisikan apakah file mengandung header atau tidak. Default file msd.out mengandung header berupa #.')
parser.add_argument('-g', '--groups', type=str, help='Molekul didefinisikan dengan Nama*JumlahAtom*JumlahMolekul. Contoh NAP*1*2:DMF*1*20.')
parser.add_argument('-r', '--range', default=10.0, type=float,help='Jarak terjauh pengukuran RDF')
parser.add_argument('-p', '--pair', type=str, help='Pasangan atom, misalnya "O-O" untuk pasangan atom O-O')
parser.add_argument('-res','--resolution',default=0.01, type=float, help='Resolusi untuk plot histogram atau RDF')
# Metadinamika
parser.add_argument('-metafreq','--metafreq',type=int,default=100, help='Dalam berapa step sekali potensial Gaussian ditambahkan?')
parser.add_argument('-fesstart','--fesstart',type=float,default=0, help='Titik minimum CV. Default = 0.')
parser.add_argument('-fesend','--fesend',type=float,default=1, help='Titik maksimum CV. Default = 1.')
parser.add_argument('-fesbin','--fesbin',type=float,default=0.01,help='Selang CV. Default = 0.01')
parser.add_argument('-cv','--cvtype',type=str,default='coordnum',help='Pilihan collective variable (CV) yang digunakan. Pilihan: coordnum, distance, angle, dihedral, distancediff, distanceadd, meandistance, pointplanedistance.')
parser.add_argument('-ngaus','--ngaus',type=int,help='Jumlah Gaussian yang didepositkan')
parser.add_argument('-nlig','--nligands',type=int,help='Jumlah ligan yang akan didocking.')
parser.add_argument('-nr','--nroots',type=int,default=5,help='Jumlah keadaan tereksitasi.')
parser.add_argument('-pts','--points',type=int,default=300,help='Jumlah titik untuk plot UV.')
# Analisis DOS
parser.add_argument('-pdos','--pdos',type=str,help='File DOS parsial hasil perhitungan DFTB+')
parser.add_argument('-xlabel','--xlabel',type=str,default='x',help='Label sumbu-x dalam sebuah plot')
parser.add_argument('-ylabel','--ylabel',type=str,default='y',help='Label sumbu-y dalam sebuah plot')
parser.add_argument('-file','--file',type=str,help='Nama file yang akan diplot')
parser.add_argument('-np','--nproc',type=int,default=1,help='Jumlah prosesor yang digunakan.')
parser.add_argument('-grid','--grid',type=str,default='FINE',help='Jenis grid yang akan dibangun untuk perhitungan NCI. Pilihan: FINE (default), COARSE, dan ULTRAFINE.')
# Analisis terkait perhitungan Quantum Espresso
parser.add_argument('-outdir','--outdir',type=str,default='./out',help='Folder tempat menyimpan output.')
parser.add_argument('-emin','--emin',type=float,default=-9.0,help='Batas minimum nilai keadaan elektronik pada DOS.')
parser.add_argument('-emax','--emax',type=float,default=16.0,help='Batas maksimum nilai keadaan elektronik pada DOS.')

opt = parser.parse_args(sys.argv[1:])

# Berbagai faktor konversi
Ha2kcal = 627.5096080305927
ev2kcal = 23.060541945329334
Ha2kj = 2625.5
ev2kj = 96
# Analisis Output
method = opt.method
if opt.job == 'sp' and opt.software == 'orca' and not (opt.method == 'XTB' or opt.method == 'XTB2'):
    if os.path.isfile('cmmd.out'):
        with open('cmmd.out','r') as f:
            for line in f:
                if "Total Energy" in line:
                    arr= line.split()
                    Energy = float(arr[3])

        print('Energi total molekul: {} Hartree = {:.2f} kJ/mol'.format(Energy,Energy*Ha2kj))
    else:
        print("Mohon menunggu hingga perhitungan anda selesai.")
        exit

elif opt.job == 'sp' and opt.software == 'orca':
    if os.path.isfile('cmmd.out'):
        with open('cmmd.out','r') as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    arr= line.split()
                    Energy = float(arr[3])
                if "GRADIENT NORM" in line:
                    arr = line.split()
                    Gradient = float(arr[3])
                if "HOMO-LUMO GAP" in line:
                    arr = line.split()
                    Egap = float(arr[3])

        print('Energi total molekul: {} Hartree = {:.2f} kJ/mol'.format(Energy,Energy*Ha2kj))
        print('Norm Gradien: {} Hartree/a0 = {:.2f} kJ/mol/a0'.format(Gradient,Gradient*Ha2kj))
        print('HOMO-LUMO gap: {:.2f} eV = {:.2f} kJ/mol'.format(Egap,Egap*ev2kj))
    else:
        print("Mohon menunggu hingga perhitungan anda selesai.")
        exit

if opt.job == 'opt' and opt.software == 'orca':
    Energy = []
    Gradient = []
    TotalEnergy = []
    if os.path.isfile('cmmd.out'):
        with open('cmmd.out','r') as f:
            for line in f:
                if "Current Energy" in line:
                    arr = line.split()
                    Energy.append(float(arr[3]))
                if "Current gradient" in line:
                    arr = line.split()
                    Gradient.append(float(arr[4]))
                if "FINAL SINGLE POINT ENERGY" in line:
                    arr = line.split()
                    TotalEnergy.append(float(arr[4]))
        Optstep = range(0,len(Energy))
        print("Total energi elektronik = {} Hartree = {} kJ/mol".format(TotalEnergy[-1],TotalEnergy[-1]*Ha2kj))
        if method != "XTB1" and method != "xtb1" and method != "XTB2" and method != "xtb2" and method != "XTB" and method != "xtb":
            with open('cmmd.out','r') as f:
                lines = f.readlines()
                occ = []
                energy = []
                for index,line in enumerate(lines):
                    if 'Basis Dimension' in line:
                        arr = line.split()
                        NBas = int(arr[4])
                    if 'ORBITAL ENERGIES' in line:
                        for i in range(NBas):
                            arr = lines[index+i+4].split()
                            occ.append(arr[1])
                            energy.append(float(arr[3]))
                Eocc = []
                Evir = []
                for energy in energy:
                    if energy > 0:
                        Evir.append(energy)
                    if energy < 0:
                        Eocc.append(energy)
                E_HOMO = max(Eocc)
                E_LUMO = min(Evir)
                print('######INFORMASI ENERGI HOMO & LUMO######')
                print('Energi HOMO = {:.2f} eV'.format(E_HOMO))
                print('Energi LUMO = {:.2f} eV'.format(E_LUMO))
                print('Gap HOMO-LUMO = {:.2f} eV'.format(E_LUMO-E_HOMO))
        else:
            Gaps = []
            with open("cmmd.out", 'r') as f:
                print('######INFORMASI ENERGI HOMO & LUMO######')
                for line in f:
                    if "GAP" in line:
                        arr = line.split()
                        Gaps.append(float(arr[3]))
            print("Gap HOMO-LUMO = {:.2f} eV".format(Gaps[-1]))
                
        with open('optimized.dat','w') as fout:
            print("#Step Energy Gradient", file=fout)
            for optstep,energy,gradient in zip(Optstep,Energy,Gradient):
                print("{} {} {}".format(optstep,energy,gradient),file=fout)

        ## Plotting
        with open('optimized.gp','w') as fout:
            print("""set terminal pdf
set output "optimized.pdf"
set title "Energy Vs. Optimization Step"
set xlabel "Optimization Step"
set xtics 1
set ylabel "Energy [Hartree]"
set style line 1 \\
    linecolor rgb '#0060ad' \\
    linetype 1 linewidth 2 \\
    pointtype 7 pointsize 1.5

plot 'optimized.dat' using 1:2 notitle with linespoints linestyle 1""",file=fout)

        os.system('gnuplot optimized.gp')
    else:
        print("Mohon menunggu hingga perhitungan anda selesai.")
        exit

if opt.job == 'ir' and opt.software == 'orca':
    def gaussband(x,band,strength,stdev):
        bandshape=1.3062974e8 * (strength / (1e7/stdev))  * np.exp(-(((1.0/x)-(1.0/band))/(1.0/stdev))**2)
        return bandshape

    Freq = []
    Int = []
    NFreq = 0
    if os.path.isfile('cmmd.out'):
        with open('cmmd.out','r') as f:
            lines = f.readlines()

            for index,line in enumerate(lines):
                if "The total number of vibrations considered is" in line:
                    arr = line.split()
                    NFreq += int(arr[7])
            for index,line in enumerate(lines):
                if "Total enthalpy" in line:
                    arr = line.split()
                    Enthalpy = float(arr[3])
                if "Final Gibbs free energy" in line:
                    arr = line.split()
                    Gibbs = float(arr[5])
                # if "Total entropy" in line:
                #     arr = line.split()
                #     Entropy = float(arr[4])
                if "IR SPECTRUM" in line:
                    for i in range(NFreq):
                        arr = lines[index+6+i].split()
                        Freq.append(float(arr[1]))
                        Int.append(float(arr[2]))
        x = np.linspace(500,4000,5000)
        composite = 0
        for count,peak in enumerate(Freq):
            thispeak = gaussband(x,peak,Int[count],1e5)
            composite+=thispeak

        fig,ax = plt.subplots()
        ax.plot(x,composite,color='blue')
        ax.tick_params(direction='in')
        plt.grid(linestyle=':')
        plt.xlim(500,4000)
        plt.ylim(0,)
        plt.xlabel('Wavenumbers [cm$^{-1}$]')
        plt.ylabel('Absorbance')
        plt.savefig('IR.pdf',dpi=1200,format='pdf')

        with open('IR_fit.dat','w') as fout:
            for x,y in zip(x,composite):
                print(x,y,file=fout)

    else:
        print("Mohon menunggu hingga perhitungan anda selesai.")
        exit

if opt.job == 'thermo' and opt.software == 'orca':

    if os.path.isfile('cmmd.out'):
        with open('cmmd.out','r') as f:
            for line in f:
                if "Total enthalpy" in line:
                    arr = line.split()
                    Enthalpy = float(arr[3])
                if "Final Gibbs free energy" in line:
                    arr = line.split()
                    Gibbs = float(arr[5])
                # if "Total entropy" in line:
                #     arr = line.split()
                #     Entropy = float(arr[4])
                if "Non-thermal (ZPE)" in line:
                    arr = line.split()
                    ZPE = float(arr[3])
                if "Electronic energy" in line:
                    arr = line.split()
                    E_elec = float(arr[3])
                if "Thermal vibrational correction" in line:
                    arr = line.split()
                    Evib = float(arr[4])
                if "Thermal rotational correction" in line:
                    arr = line.split()
                    Erot = float(arr[4])
                if "Thermal translational correction" in line:
                    arr = line.split()
                    Etrans = float(arr[4])
        if opt.method == 'XTB2' or opt.method == 'XTB':
            E_HOMO = 0
            E_LUMO = 0
            with open('cmmd.out','r') as f:
                for line in f:
                    if '(HOMO)' in line:
                        arr = line.split()
                        E_HOMO += float(arr[3])
                    if "(LUMO)" in line:
                        arr = line.split()
                        E_LUMO += float(arr[2])
            print('######INFORMASI ENERGI HOMO & LUMO######')
            print('Energi HOMO = {:.2f} eV'.format(E_HOMO))
            print('Energi LUMO = {:.2f} eV'.format(E_LUMO))
            print('Gap HOMO-LUMO = {:.2f} eV'.format(E_LUMO-E_HOMO))
            print('')

        else:
            with open('cmmd.out','r') as f:
                lines = f.readlines()
                occ = []
                energy = []
                for index,line in enumerate(lines):
                    if 'Basis Dimension' in line:
                        arr = line.split()
                        NBas = int(arr[4])
                    if 'ORBITAL ENERGIES' in line:
                        for i in range(NBas):
                            arr = lines[index+i+4].split()
                            occ.append(arr[1])
                            energy.append(float(arr[3]))
                Eocc = []
                Evir = []
                for energy in energy:
                    if energy > 0:
                        Evir.append(energy)
                    if energy < 0:
                        Eocc.append(energy)
                E_HOMO = max(Eocc)
                E_LUMO = min(Evir)
                print('######INFORMASI ENERGI HOMO & LUMO######')
                print('Energi HOMO = {:.2f} eV'.format(E_HOMO))
                print('Energi LUMO = {:.2f} eV'.format(E_LUMO))
                print('Gap HOMO-LUMO = {:.2f} eV'.format(E_LUMO-E_HOMO))
                print('')


        print('######INFORMASI ENERGI ELEKTRONIK & KOREKSI TERMAL######')
        print("Energi elektronik = {} Hartree = {} kJ/mol".format(E_elec,E_elec*Ha2kj))
        print("Zero point energy (ZPE) = {} Hartree = {} kJ/mol".format(ZPE,ZPE*Ha2kj))
        print("Koreksi termal vibrasi = {} Hartree = {} kJ/mol".format(Evib,Evib*Ha2kj))
        print("Koreksi termal rotasi = {} Hartree = {} kJ/mol".format(Erot,Erot*Ha2kj))
        print("Koreksi termal translasi = {} Hartree = {} kJ/mol".format(Etrans,Etrans*Ha2kj))
        print('Total koreksi termal = {} Hartree = {} kJ/mol'.format(Etrans+Evib+Erot,(Etrans+Evib+Erot)*Ha2kj))
        print('')
        print('######INFORMASI BESARAN TERMOKIMIA######')
        print("Entalpi (H) = {} Hartree = {} kj/mol".format(Enthalpy,Enthalpy*Ha2kj))
        # print("Entropi (S) = {} Hartree/K = {} J/(mol K)".format(Entropy/Temperature,Entropy*Ha2kj/Temperature*1000))
        print("Energi bebas Gibbs (G) = {} Hartree = {} kJ/mol".format(Gibbs,Gibbs*Ha2kj))

    else:
        print("Mohon menunggu hingga perhitungan anda selesai.")
        exit

### Opsi untuk print hanya HOMO-LUMO energi
if opt.job == 'gap' and ('XTB' not in opt.method and 'XTB2' not in opt.method) and opt.software == 'orca':
    with open('cmmd.out','r') as f:
            lines = f.readlines()
            occ = []
            energy = []
            for index,line in enumerate(lines):
                if 'Basis Dimension' in line:
                    arr = line.split()
                    NBas = int(arr[4])
                if 'ORBITAL ENERGIES' in line:
                    for i in range(NBas):
                        arr = lines[index+i+4].split()
                        occ.append(arr[1])
                        energy.append(float(arr[3]))
            Eocc = []
            Evir = []
            for energy in energy:
                if energy > 0:
                    Evir.append(energy)
                if energy < 0:
                    Eocc.append(energy)
            E_HOMO = max(Eocc)
            E_LUMO = min(Evir)
            print('######INFORMASI ENERGI HOMO & LUMO######')
            print('Energi HOMO = {:.2f} eV'.format(E_HOMO))
            print('Energi LUMO = {:.2f} eV'.format(E_LUMO))
            print('Gap HOMO-LUMO = {:.2f} eV'.format(E_LUMO-E_HOMO))
if opt.job == 'gap' and ('XTB' in opt.method or 'XTB2' in opt.method) and opt.software == 'orca':
    E_HOMO = []
    E_LUMO = []
    with open('cmmd.out','r') as f:
        for line in f:
            if '(HOMO)' in line:
                arr = line.split()
                E_HOMO.append(float(arr[3]))
            if '(LUMO)' in line:
                arr = line.split()
                E_LUMO.append(float(arr[2]))
    E_gap = E_LUMO[-1] - E_HOMO[-1]
    print('######INFORMASI ENERGI HOMO & LUMO######')
    print('Energi HOMO = {:.2f} eV'.format(E_HOMO[-1]))
    print('Energi LUMO = {:.2f} eV'.format(E_LUMO[-1]))
    print('Gap HOMO-LUMO = {:.2f} eV'.format(E_gap))

if opt.job == 'td' and opt.software == 'orca':
    with open('cmmd.out','r') as f:
            lines = f.readlines()
            occ = []
            energy = []
            for index,line in enumerate(lines):
                if 'Basis Dimension' in line:
                    arr = line.split()
                    NBas = int(arr[4])
                if 'ORBITAL ENERGIES' in line:
                    for i in range(NBas):
                        arr = lines[index+i+4].split()
                        occ.append(arr[1])
                        energy.append(float(arr[3]))
            Eocc = []
            Evir = []
            for energy in energy:
                if energy > 0:
                    Evir.append(energy)
                if energy < 0:
                    Eocc.append(energy)
            E_HOMO = max(Eocc)
            E_LUMO = min(Evir)
            print('######INFORMASI ENERGI HOMO & LUMO######')
            print('Energi HOMO = {:.2f} eV'.format(E_HOMO))
            print('Energi LUMO = {:.2f} eV'.format(E_LUMO))
            print('Gap HOMO-LUMO = {:.2f} eV'.format(E_LUMO-E_HOMO))


    # A sqrt(2) * standard deviation of 0.4 eV is 3099.6 nm. 0.1 eV is 12398.4 nm. 0.2 eV is 6199.2 nm.
    stdev = 4132.806433333334
    # For Lorentzians, gamma is half bandwidth at half peak height (nm)
    gamma = 12.5
    # Excitation energies dalam nm
    #bands = [330,328,328,308,290,290,288,283,276,270,268]
    bands = []
    fosc = []
    with open("cmmd.out", 'r') as f:
        for line in f:
            if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                next(f)
                next(f)
                next(f)
                next(f)
                for i in range(opt.nroots):
                    arr = next(f).split()
                    bands.append(float(arr[2]))
                    fosc.append(float(arr[3]))

    #f = [7.90e-7,0.00,7.16e-4,1.02e-2,1.38e-6,2.94e-7,0.00,8.86e-4,1.54e-5,1.25e-2,9.31e-3]
    if len(bands) != len(fosc):
        print('Jumlah panjang gelombang tidak sama dengan jumlah oscillator strength')
        sys.exit()

    def KurvaGaussian(x, band, strength, stdev):
        "Memproduksi kurva Gaussian"
        bandshape = 1.3062974e8 * (strength / (1e7/stdev))  * np.exp(-(((1.0/x)-(1.0/band))/(1.0/stdev))**2)
    # Definisi di atas diambil dari P. J. Stephens, N. Harada, Chirality 22, 229 (2010)
        return bandshape

    def lorentzBand(x, band, strength, stdev, gamma):
        "Memproduksi kurva Lorentzian"
        bandshape= 1.3062974e8 * (strength / (1e7/stdev)) * ((gamma**2)/((x - band)**2 + gamma**2))
        return bandshape

    x = np.linspace(opt.start, opt.end, opt.points)

    composite = 0
    for count, peak in enumerate(bands):
        PuncakIni = KurvaGaussian(x, peak, fosc[count], stdev)
        composite += PuncakIni

    composite = np.array(composite)

    fig,ax = plt.subplots()
    ax.plot(x, composite,color='blue')
    plt.xlim(opt.start,opt.end)
    plt.xlabel('$\lambda$ [nm]')
    plt.ylabel('$\epsilon$ [L mol$^{-1}$cm$^{-1}$]')
    # plt.ticklabel_format(style="sci",scilimits=(0,0))
    plt.savefig('plot_uv.pdf')

if opt.job == 'td' and opt.software == 'dcdftb':
    nroots = 0
    with open('cmmd.out','r') as f:
            for line in f:
                if 'Excitation Energy' in line:
                    arr = line.split()
                    print("Gap HOMO-LUMO = {:.2f} eV".format(float(arr[3])*27.2114))
                    from scipy.constants import h, c
                    wavelength = h*c/(float(arr[3])*4.35974e-18)*1e9
                    print("Panjang gelombang serapan = {:.0f} nm".format(wavelength))
                if 'Number of states' in line:
                    arr = line.split()
                    nroots+=int(arr[4])
    print('Jumlah keadaan tereksitasi = {}'.format(nroots))

    # A sqrt(2) * standard deviation of 0.4 eV is 3099.6 nm. 0.1 eV is 12398.4 nm. 0.2 eV is 6199.2 nm.
    stdev = 4132.806433333334
    # For Lorentzians, gamma is half bandwidth at half peak height (nm)
    gamma = 12.5
    # Excitation energies dalam nm
    #bands = [330,328,328,308,290,290,288,283,276,270,268]
    bands = []
    fosc = []
    with open("cmmd.out", 'r') as f:
        for line in f:
            if 'Energies' in line:
                next(f)
                next(f)
                next(f)
                for i in range(nroots):
                    arr = next(f).split()
                    bands.append(float(arr[1]))
            if 'Strength' in line:
                next(f)
                for i in range(nroots):
                    arr = next(f).split()
                    fosc.append(float(arr[1]))
    #f = [7.90e-7,0.00,7.16e-4,1.02e-2,1.38e-6,2.94e-7,0.00,8.86e-4,1.54e-5,1.25e-2,9.31e-3]
    eVtoJoule = 1.60218e-19
    bands = np.array(bands)
    from scipy.constants import h,c
    bands = h*c/(bands*eVtoJoule)*1e9
    if len(bands) != len(fosc):
        print('Jumlah panjang gelombang tidak sama dengan jumlah oscillator strength')
        sys.exit()

    def KurvaGaussian(x, band, strength, stdev):
        "Memproduksi kurva Gaussian"
        bandshape = 1.3062974e8 * (strength / (1e7/stdev))  * np.exp(-(((1.0/x)-(1.0/band))/(1.0/stdev))**2)
    # Definisi di atas diambil dari P. J. Stephens, N. Harada, Chirality 22, 229 (2010)
        return bandshape

    def lorentzBand(x, band, strength, stdev, gamma):
        "Memproduksi kurva Lorentzian"
        bandshape= 1.3062974e8 * (strength / (1e7/stdev)) * ((gamma**2)/((x - band)**2 + gamma**2))
        return bandshape

    x = np.linspace(opt.start, opt.end, opt.points)

    composite = 0
    for count, peak in enumerate(bands):
        PuncakIni = KurvaGaussian(x, peak, fosc[count], stdev)
        composite += PuncakIni

    composite = np.array(composite)

    with open("uv.dat",'w') as f:
        for i, j in zip(x,composite):
            print(i,j,file=f)

    fig,ax = plt.subplots()
    ax.plot(x, composite,color='blue')
    plt.xlim(opt.start,opt.end)
    plt.ylim(0,)
    plt.xlabel('$\lambda$ [nm]')
    plt.ylabel('$\epsilon$ [L mol$^{-1}$cm$^{-1}$]')
    plt.ticklabel_format(style="sci",scilimits=(0,3))
    plt.savefig('plot_uv.pdf')

if opt.job == 'msdcom':
    msd_com(opt.groups,opt.traject,opt.start,opt.latt,opt.dt)

if opt.job == 'msdfit':
    msd_fit(opt.msd, opt.noheader,opt.start)

if opt.job == 'reax' or opt.inputreaction != 'None':
    Input = opt.inputreaction
    Reaksi = Input.split('->')
    Reaktan = Reaksi[0]
    Reaktan = Reaktan.split('+')

    KoefisienReaktan = []
    reak = []
    for i in Reaktan:
    	coeff = re.split(r"(\d+)", i)
    	KoefisienReaktan.append(coeff[1])
    	reak.append("".join(coeff[2:]))
    KoefisienReaktan = [-int(i) for i in KoefisienReaktan]
    Reaktan = reak
    Produk = Reaksi[1]
    Produk = Produk.split('+')
    KoefisienProduk = []
    prod = []
    for i in Produk:
    	coeff = re.split(r"(\d+)", i)
    	KoefisienProduk.append(coeff[1])
    	prod.append("".join(coeff[2:]))
    KoefisienProduk = [int(i) for i in KoefisienProduk]
    Produk = prod
    E_reaktan = []
    zpe_reaktan = []
    E_cor_reaktan = []
    H_cor_reaktan = []
    G_cor_reaktan = []
    E_produk = []
    zpe_produk = []
    E_cor_produk = []
    H_cor_produk = []
    G_cor_produk = []
    H_reaktan = []
    G_reaktan = []
    H_produk = []
    G_produk = []
    for i in Reaktan:
        with open('{}/cmmd.out'.format(i),'r') as f:
            if opt.software == 'orca':
                for line in f:
                    if "Total enthalpy" in line:
                        arr = line.split()
                        H_reaktan.append(float(arr[3]))
                    if "Final Gibbs free energy" in line:
                        arr = line.split()
                        G_reaktan.append(float(arr[5]))
                    if "Total thermal energy" in line:
                        arr = line.split()
                        E_reaktan.append(float(arr[3]))
            if opt.software == 'dcdftb':
            	temp_ener = []
            	for line in f:
            		if "Final" in line and "Energy" in line:
            			arr = line.split()
            			temp_ener.append(float(arr[4]))
            		if "Zero" in line:
            			arr = line.split()
            			zpe_reaktan.append(float(arr[4]))
            		if "Thermal correction to energy" in line:
            			arr = line.split()
            			E_cor_reaktan.append(float(arr[5]))
            		if "Thermal correction to enthalpy" in line:
            			arr = line.split()
            			H_cor_reaktan.append(float(arr[5]))
            		if "Thermal correction to Gibbs free energy" in line:
            			arr = line.split()
            			G_cor_reaktan.append(float(arr[7]))
            	E_reaktan.append(temp_ener[-1])


    for i in Produk:
        with open('{}/cmmd.out'.format(i),'r') as f:
            if opt.software == 'orca':
                for line in f:
                    if "Total enthalpy" in line:
                        arr = line.split()
                        H_produk.append(float(arr[3]))
                    if "Final Gibbs free energy" in line:
                        arr = line.split()
                        G_produk.append(float(arr[5]))
                    if "Total thermal energy" in line:
                        arr = line.split()
                        E_produk.append(float(arr[3]))
            if opt.software == 'dcdftb':
            	temp_ener = []
            	for line in f:
            		if "Final" in line and "Energy" in line:
            			arr = line.split()
            			temp_ener.append(float(arr[4]))
            		if "Zero" in line:
            			arr = line.split()
            			zpe_produk.append(float(arr[4]))
            		if "Thermal correction to energy" in line:
            			arr = line.split()
            			E_cor_produk.append(float(arr[5]))
            		if "Thermal correction to enthalpy" in line:
            			arr = line.split()
            			H_cor_produk.append(float(arr[5]))
            		if "Thermal correction to Gibbs free energy" in line:
            			arr = line.split()
            			G_cor_produk.append(float(arr[7]))
            	E_produk.append(temp_ener[-1])

    Koefisien = KoefisienReaktan + KoefisienProduk
    Koefisien = np.array(Koefisien,dtype=int)

    if opt.software == 'orca':
        H = H_reaktan + H_produk
        H = np.array(H,dtype=float)
        delta_H = (np.dot(H,Koefisien))*Ha2kj
        E = E_reaktan + E_produk
        E = np.array(E,dtype=float)
        delta_E = (np.dot(E,Koefisien))*Ha2kj
        print('Delta_E = {:.2f} kJ/mol'.format(delta_E))
        print('Delta_H = {:.2f} kJ/mol'.format(delta_H))
        G = G_reaktan + G_produk
        G = np.array(G,dtype=float)
        delta_G = (np.dot(G,Koefisien))*Ha2kj
        print('Delta_G = {:.2f} kJ/mol'.format(delta_G))
        delta_S = (delta_H - delta_G)/opt.temp*1000
        print('Delta_S = {:.2f} J/(mol K)'.format(delta_S))
    if opt.software == 'dcdftb':
    	E = np.array(E_reaktan + E_produk, dtype=float)

    	E_cor = np.array(E_cor_reaktan + E_cor_produk, dtype=float)
    	H_cor = np.array(H_cor_reaktan + H_cor_produk, dtype=float)
    	G_cor = np.array(G_cor_reaktan + G_cor_produk, dtype=float)
    	zpe = np.array(zpe_reaktan + zpe_produk, dtype=float)
    	delta_zpe = (np.dot(zpe,Koefisien))*Ha2kj
    	delta_El = (np.dot(E,Koefisien))*Ha2kj
    	delta_Ecor = (np.dot(E_cor,Koefisien))*Ha2kj

    	delta_Hcor = (np.dot(H_cor,Koefisien))*Ha2kj
    	delta_Gcor = (np.dot(G_cor,Koefisien))*Ha2kj
    	delta_E = delta_El + delta_zpe + delta_Ecor
    	delta_H = delta_El + delta_zpe + delta_Hcor
    	delta_G = delta_El + delta_zpe + delta_Gcor

    	print('Delta_E = {:.2f} kJ/mol'.format(delta_E))
    	print('Delta_H = {:.2f} kJ/mol'.format(delta_H))
    	print('Delta_G = {:.2f} kJ/mol'.format(delta_G))
    	delta_S = (delta_H - delta_G)/opt.temp*1000
    	print('Delta_S = {:.2f} J/(mol K)'.format(delta_S))

if opt.job == 'ts' and opt.software == 'orca':
    Input = opt.inputreaction
    Reaksi = Input.split('->')
    Reaktan = Reaksi[0]
    Reaktan = Reaktan.split('+')
    KoefisienReaktan = [i[0] for i in Reaktan]
    KoefisienReaktan = [-int(i) for i in KoefisienReaktan]
    Reaktan = [i[1:] for i in Reaktan]
    TS = Reaksi[1]
    TS = TS.split('+')
    KoefisienTS = [i[0] for i in TS]
    KoefisienTS = [int(i) for i in KoefisienTS]
    TS = [i[1:] for i in TS]
    Produk = Reaksi[2]
    Produk = Produk.split('+')
    KoefisienProduk = [i[0] for i in Produk]
    KoefisienProduk = [-int(i) for i in KoefisienProduk]
    Produk = [i[1:] for i in Produk]

    H_reaktan = []
    H_TS = []
    H_produk = []

    G_reaktan = []
    G_TS = []
    G_produk = []

    E_reaktan = []
    E_TS = []
    E_produk = []

    for i in Reaktan:
        with open('{}/cmmd.out'.format(i),'r') as f:
            for line in f:
                if "Total enthalpy" in line:
                    arr = line.split()
                    H_reaktan.append(float(arr[3]))
                if "Final Gibbs free energy" in line:
                    arr = line.split()
                    G_reaktan.append(float(arr[5]))
                if "Total thermal energy" in line:
                    arr = line.split()
                    E_reaktan.append(float(arr[3]))

    for i in TS:
        with open('{}/cmmd.out'.format(i),'r') as f:
            for line in f:
                if "Total enthalpy" in line:
                    arr = line.split()
                    H_TS.append(float(arr[3]))
                if "Final Gibbs free energy" in line:
                    arr = line.split()
                    G_TS.append(float(arr[5]))
                if "Total thermal energy" in line:
                    arr = line.split()
                    E_TS.append(float(arr[3]))
    for i in Produk:
        with open('{}/cmmd.out'.format(i),'r') as f:
            for line in f:
                if "Total enthalpy" in line:
                    arr = line.split()
                    H_produk.append(float(arr[3]))
                if "Final Gibbs free energy" in line:
                    arr = line.split()
                    G_produk.append(float(arr[5]))
                if "Total thermal energy" in line:
                    arr = line.split()
                    E_produk.append(float(arr[3]))

    Koefisien_forward = KoefisienReaktan + KoefisienTS
    Koefisien_forward = np.array(Koefisien_forward,dtype=int)
    Koefisien_backward = KoefisienProduk + KoefisienTS
    Koefisien_backward = np.array(Koefisien_backward,dtype=int)
    E_forward = E_reaktan + E_TS
    E_forward = np.array(E_forward,dtype=float)
    E_backward = E_produk + E_TS
    E_backward = np.array(E_backward,dtype=float)
    delta_E_forward = (np.dot(E_forward,Koefisien_forward))*Ha2kj
    delta_E_backward = (np.dot(E_backward,Koefisien_backward))*Ha2kj
    print('Ea_forward* = {:.2f} kJ/mol'.format(delta_E_forward))
    print('Ea_backward* = {:.2f} kJ/mol'.format(delta_E_backward))

    H_forward = H_reaktan + H_TS
    H_forward = np.array(H_forward,dtype=float)
    H_backward = H_produk + H_TS
    H_backward = np.array(H_backward,dtype=float)
    delta_H_forward = (np.dot(H_forward,Koefisien_forward))*Ha2kj
    delta_H_backward = (np.dot(H_backward,Koefisien_backward))*Ha2kj
    print('Delta_H_forward* = {:.2f} kJ/mol'.format(delta_H_forward))
    print('Delta_H_backward* = {:.2f} kJ/mol'.format(delta_H_backward))

    G_forward = G_reaktan + G_TS
    G_forward = np.array(G_forward,dtype=float)
    G_backward = G_produk + G_TS
    G_backward = np.array(G_backward,dtype=float)
    delta_G_forward = (np.dot(G_forward,Koefisien_forward))*Ha2kj
    delta_G_backward = (np.dot(G_backward,Koefisien_backward))*Ha2kj
    print('Delta_G_forward* = {:.2f} kJ/mol'.format(delta_G_forward))
    print('Delta_G_backward* = {:.2f} kJ/mol'.format(delta_G_backward))
    delta_S_forward = (delta_H_forward - delta_G_forward)/opt.temp*1000
    delta_S_backward = (delta_H_backward - delta_G_backward)/opt.temp*1000
    print('Delta_S_forward* = {:.2f} J/(mol K)'.format(delta_S_forward))
    print('Delta_S_backward* = {:.2f} J/(mol K)'.format(delta_S_backward))

if opt.job == 'force' and opt.software == 'orca':
    with open('cmmd.xyz') as f:
        x = int(f.readlines()[0])
    x = round(x*3/4)
    with open('cmmd.opt') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "$gradients" in line:
                data = lines[i+2:i+x+2]
    data = [x.strip("\n").strip(" ") for x in data]
    data_clean = []
    for i in data:
        x = i.split(" ")
        for i in x:
            data_clean.append(i.strip(""))
    # data_clean = [float(x) for x in data_clean]
    data = []
    for i in data_clean:
        if i != '':
            data.append(float(i))
    
    j = 0
    for j in range(int(len(data)/3)):
        if j < (len(data) - 2): 
            print(data[j],data[j+1],data[j+2])
            j += 3
        
        
    
    
    
    # data = [x.strip(" ") for x in data]
    # data = [x.strip(",") for x in data]
    
            



if opt.job == 'opt' and opt.software == 'dcdftb':
    Energy = []
    Gradient = []
    if os.path.isfile('cmmd.out'):
        with open(opt.input,'r') as f:
            Natom = int(next(f))
        with open('cmmd.out', 'r') as f:
            koord = []
            lines = f.readlines()
            for i,line in enumerate(lines):
                if 'Final molecular coordinate [Angstrom]' in line:
                    for j in range(Natom):
                        koord.append(lines[i+5+j].strip('\n'))
        with open('cmmd.xyz', 'w') as fout:
            print(Natom,file=fout)
            print('File generated by CMMDE output parser',file=fout)
            for atom in koord:
                print(atom, file=fout)
    else:
        print("Mohon menunggu hingga perhitungan anda selesai.")
        exit

if opt.job == 'rdf':
    rdf(opt.traject,opt.latt,opt.pair,opt.range,opt.resolution)

if opt.job == 'cv' and opt.software == 'dcdftb':
    time = []
    cv = []
    gaus = []
    with open('biaspot','r') as f:
        for line in f:
            if "THIS RUN'S STEP NO.=" in line:
                arr = line.split()
                time.append(int(arr[9])/1000.)
            if "GAUSSIAN BIAS POTENTIAL:" in line:
                arr = line.split()
                gaus.append(arr[3])
            if "Coordinate" in line:
                arr = line.split()
                cv.append(float(arr[2]))

    with open('cv.dat','w') as f:
        for i,j in zip(time,cv):
            print(i,j,file=f)
    time = np.array(time)
    # Plot Collective Variable
    fig,ax = plt.subplots()
    ax2 = ax.twiny()
    ax.plot(time,cv,color='blue')
    ax.tick_params(direction='in')
    ax2.tick_params(direction='in')

    def tick_func(metafreq,step):
        gauss = step/metafreq
        return [int(x)-1 for x in gauss]
    ax.set_xlim(time[0],time[-1])
    ax2.set_xlim(ax.get_xlim())
    newtics = np.arange(min(time),max(time),time[-1]/10)
    ax2.set_xticks(newtics)
    ax2.set_xticklabels(tick_func(opt.metafreq,newtics*1000))
    ax2.set_xlabel('Number of deposited Gaussian')
    plt.grid(linestyle=':')
    plt.yticks(np.arange(min(cv),max(cv),0.1))
    ax.set_xlabel('Time [ps]')
    if opt.cvtype == 'coordnum':
        ax.set_ylabel('Coordination Number')
    if opt.cvtype == 'distance':
        ax.set_ylabel('Distance [Angstroms]')
    if opt.cvtype == 'angle':
        ax.set_ylabel('Angle [degree]')
    if opt.cvtype == 'dihedral':
        ax.set_ylabel('Dihedral [degree]')
    if opt.cvtype == 'distancediff':
        ax.set_ylabel('Distance diff. [Angstroms]')
    if opt.cvtype == 'distanceadd':
        ax.set_ylabel('Distance add. [Angstroms]')
    if opt.cvtype == 'meandistance':
        ax.set_ylabel('Mean distance [Angstroms]')
    if opt.cvtype == 'pointplanedistance':
        ax.set_ylabel('Point plane distance [Angstroms]')
    plt.savefig('cv.pdf',dpi=1200,format='pdf')

if opt.job == 'fes' and opt.software == 'dcdftb':
    fesdata = []
    with open('fes.dat','r') as f:
        lines = f.readlines()
        for index,line in enumerate(lines):
            if 'FREE ENERGY SURFACE CONSISTING OF     {} GAUSSIANS'.format(opt.ngaus) in line:
                for i in range(int((opt.fesend-opt.fesstart)/opt.fesbin)+1):
                    fesdata.append(lines[index+1+i].strip())
    with open('fes_{}.dat'.format(opt.ngaus),'w') as f:
        for i in fesdata:
            print(i,file=f)

    with open('fes_{}.dat'.format(opt.ngaus),'r') as f:
        cv = []
        fes = []
        for line in f:
            arr = line.split()
            cv.append(float(arr[0]))
            fes.append(float(arr[1]))
    fes = np.array(fes)
    # Plot free energy surface
    fig,ax = plt.subplots()
    ax.plot(cv,fes*2625.5-min(fes*2625.5),color='blue')
    ax.tick_params(direction='in')

    ax.set_xlim(cv[0],cv[-1])
    plt.grid(linestyle=':')
    plt.ylabel('Free energy [kJ/mol]')
    if opt.cvtype == 'coordnum':
        plt.xlabel('Coordination Number')
    if opt.cvtype == 'distance':
        plt.xlabel('Distance [Angstroms]')
    if opt.cvtype == 'angle':
        plt.xlabel('Angle [degree]')
    if opt.cvtype == 'dihedral':
        plt.xlabel('Dihedral [degree]')
    if opt.cvtype == 'distancediff':
        plt.xlabel('Distance diff. [Angstroms]')
    if opt.cvtype == 'distanceadd':
        plt.xlabel('Distance add. [Angstroms]')
    if opt.cvtype == 'meandistance':
        plt.xlabel('Mean distance [Angstroms]')
    if opt.cvtype == 'pointplanedistance':
        plt.xlabel('Point plane distance [Angstroms]')
    plt.savefig('fes.pdf',dpi=1200,format='pdf')

if opt.job == 'barrier' and opt.software == 'dcdftb':
    from scipy.signal import argrelextrema
    kJ2eV = 1.0364e-2
    x = []
    y = []
    with open('fes_{}.dat'.format(opt.ngaus),'r') as f:
        for line in f:
            arr = line.split()
            x.append(float(arr[0]))
            y.append(float(arr[1]))
    x = np.array(x)
    y = np.array(y)
    sortid = np.argsort(x)
    x = x[sortid]
    y = y[sortid]

    maxm = argrelextrema(y, np.greater)
    minm = argrelextrema(y, np.less)
    maxval = [float(i) for i in y[maxm]]
    minval = [float(i) for i in y[minm]]

    # print('Minimum absis: {}'.format(x[minm]))
    # print('Maximum absis: {}'.format(x[maxm]))
    # print('Minimum ordinate: {}'.format(y[minm]))
    # print('Maximum ordinate: {}'.format(y[maxm]))

    ###Activation barrier for forward and backward reactions

    if (len(minval) == 3):
        # If exists two transition states
        Ef1 = (maxval[1]-minval[2])*2625.5 # in kJ/mol
        Ef2 = (maxval[0]-minval[1])*2625.5
        Eb1 = (maxval[0]-minval[0])*2625.5
        Eb2 = (maxval[1]-minval[1])*2625.5
        ##Now calculating delta G
        deltaG = (minval[-1]-minval[0])*2625.5 # in kJ/mol
        print('delta_G: {:.2f} kJ/mol ({:.2f} eV)'.format(deltaG, deltaG*kJ2eV))
        print('delta_G*(backward1): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef1, Ef1*kJ2eV))
        print('delta_G*(backward2): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef2, Ef2*kJ2eV))
        print('delta_G*(forward1): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb1, Eb1*kJ2eV))
        print('delta_G*(forward2): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb2, Eb2*kJ2eV))
    elif (len(minval) == 1):
        Ef1 = (maxval[0]-minval[0])*2625.5
        Eb1 = (maxval[0]-y[0])*2625.5
        print('delta_G*(forward): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef1, Ef1*kJ2eV))
        print('delta_G*(backward): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb1, Eb1*kJ2eV))
        deltaG = (y[0]-minval[0])*2625.5
        print('delta_G: {:.2f} kJ/mol ({:.2f} eV)'.format(deltaG, deltaG*kJ2eV))
    else:
        Ef1 = (maxval[0]-minval[1])*2625.5 # in kJ/mol
        Ef2 = 0
        Eb1 = (maxval[0]-minval[0])*2625.5
        Eb2 = 0
        ##Now calculating delta G
        deltaG = (minval[-1]-minval[0])*2625.5 # in kJ/mol
        print('delta_G: {:.2f} kJ/mol ({:.2f} eV)'.format(deltaG, deltaG*kJ2eV))
        print('delta_G*(backward1): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef1, Ef1*kJ2eV))
        print('delta_G*(backward2): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef2, Ef2*kJ2eV))
        print('delta_G*(forward1): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb1, Eb1*kJ2eV))
        print('delta_G*(forward2): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb2, Eb2*kJ2eV))

if opt.job == 'rigiddock':
    path = './'
    files = []
    ligands = []
    scores = []
    vdw = []
    es = []
    rep = []
    rmsd = []
    mw = []
    for i in os.listdir(path):
        if os.path.isdir(os.path.join(path,i)) and 'RigidDock' in i:
            files.append(i)
    for i in files:
        ligands.append(i.split("RigidDock_")[1])
        with open('{}/rigid_scored.mol2'.format(i),'r') as f:
            for line in f:
                if 'Grid_Score:' in line:
                    arr = line.split()
                    scores.append(round(float(arr[2]),2))
                if 'Grid_vdw_energy:' in line:
                    arr = line.split()
                    vdw.append(round(float(arr[2]),2))
                if 'Grid_es_energy:' in line:
                    arr = line.split()
                    es.append(round(float(arr[2]),2))
                if 'Internal_energy_repulsive:' in line:
                    arr = line.split()
                    rep.append(round(float(arr[2]),2))
                if 'HA_RMSDh' in line:
                    arr = line.split()
                    rmsd.append(round(float(arr[2]),2))
                if 'Molecular_Weight' in line:
                    arr = line.split()
                    mw.append(round(float(arr[2]),2))
    if len(scores) < 20:
        t = PrettyTable()
        t.add_column("Ligands",ligands)
        # t.add_column("Mw [g/mol]",mw)
        t.add_column("Scores [kJ/mol]",scores)
        t.add_column("VDW [kJ/mol]",vdw)
        t.add_column("Elec. [kJ/mol]",es)
        # t.add_column("Rep.  [kJ/mol]",rep)
        if len(rmsd) > 0:
            t.add_column("RMSD [Angs.]",rmsd)
        print(t.get_string(sortby="Scores [kJ/mol]"))
    else:
        with open('DockTable.dat','w') as f:
            t = PrettyTable()
            t.add_column("Ligands",ligands)
            # t.add_column("Mw [g/mol]",mw)
            t.add_column("Scores [kJ/mol]",scores)
            t.add_column("VDW [kJ/mol]",vdw)
            t.add_column("Elec. [kJ/mol]",es)
            # t.add_column("Rep.  [kJ/mol]",rep)
            if len(rmsd) > 0:
                t.add_column("RMSD [Angs.]",rmsd)
            print(t.get_string(sortby="Scores [kJ/mol]"),file=f)
if opt.job == 'flexdock':
    path = './'
    files = []
    ligands = []
    scores = []
    vdw = []
    es = []
    rep = []
    mw = []
    rmsd = []
    for i in os.listdir(path):
        if os.path.isdir(os.path.join(path,i)) and 'FlexDock' in i:
            if os.path.exists('{}/flex_scored.mol2'.format(i)):
                if os.stat('{}/flex_scored.mol2'.format(i)).st_size !=0:
                    files.append(i)
    for i in files:
        ligands.append(i.split("FlexDock_")[1])
        with open('{}/flex_scored.mol2'.format(i),'r') as f:
            for line in f:
                if 'Grid_Score:' in line:
                    arr = line.split()
                    scores.append(round(float(arr[2]),2))
                if 'Grid_vdw_energy:' in line:
                    arr = line.split()
                    vdw.append(round(float(arr[2]),2))
                if 'Grid_es_energy:' in line:
                    arr = line.split()
                    es.append(round(float(arr[2]),2))
                if 'Internal_energy_repulsive:' in line:
                    arr = line.split()
                    rep.append(round(float(arr[2]),2))
                if 'HA_RMSDh' in line:
                    arr = line.split()
                    rmsd.append(round(float(arr[2]),2))
                if 'Molecular_Weight' in line:
                    arr = line.split()
                    mw.append(round(float(arr[2]),2))
    if len(scores) < 20:
        t = PrettyTable()
        t.add_column("Ligands",ligands)
        # t.add_column("Mw [g/mol]",mw)
        t.add_column("Scores [kJ/mol]",scores)
        t.add_column("VDW [kJ/mol]",vdw)
        t.add_column("Elec. [kJ/mol]",es)
        # t.add_column("Rep.  [kJ/mol]",rep)
        if len(rmsd) > 0:
            t.add_column("RMSD [Angs.]",rmsd)
        print(t.get_string(sortby="Scores [kJ/mol]"))
    else:
        with open('DockTable.dat','w') as f:
            t = PrettyTable()
            t.add_column("Ligands",ligands)
            # t.add_column("Mw [g/mol]",mw)
            t.add_column("Scores [kJ/mol]",scores)
            t.add_column("VDW [kJ/mol]",vdw)
            t.add_column("Elec. [kJ/mol]",es)
            # t.add_column("Rep.  [kJ/mol]",rep)
            if len(rmsd) > 0:
                t.add_column("RMSD [Angs.]",rmsd)
            print(t.get_string(sortby="Scores [kJ/mol]"),file=f)

if opt.job == 'checkopt' and opt.software == 'dock':
    checkopt(opt.nligands)

# Rescoring docking atau perhitungan lainnya yang melibatkan perhitungan satu titik menggunakan DCDFTB
if opt.job == 'dock' and opt.software == 'dcdftb':
    Input = opt.inputreaction
    Reaksi = Input.split('->')
    Reaktan = Reaksi[0]
    Reaktan = Reaktan.split('+')
    KoefisienReaktan = [i[0] for i in Reaktan]
    KoefisienReaktan = [-int(i) for i in KoefisienReaktan]
    Reaktan = [i[1:] for i in Reaktan]
    Produk = Reaksi[1]
    Produk = Produk.split('+')
    KoefisienProduk = [i[0] for i in Produk]
    KoefisienProduk = [int(i) for i in KoefisienProduk]
    Produk = [i[1:] for i in Produk]

    E_reaktan = []
    E_produk = []

    for i in Reaktan:
        with open('{}/cmmd.out'.format(i),'r') as f:
            for line in f:
                if "Final DC-DFTB-3rd Energy"in line:
                    arr = line.split()
                    E_reaktan.append(float(arr[4]))
                elif "Final DFTB-3rd Energy" in line:
                    arr = line.split()
                    E_reaktan.append(float(arr[4]))
                elif "Final SCC-DFTB Energy" in line:
                    arr = line.split()
                    E_reaktan.append(float(arr[4]))

    for i in Produk:
        with open('{}/cmmd.out'.format(i),'r') as f:
            for line in f:
                if "Final DC-DFTB-3rd Energy" in line:
                    arr = line.split()
                    E_produk.append(float(arr[4]))
                elif "Final DFTB-3rd Energy" in line:
                    arr = line.split()
                    E_produk.append(float(arr[4]))
                elif "Final SCC-DFTB Energy" in line:
                    arr = line.split()
                    E_produk.append(float(arr[4]))

    Koefisien = KoefisienReaktan + KoefisienProduk
    Koefisien = np.array(Koefisien,dtype=int)
    E = E_reaktan + E_produk
    E = np.array(E,dtype=float)
    delta_E = (np.dot(E,Koefisien))*Ha2kj
    print('Delta_E = {:.2f} kJ/mol'.format(delta_E))

if opt.job == 'opt' and opt.software == 'dftb':
    os.system('cmmde_gen2poscar.py cmmd.gen > cmmd.vasp')
if opt.job == 'reax' and opt.software == 'dftb':
    Input = opt.inputreaction
    Reaksi = Input.split('->')
    Reaktan = Reaksi[0]
    Reaktan = Reaktan.split('+')
    KoefisienReaktan = [i[0] for i in Reaktan]
    KoefisienReaktan = [-int(i) for i in KoefisienReaktan]
    Reaktan = [i[1:] for i in Reaktan]
    Produk = Reaksi[1]
    Produk = Produk.split('+')
    KoefisienProduk = [i[0] for i in Produk]
    KoefisienProduk = [int(i) for i in KoefisienProduk]
    Produk = [i[1:] for i in Produk]

    E_reaktan = []
    E_produk = []

    for i in Reaktan:
        with open('{}/detailed.out'.format(i),'r') as f:
            for line in f:
                if "Total energy:"in line:
                    arr = line.split()
                    E_reaktan.append(float(arr[2]))

    for i in Produk:
        with open('{}/detailed.out'.format(i),'r') as f:
            for line in f:
                if "Total energy:" in line:
                    arr = line.split()
                    E_produk.append(float(arr[2]))

    Koefisien = KoefisienReaktan + KoefisienProduk
    Koefisien = np.array(Koefisien,dtype=int)
    E = E_reaktan + E_produk
    E = np.array(E,dtype=float)
    delta_E = (np.dot(E,Koefisien))*Ha2kj
    print('Delta_E = {:.2f} kJ/mol'.format(delta_E))

if opt.job == 'cell':
    os.system('aflow --data < {}'.format(opt.input)) # Untuk menggunakan fitur ini, install aflow, wget http://aflow.org/install-aflow/install-aflow.sh


if opt.job == 'dos' and opt.software == 'dftb':
    os.system('$DPTOOLS_DIR/dp_dos band.out dos_total.dat')

if opt.job == 'pdos' and opt.software == 'dftb':
    filename = opt.pdos.split('.out')[0]
    os.system('$DPTOOLS_DIR/dp_dos -w {} {}.dat'.format(opt.pdos,filename))
# Plot file x y sembarang dengan ekstensi .dat
if opt.job == 'plot':
    x = []
    y = []
    with open(opt.file, 'r') as f:
        for line in f:
            arr = line.split()
            x.append(float(arr[0]))
            y.append(float(arr[1]))
    x = np.array(x)
    y = np.array(y)
    fig,ax = plt.subplots()
    ax.plot(x,y,color='blue')
    ax.tick_params(direction='in')
    # fig,ax = plt.subplots()
    plt.grid(linestyle=':')
    plt.ylim(0,max(y))
    # plt.xticks(np.arange(0,len(x)+1,25))
    # plt.yticks(np.arange(0,max(y),2))
    # ax.tick_params(direction='in')

    # ax.set_xlim(min(x),max(x))
    # ax.set_ylim(min(y),max(y))
    plt.ylabel(opt.ylabel)
    plt.xlabel(opt.xlabel)
    filename = opt.file.split('.dat')[0]
    plt.savefig('{}.pdf'.format(filename),dpi=300,format='pdf')

if opt.job == 'nci2d':
    with open('cmmd.nci','w') as f:
        print("""1
cmmd.wfn
{}
RADIUS 0. 0. 0. 2.
RANGE 3
-0.5 -0.02
-0.02 0.02
 0.02 0.5""".format(opt.grid),file=f)
    with open('run_nci.sh','w') as f:
        print("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:0:0
export OMP_NUM_THREADS={}
cd $PWD
$NCI_COMMAND cmmd.nci""".format(opt.nproc),file = f)
    os.system('sbatch run_nci.sh')

if opt.job == 'nciplot':
    with open('cmmd.gp','w') as f:
        print("""set terminal pngcairo size 1000,1000 enhanced font 'Helvetica,20'
set encoding iso_8859_1
set output 'nci.png'
set key
set ylabel 's(a.u.)' font "Helvetica, 30"
set xlabel 'sign({/Symbol l}_2){/Symbol r}(a.u.)' font "Helvetica, 30"
set pm3d map
# Define a color gradient palette used by pm3d
set palette defined (-0.04 "blue",0.00 "green", 0.04 "red")
set format y "% .2f"
set format x "% .2f"
set format cb "% -.2f"
set border lw 4
set xtic  -0.06,0.01,0.06 nomirror rotate font "Helvetica"
set ytic   0.0,0.25,1.0 nomirror font "Helvetica"
# set the color bar tics
set cbtic  -0.06,0.01,0.06 nomirror font "Helvetica"
set xrange [-0.06:0.06]
set yrange [0.0:1.0]
# set the range of values which are colored using the current palette
set cbrange [-0.06:0.06]
plot 'cmmd.dat' u 1:2:1 w p lw 6 palette t '' """,file = f)
    os.system('gnuplot cmmd.gp')

if 'dos' in opt.job and opt.software == 'qe':
    with open('cmmd.in','w') as f:
            print("""DOS
 prefix = 'cmmd_dos',
 outdir = {},
 fildos = 'cmmd.dos',
 emin = {},
 emax = {},
 /""".format(opt.outdir,opt.emin,opt.emax),file=f)
    os.system("$QE_DOS_COMMAND < cmmd.in > cmmd.out") # Hasil perhitungan berupa cmmd.dos

if 'charge' in opt.job and opt.software == 'qe':
    with open("cmmd.in",'w') as f:
            print("""&INPUTPP
 outdir = {},
 prefix = 'cmmd_charge',
 plot_num = 0,
/
&PLOT
 iflag = 3,
 output_format = 6,
 fileout = 'cmmd_rho.cube',
 nx = 64, ny = 64, nz = 64,
 /""".format(opt.outdir),file=f)
    with open('run.sh','w') as fout:
        print("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:0:0
export OMP_NUM_THREADS=1
cd $PWD
$QE_PP_COMMAND < cmmd.in > cmmd.out""",file=fout)
    os.system('sbatch run.sh')

if opt.job == 'opt' and opt.software == 'qe':
    with open('cmmd.out', 'r') as f:
        lines = f.readlines()
        for i,line in enumerate(lines):
            if "ATOMIC_POSITIONS" in line:
                print(lines[i+1])
