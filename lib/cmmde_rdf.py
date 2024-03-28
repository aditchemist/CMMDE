import warnings
warnings.filterwarnings("ignore")
import MDAnalysis
import MDAnalysis.analysis.rdf
import numpy as np
from pymatgen.core import lattice
import builtins
import os

def rdf(traject,latt,pair,endrange,resolution):
    cell = []
    with open(latt, 'r') as f:
        for line in f:
            if 'TV' in line:
                cell.append(list(map(float, line.split()[1:4])))

    box = lattice.Lattice(cell)
    # (lx, ly, lz), (a, b, c) = box.lengths_and_angles
    lx, ly, lz = box.lengths
    a, b, c = box.angles
    u = MDAnalysis.Universe(traject, format='XYZ')
    u.dimensions = np.array([lx, ly, lz, a, b, c], dtype=np.float32)

    builtins.u = u

    ag1 = None
    ag2 = None
    elems = pair.split('-')

    # Definisi file keluaran
    outfile = open('rdf_{}-{}.dat'.format(elems[0],elems[1]), 'w')

    ag1 = u.select_atoms(f'name {elems[0]}')
    ag2 = u.select_atoms(f'name {elems[1]}')

    nbins = int((endrange - 0.05)/resolution)

    rdf = MDAnalysis.analysis.rdf.InterRDF(ag1, ag2, nbins=nbins, range=(0.05, endrange),   verbose=True)
    rdf.run()

    vol = np.power(rdf.edges[1:], 3) - np.power(rdf.edges[:-1], 3)
    vol *= 4/3.0 * np.pi   
    box_vol = box.volume
    pair_den = box_vol / (len(ag1.atoms) * len(ag2.atoms))
    integral = rdf.rdf / (pair_den / vol)

    total = 0.0
    step = 0
    for r, v, d in zip(rdf.bins, rdf.rdf, integral):
        step += 1
        total += d/len(ag1.atoms)
        if (step == 1):
            continue
        print('{:8.4f} {:8.4f} {:8.4f}'.format(r, v, total), file=outfile)

    outfile.close()

    with open('rdf_{}-{}.gp'.format(elems[0],elems[1]),'w') as fout:
        print("""set terminal pdf
        set output "rdf_{}-{}.pdf"
        set key
        set yrange [0:10]
        set ytics (0,1,2,3,4,5,6,7,8,9,10)
        set xlabel "Distance [Angstroms]"
        set ylabel "RDF [arb.]"
        set style line 1 lc rgb '#00008b' lt 1 lw 2 pt 7
        set style line 2 lc rgb '#FF0000' lt 1 lw 2 pt 7
        plot 'rdf_{}-{}.dat' u 1:2 with line title "RDF" linestyle 1,\
            '' u 1:3 with line title "Coordination Number" linestyle 2""".format(elems[0],elems[1],elems[0],elems[1]),file=fout)
    os.system('gnuplot rdf_{}-{}.gp'.format(elems[0],elems[1]))