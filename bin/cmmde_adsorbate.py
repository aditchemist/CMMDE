#!/usr/bin/env python3
import pymatgen.analysis.adsorption as pa
import pymatgen.core.structure as st
from pymatgen.core import Structure
import argparse
import sys
import warnings
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(description='Build possible adsorbate-adsorbent interactions for the surface reaction')
parser.add_argument('-s', '--slab', type=str, help='structure POSCAR of the surface' )
parser.add_argument('-height','--height', type=float, help='height criteria for selection of surface sites' )
parser.add_argument('-d','--distance', type=float, help='distance from the coordinating ensemble of atoms along the miller index for the site (i. e. the distance from the slab itself)' )
parser.add_argument('-pi','--put_inside', default=True, type=bool, help='whether to put the site inside the cell ' )
parser.add_argument('-sr','--symm_reduce', type=float, default=0.1, help='symmetry reduction threshold' )
parser.add_argument('-nr','--near_reduce', type=float, default=0.1, help='near reduction threshold' )
parser.add_argument('-noh','--no_obtuse_hollow', type=bool, default=True, help='flag to indicate whether to include obtuse triangular ensembles in hollow sites' )
parser.add_argument('-sel','--selective_dynamics', type=bool, default=True, help='create the selective dynamics options for VASP calculations')
parser.add_argument('-dyn','--dyn', type=float, default=3.50, help='Height for dynamical atoms selections')
parser.add_argument('-t','--type', type=str, help='Type of adsorption site based on the list' )
parser.add_argument('-ad','--ad', type=str, help='Structure of adsorbate from file' )
parser.add_argument('-all','--all', type=bool, default=False, help='Print all possible structures')

opt = parser.parse_args(sys.argv[1:])

print('Begin searching for adsorption sites...')
# Defining the slab coordinate
slab_coord = Structure.from_file(opt.slab)

# Searching the active sites
check_site = pa.AdsorbateSiteFinder(slab_coord, height=opt.height)
# Preparing the adsorbate molecules in xyz cartesian format
ad = st.Molecule.from_file(opt.ad)

# Defining thke protocol for selective dynamics. Here you need to manually specify the height measured from the top of your slab structure. The default value is 3.50 Angstroms.
adsorbent = pa.AdsorbateSiteFinder(slab_coord, height=opt.dyn, selective_dynamics=opt.selective_dynamics)
adsorbent.slab.to(filename='adsorbent.vasp', fmt='poscar')
# Here, I assume that you want to generate one by one the adsorption types. The choices are 'ontop', 'hollow', and 'bridge' types.
if (opt.all == False):
    sites = check_site.find_adsorption_sites(distance=opt.distance, put_inside=opt.put_inside,  symm_reduce=opt.symm_reduce, near_reduce=opt.near_reduce, positions=[opt.type],      no_obtuse_hollow=opt.no_obtuse_hollow)
    site_coord = sites.get(opt.type)
    index = 0
    site_all = []
    for site in site_coord:
        adsorbed = adsorbent.add_adsorbate(ad, site, repeat=None, reorient=True)
        adsorbed.to(filename='{}_{}.vasp'.format(opt.type,index), fmt='poscar')
        index += 1
        site_all.append(site.tolist())
    # adsorbed = adsorbent.add_adsorbate(ad, site_all, repeat=None, reorient=True)
    # adsorbed.to(filename='{}_all.vasp'.format(opt.type), fmt='poscar')
    print('{}: {}'.format(opt.type, len(site_coord)))
    print('Done for all possibilities of {} adsorption mode!'.format(opt.type))
# The following if in case you want to generate all possible adsorption sites.
else:
    sites = check_site.find_adsorption_sites(distance=opt.distance, put_inside=opt.put_inside,  symm_reduce=opt.symm_reduce, near_reduce=opt.near_reduce, positions=['hollow','ontop','bridge'],      no_obtuse_hollow=opt.no_obtuse_hollow)
    hollow_coord = sites.get('hollow')
    bridge_coord = sites.get('bridge')
    ontop_coord = sites.get('ontop')
    hollow_ind = 0
    bridge_ind = 0
    ontop_ind = 0

    for site in hollow_coord:
        adsorbed = adsorbent.add_adsorbate(ad, site, repeat=None, reorient=True)
        elements = adsorbed.types_of_specie
        element = [str(x) for x in elements]
        sorted_ads = adsorbed.get_sorted_structure(key=lambda x: element.index(str(x.specie)))
        sorted_ads.to(filename='hollow_{}.vasp'.format(hollow_ind), fmt='poscar')
        hollow_ind += 1
    for site in ontop_coord:
        adsorbed = adsorbent.add_adsorbate(ad, site, repeat=None, reorient=True)
        elements = adsorbed.types_of_specie
        element = [str(x) for x in elements]
        sorted_ads = adsorbed.get_sorted_structure(key=lambda x: element.index(str(x.specie)))
        sorted_ads.to(filename='ontop_{}.vasp'.format(ontop_ind), fmt='poscar')
        ontop_ind += 1
    for site in bridge_coord:
        adsorbed = adsorbent.add_adsorbate(ad, site, repeat=None, reorient=True)
        elements = adsorbed.types_of_specie
        element = [str(x) for x in elements]
        sorted_ads = adsorbed.get_sorted_structure(key=lambda x: element.index(str(x.specie)))
        sorted_ads.to(filename='bridge_{}.vasp'.format(bridge_ind), fmt='poscar')
        bridge_ind += 1
    print('###############THE LIST OF ADSORPTION SITES####################')
    print('hollow:{}, bridge:{}, ontop:{}'.format(len(hollow_coord), len(bridge_coord), len(ontop_coord)))
    print('Done for all possibilities of adsorption modes!')
