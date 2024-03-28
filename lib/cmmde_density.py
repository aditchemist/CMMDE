from scipy.constants import N_A

def density(geom):
	a1 = 0
	a2 = 0
	a3 = 0 
	with open(geom, 'r') as f:
		for line in f:
			if "CRYST1" in line:
				arr = line.split()
				a1 = float(arr[1])
				a2 = float(arr[2])
				a3 = float(arr[3])		
	V = a1*a2*a3*10**-24 # Volume dalam cm^3
	dens = mass/V # Massa jenis dalam g/cm^3

	print("Density = {}".format(dens))