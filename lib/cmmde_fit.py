import sys
import argparse
import matplotlib
matplotlib.use('PS')
#matplotlib.rc('text', usetex=True)
import numpy as np
import statsmodels.api as sm
#matplotlib.use('Qt5agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Compute Diffusion Coeffients from MSD data (nxy format)')
parser.add_argument('data', metavar='data', type=str)
parser.add_argument('-n', '--noheader', action="store_true")
parser.add_argument('-s', '--start', type=int, default=0, help='StartStep')
parser.add_argument('-e', '--end', type=int, default=-1, help='EndStep')
opt = parser.parse_args(sys.argv[1:])


data_file = opt.data

x = []
ys = []
ncol = -1
lineno = 0
headers = []
with open(data_file, 'r') as f:
    for line in f:
        lineno +=1
        if (line.startswith('#')):
            if opt.noheader:
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


print('There are {} columns of data.'.format(ncol))
print()

print(headers)
if (len(headers) != ncol):
    headers = []
    for col in range(ncol):
        headers.append('Column {}'.format(col+1))

for col in range(ncol):
    print ('Column {}:'.format(col+1))

    all_data_x = np.array(x, dtype=np.double).reshape(-1, 1)
    data_x = all_data_x[opt.start:opt.end]

    all_data_y = [d[col] for d in ys]
    data_y = np.array(all_data_y[opt.start:opt.end], dtype=np.double).reshape(-1, 1)
    X = sm.add_constant(data_x)
    ols = sm.OLS(data_y, X)
    ols_result = ols.fit()

    summ = ols_result.summary()
    print(summ)
    print('*'*80)
    diff_coeff = ols_result.params[1]/6.0
    print(' Diffusion coefficient (1/6): {:16.10F} A^2/ps, {:17.12E} m^2/s'.format(diff_coeff, diff_coeff*1.0e-8))
    diff_coeff = ols_result.params[1]/2.0
    print(' Diffusion coefficient (1/2): {:16.10F} A^2/ps, {:17.12E} m^2/s'.format(diff_coeff, diff_coeff*1.0e-8))
    diff_coeff = ols_result.params[1]/4.0
    print(' Diffusion coefficient (1/4): {:16.10F} A^2/ps, {:17.12E} m^2/s'.format(diff_coeff, diff_coeff*1.0e-8))
    print('*'*80)
    print()

    predicted_y = ols_result.predict(X)

    plt.plot(x, all_data_y, label=headers[col])
    plt.plot(data_x, predicted_y, 'k--')

plt.legend(loc=2)
plt.xlabel('Time (ps)')
plt.ylabel('MSD (Angstrom$^2$)')
plt.savefig('msd.pdf', dpi=1200, format='pdf')
# plt.show()
