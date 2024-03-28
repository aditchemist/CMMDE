# Computational Molecular and Material Design Environment (CMMDE)

## Background
This software is developed for decreasing the barrier in using popular open-source computational chemistry software. Currently, CMMDE can use the following software:

1. [Orca](https://orcaforum.kofo.mpg.de/app.php/portal)
2. [DFTB+](https://dftbplus.org/) 
3. [GROMACS](https://manual.gromacs.org/current/download.html)
4. [DOCK6](https://dock.compbio.ucsf.edu/DOCK_6/index.htm)
5. [DCDFTBMD](http://www.chem.waseda.ac.jp/dcdftbmd/?lang=en)
6. [Quantum Espresso](https://www.quantum-espresso.org/)
7. [xTB](https://github.com/grimme-lab/xtb)
8. [Open Babel](http://openbabel.org/wiki/Main_Page)

## About CMMDE
CMMDE is a set of tools based on Python for running computational jobs, as well as analyzing, visualizing, and post-processing the results in free/libre and open source applications for computational molecular & material design. 

- Core developers: Aditya Wibawa Sakti, Atthar Luqman Ivansyah, Hasan Al Rasyid, Muhamad Abdulkadir Martoprawiro, Aulia Sukma Hutama
- Contributors: Athiya Mahmud Hanna, Arifin, Daniel Sethio
- Core reviewers: Rahmat Gunawan, Imam Siswanto, Parsaoran Siahaan, Nova Pratiwi Indriyani
- Committed users: Yusthinus Thobias Male, Veliyana Londong, Mirella Fonda, Riyanto, Badra, Hilda, Rustaman, Edu

CMMDE telah diluncurkan dan akan diikuti oleh serangkaian lokakarya:

**Sabtu, 17 September 2022 (10:00 WIB):**</br>
Penyiapan Server Komputasi</br>
Tautan: [https://mki.ac/CMMD-server](https://mki.ac/CMMD-server)

**Sabtu, 24 September 2022 (10:00 WIB):**</br>
Workshop on Text-Based CMMDE</br>
Tautan:  [https://mki.ac/CMMDE-text]( https://mki.ac/CMMDE-text) 

**Sabtu, 1 Oktober 2022 (10:00 WIB):**</br>
Workshop on Web-Based CMMDE</br>
Tautan: [https://mki.ac/CMMDE-web]( https://mki.ac/CMMDE-web)

Informasi lebih lanjut: cmmde@mki.or.id

## Requirements
Sebelum menginstall CMMDE, perhatikan bahwa anda harus terlebih dahulu menginstall `slurm` untuk mengatur pengiriman pekerjaan komputasi anda ke suatu sistem antrian. Salah satu metode yang paling mudah untuk menginstall `slurm` dapat ditemukan di [sini](https://drtailor.medium.com/how-to-setup-slurm-on-ubuntu-20-04-for-single-node-work-scheduling-6cc909574365). 

## Installation
1. Cloning the repository:
```bash
git clone https://git.mki.or.id/coredev/cmmde
```
2. Change directory to CMMDE:
```bash
cd cmmde
```
3. Install CMMDE dengan mengetik:
```bash
./install.sh
```
Jika ditemukan _error_:
```bash
FileNotFoundError: [Errno 2] No such file or directory: '/usr/bin/pip3.9'
```
Lakukan: 
```bash
ln -s /usr/bin/pip /usr/bin/pip3.9
```
4. Jika ingin menjalankan program cmmde-gui:
```bash
cd cmmde_gui
panel serve gui.py --autoreload
```

## Catatan
Jika hanya ingin menginstall CMMDE versi _text-based_, cukup lakukan langkah 1-3. Langkah 4 hanya dilakukan jika ingin menjalankan CMMDE versi _graphical user interface_. Namun versi ini perlu menginstall `panel-chemistry versi 0.0.11`.
```bash
 pip install -Iv panel-chemistry==0.0.11
```
