import os
import panel as pn
import subprocess
def path(workdir,id_input,textarea):
	reaktan_input = pn.widgets.FileInput(name="Struktur reaktan",accept='.xyz')
	produk_input = pn.widgets.FileInput(name="Struktur produk", accept='.xyz')
	path_name = pn.widgets.TextInput(name="Directory name")
	alp = pn.widgets.TextInput(name="Bias energy",value="1.2")
	kpush = pn.widgets.TextInput(name="Push force",value="0.003")
	kpull = pn.widgets.TextInput(name="Pull force",value="-0.015")
	ppull = pn.widgets.TextInput(name="Histogram bin",value="0.05")
	def run_path(event):
		reaktan_input.save("reaktan.xyz")
		reaktan = "reaktan.xyz"
		produk_input.save("produk.xyz")
		produk = "produk.xyz"
		Folder = workdir + "/" + id_input + "/" + path_name
		if os.path.exists(Folder):
			textarea = textarea + "Directory exists! Please change!"
		else:
			os.makedirs(Folder)
		os.chdir(Folder)
		cmd = subprocess.run(["cmmde.py","-s","xtb","-i","{}".format(reaktan),"-j","path","-produk","{}".format(produk),"-alp","{}".format(alp.value),"-kpush",])

