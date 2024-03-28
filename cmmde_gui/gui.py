#!/usr/bin/env python3
# pylint: disable=missing-function-docstring
from __future__ import print_function
import panel as pn

from panel_chemistry.widgets import JSMEEditor
from panel.interact import interact
import os
import sys
from panel_chemistry.pane import NGLViewer
from panel_chemistry.pane.ngl_viewer import EXTENSIONS
import py3Dmol
from panel_chemistry.pane import Py3DMol
import subprocess
import uuid
# from iodata import IOData


# def test_can_construct():
#     JSMEEditor()


def cmmde_gui():
    geom = ''
    pn.extension("jsme", sizing_mode="stretch_width")
    
    editor = JSMEEditor(value = " ",height=500,format="smiles",subscriptions=['smiles'])
    
     
# Name your molecule
    Molecule_input = pn.widgets.TextInput(name="Molecule name")
    id_input = pn.widgets.TextInput(name="Input your name")
    charge = pn.widgets.TextInput(name="Charge",value="0")
    mult = pn.widgets.TextInput(name="Spin multiplicity",value="1")
    workdir = os.getenv("HOME") + "/" + "scr"
# Terminal widget
    TextArea= pn.widgets.TextAreaInput(value = "Computational Molecular and Material Design Environment\n Authors:\n Universitas Pertamina\n Institut Teknologi Sumatera\n Institut Teknologi Bandung\n Masyarakat Komputasi Indonesia\n\nSupported by:\n Konsorsium Pengembangan Sains Komputasi\n",
        height = 500, disabled=True
        )
# CMMDE software options
    software_main = pn.widgets.Select(name="Software selections",value='Orca',options=['Orca','XTB','Dcdftbmd','Quantum Espresso','GROMACS'])
    software = {'Orca':'orca','GROMACS':'gromacs','Dcdftbmd':'dcdftb','Quantum Espresso':'qe','XTB':'xtb'}

# Orca card
    job_orca = pn.widgets.Select(name="Job selections",value=['Single point calculation'],options=['Single point calculation','Geometry optimization','Frequency calculation','TS optimizer','Nudged elastic band'])
    job_orca_dict = {'Single point calculation':'sp','Geometry optimization':'opt','Frequency calculation':'freq','TS optimizer':'ts','Nudged elastic band':'neb'} 
    method_orca = pn.widgets.Select(name="Method selections",value='GFN2-xTB',options=['GFN2-xTB','GFN1-xTB','B3LYP/def2-svp','M06/def2-svp'])
    method_orca_dict = {'GFN2-xTB':'XTB2', 'GFN1-xTB':'XTB1','DFTB2':'DFTB2','DFTB2-gammah':'DFTB2_gammah','DFTB3':'DFTB3','DFTB3-diag':'DFTB3-diag','B3LYP/def2-svp':'B3LYP def2-svp','M06/def2-svp':'M06 def2-svp'}
    dispersion_cor = pn.widgets.Select(name="Dispersion corrections",value='None',options=['None','D3','D3BJ','D4'])
    solvent = pn.widgets.Select(name='Solvent',value='None',options=['None','water','acetone','acetonitrile','aniline','benzaldehyde','benzene','CH2Cl2','CHCl3','CS2','dioxane','DMF','ethanol','ether','ethylacetate','furane','hexadecane','hexane','octanol','octanol(wet)','phenol','toluene','THF'])
    run_orca_btn = pn.widgets.Button(name="Run Orca!",button_type='primary')
    def run_orca(event):
        # RunMessage.value = " "
        # terminal.clear()
        # unik = str(uuid.uuid4().hex)
        
        Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value 
        
        os.chdir(Folder)

        
        job_list = [job[i] for i in job_orca.value] 
        jobs = ",".join(job_list)
        
        if editor.value == "":
            FileInput.save("geom.xyz")
            geom = "geom.xyz"
        else:
            geom = editor.value
            
        TextArea.value = TextArea.value + "\n" + "Mempersiapkan Struktur 3 Dimensi!"
        cmd = subprocess.run(["cmmde.py","-i","{}".format(geom),"-s","orca","-j","{}".format(jobs),"-m","{}".format(method[method_orca.value]),"-c","{}".format(charge.value),"-mult","{}".format(mult.value)],capture_output=True,text=True)
        #terminal.subprocess.run("cmmde.py","-i{}".format(editor.value),"-s{}".format(software[software_main.value]), "-j{}".format(jobs), "-m{}".format(method[method_btn.value]))
        TextArea.value = TextArea.value + "\n" + "Perhitungan anda telah tersubmit!"
    # RunMessage = pn.widgets.StaticText()
    run_orca_btn.on_click(run_orca)

# CMMDE job options
    job_main = pn.widgets.Select(name="Job selections",value='Single point calculation',options=['Single point calculation','Geometry optimization','Frequency calculation'])
    job = {'Single point calculation':'sp','Geometry optimization':'opt','Frequency calculation':'freq'}


# CMMDE method options
    method_btn = pn.widgets.Select(name="Method selections",value='GFN2-xTB',options=['GFN2-xTB','GFN1-xTB','B3LYP/def2-svp','M06/def2-svp','DFTB2','DFTB2-gammah','DFTB3','DFTB3-diag'])
    method = {'GFN2-xTB':'XTB2', 'GFN1-xTB':'XTB1','DFTB2':'DFTB2','DFTB2-gammah':'DFTB2_gammah','DFTB3':'DFTB3','DFTB3-diag':'DFTB3-diag','B3LYP/def2-svp':'B3LYP def2-svp','M06/def2-svp':'M06 def2-svp'}


    

# File input (if you don't want to draw the structure)
    FileInput = pn.widgets.FileInput(title='Input structure',accept='.xyz,.vasp,.pdb')
    # def fileinput(event):
    #     if FileInput.value is not None:
    #         FileInput.save("geom.xyz")
    # FileInput.param.watch(fileinput,'value')

# CMMDE running button
    Run_btn = pn.widgets.Button(name="Run CMMDE!",button_type='primary')
    # RunMessage = pn.widgets.StaticText()
    
    def run(event):
        # RunMessage.value = " "
        # terminal.clear()
        # unik = str(uuid.uuid4().hex)
        
        # job_list = [job[i] for i in job_main.value] 
        # jobs = ",".join(job_list)
        Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value + "/" + job[job_main.value] 
        if not os.path.exists(Folder):
            os.makedirs(Folder)
        if os.path.exists("{}/geom.smi".format(Folder)):
            os.system("rm {}/geom.smi".format(Folder))
        os.chdir(Folder)

        
        if editor.value == "":
            FileInput.save("geom.xyz")
            geom = "geom.xyz"
        else:
            geom = editor.value
            
        TextArea.value = TextArea.value + "\n" + "Mempersiapkan Struktur 3 Dimensi!"
        list_commands = ["cmmde.py","-i","{}".format(geom),"-s","{}".format(software[software_main.value]),"-j","{}".format(job[job_main.value]),"-m","{}".format(method[method_btn.value]),"-c","{}".format(charge.value),"-mult","{}".format(mult.value)]
        if dispersion_cor.value != 'None':
            new_commands = ["-disp","{}".format(dispersion_cor.value)]
            for i in new_commands:  
                list_commands.append(i)
        if solvent.value != 'None':
            new_commands = ["-solvent","{}".format(solvent.value)]
            for i in new_commands:
                list_commands.append(i)  
        cmd = subprocess.run(list_commands,capture_output=True,text=True)
        #terminal.subprocess.run("cmmde.py","-i{}".format(editor.value),"-s{}".format(software[software_main.value]), "-j{}".format(jobs), "-m{}".format(method[method_btn.value]))
        TextArea.value = TextArea.value + "\n" + "Perhitungan anda telah tersubmit!"

        # RunMessage.value = "Perhitungan telah tersubmit!"

    download_opt = pn.widgets.FileDownload(file="cmmd.xyz",filename="optimized.xyz")

    Run_btn.on_click(run)

# Check directory button
    checkdir_btn = pn.widgets.Button(name="Generate work directory",type="primary")
    TextWarning = pn.widgets.StaticText()
    def checkdir(event):
        Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value 
        isExist = os.path.exists(Folder)
        TextWarning.value = ""
        if isExist:
            TextWarning.value = "Directory exists! Change the molecule name!"
        else:
            os.makedirs(Folder)
            TextWarning.value = "Successfully create the directory!"

    checkdir_btn.on_click(checkdir)






# Post calculations
    post_calc = {'Frequency calculation':'freq', 'Radial distribution function':'rdf','Mean Square Displacement':'msd','Time-dependent calculation':'td','Thermochemistry calculation':'thermo','Optimized energy':'opt','IR plot':'ir'}
    post_btn = pn.widgets.Select(name="Job Selection",value='Frequency calculation', options=['Frequency calculation','Radial distribution function', 'Mean Square Displacement','Time-dependent calculation','Thermochemistry calculation','Optimized energy','IR plot'])
# Post Calculation CMMDE software options
    post_software_main = pn.widgets.Select(name="Software selections for post calculations",value='Orca',options=['Orca','Dcdftbmd','Quantum Espresso'])
    post_software = {'Orca':'orca','Dcdftbmd':'dcdftb','Quantum Espresso':'qe'}
    # Post CMMDE method options
    post_method_btn = pn.widgets.Select(name="Method selections for post calculations",value='GFN2-xTB',options=['GFN2-xTB','GFN1-xTB','DFTB2','DFTB2-gammah','DFTB3','DFTB3-diag','B3LYP/def2-svp'])
    post_method = {'GFN2-xTB':'XTB2', 'GFN1-xTB':'XTB1','DFTB2':'DFTB2','DFTB2-gammah':'DFTB2_gammah','DFTB3':'DFTB3','DFTB3-diag':'DFTB3-diag','B3LYP/def2-svp':'B3LYP def2-svp'}
    def post_calculation(event):
        Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value
        # terminal.clear() 
        if post_calc[post_btn.value] == 'thermo' and post_software[post_software_main.value] == 'orca':
            # Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value 
            # os.chdir(Folder)
            Post_folder = Folder + "/" + post_calc[post_btn.value] 
            if not os.path.exists(Post_folder):
                os.makedirs(Post_folder)
            os.chdir(Post_folder)
            os.system("cp {}/freq/cmmd.out .".format(Folder))
            cmd = subprocess.run(["cmmdepost.py","-j","{}".format(post_calc[post_btn.value]),"-s","{}".format(post_software[post_software_main.value])],capture_output=True,text=True)
            TextArea.value = TextArea.value + "\n" + cmd.stdout
        elif post_calc[post_btn.value] == 'ir' and post_software[post_software_main.value] == 'orca':
            # Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value 
            # os.chdir(Folder) 
            Post_folder = Folder + "/" + post_calc[post_btn.value] 
            if not os.path.exists(Post_folder):
                os.makedirs(Post_folder)
            os.chdir(Post_folder)
            os.system("cp {}/freq/cmmd.out .".format(Folder))
            cmd = subprocess.run(["cmmdepost.py","-j","{}".format(post_calc[post_btn.value]),"-s","{}".format(post_software[post_software_main.value])])
            TextArea.value = TextArea.value + "\n" + "Plot spektrum IR berhasil dilakukan"
            # terminal.subprocess.run("cmmdepost.py","-j{}".format(post_calc[post_btn.value]),"-s{}".format(post_software[post_software_main.value]))
        elif post_calc[post_btn.value] == 'opt' and post_software[post_software_main.value] == 'orca':
            # Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value 
            # os.chdir(Folder)
            os.chdir(Folder)
            cmd = subprocess.run(["cmmdepost.py","-j","{}".format(post_calc[post_btn.value]),"-s","{}".format(post_software[post_software_main.value]),"-m","{}".format(post_method[post_method_btn.value])],capture_output=True,text=True)
            TextArea.value = TextArea.value + "\n" + cmd.stdout  
        else:
            # Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value 
            # os.chdir(Folder) 
            Post_folder = Folder + "/" + post_calc[post_btn.value] 
            if not os.path.exists(Post_folder):
                os.makedirs(Post_folder)
            os.chdir(Post_folder)
            list_commands = ["cmmde.py","-i","{}/cmmd.xyz".format(Folder+"/"+job[job_main.value]),"-s","{}".format(post_software[post_software_main.value]),"-j","{}".format(post_calc[post_btn.value]),"-m","{}".format(post_method[post_method_btn.value])]
            if dispersion_cor.value != 'None':
                new_commands = ["-disp","{}".format(dispersion_cor.value)]
                for i in new_commands:  
                    list_commands.append(i)
            if solvent.value != 'None':
                new_commands = ["-solvent","{}".format(solvent.value)]
                for i in new_commands:
                    list_commands.append(i)
            cmd = subprocess.run(list_commands,capture_output=True,text=True)
            # terminal.subprocess.run("cmmde.py","-i{}".format("../cmmd.xyz"),"-s{}".format(post_software[post_software_main.value]), "-j{}".format(post_calc[post_btn.value]), "-m{}".format(post_method[post_method_btn.value]))
            TextArea.value = TextArea.value + "\n" + cmd.stdout
    runpost_btn = pn.widgets.Button(name="Run post calculation!",button_type='primary')
    runpost_btn.on_click(post_calculation)
# Check the queue progress
    def progress(event):
        # terminal.clear()
        # terminal.subprocess.run("squeue")
        cmd = subprocess.run(["squeue"],capture_output=True,text=True)
        TextArea.value = TextArea.value + "\n" + cmd.stdout


    Progress_btn = pn.widgets.Button(name="Check queue",button_type='primary')    
    Progress_btn.on_click(progress)

# Check the calculation progress
    def calc_progress(event):
        Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value + "/" + job[job_main.value] 
        os.chdir(Folder)
        cmd = subprocess.run(["tail", "-n", "10", "cmmd.out"],capture_output=True,text=True)
        TextArea.value = TextArea.value + "\n" + cmd.stdout 
    Checkcalc_btn_main = pn.widgets.Button(name="Check calculation",button_type='primary')
    Checkcalc_btn_main.on_click(calc_progress)

    def calc_progress_post(event):
        Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value + "/" + post_calc[post_btn.value]
        os.chdir(Folder)
        cmd = subprocess.run(["tail", "-n", "10", "cmmd.out"],capture_output=True,text=True)
        TextArea.value = TextArea.value + "\n" + cmd.stdout 
    Checkcalc_btn_post = pn.widgets.Button(name="Check post-calculation",button_type='primary')
    Checkcalc_btn_post.on_click(calc_progress_post)
# Slab Builder
    hkl_input = pn.widgets.TextInput(name="Miller index (hkl)",placeholder="Example: 100")
    size_input = pn.widgets.TextInput(name="Dimension",placeholder="Example: 2x2")
    layer_input = pn.widgets.TextInput(name="Layer",placeholder="Example: 2")
    slabbuilder_btn = pn.widgets.Button(name="Build it!",button_type="primary")
    Material_input = pn.widgets.TextInput(name="Material name")
    
    Material_upload = pn.widgets.FileInput(title='Input structure')
    # def materialinput(event):
    #     if Material_upload.value is not None:
    #         Material_upload.save("POSCAR")
    # Material_upload.param.watch(materialinput,'value')

    # Generate material folder button
    materialdir_btn = pn.widgets.Button(name="Generate work directory",type="primary")
    
# # Download Button
    
    download_slab = pn.widgets.FileDownload(file="slab.vasp",filename="slab.vasp")
    download_slab_xyz = pn.widgets.FileDownload(file="cmmd.xyz",filename="slab.xyz")
    download_spec_plot = pn.widgets.FileDownload(file="IR.pdf",filename="IR.pdf")
    download_spec_raw = pn.widgets.FileDownload(file="IR_fit.dat".format(workdir+"/"+id_input.value+"/"+Molecule_input.value+"/"+"ir"),filename="IR_fit.dat")
    download_opt_plot = pn.widgets.FileDownload(file="optimized.pdf",filename="optimized.pdf")
    download_opt_raw = pn.widgets.FileDownload(file="optimized.dat",filename="optimized.dat")
    def materialgen(event):
        Folder = workdir + "/" + id_input.value + "/" + Material_input.value 
        isExist = os.path.exists(Folder)
        TextWarning.value = ""
        if isExist:
            TextWarning.value = "Directory exists! Change the material name!"
        else:
            os.makedirs(Folder)
            TextWarning.value = "Successfully create the directory!"

    materialdir_btn.on_click(materialgen)
    
    def slab_builder(event):
        Folder = workdir + "/" + id_input.value + "/" + Material_input.value 
        os.chdir(Folder)
        Material_upload.save("POSCAR")
        # os.system("mv geom.xyz POSCAR")
        hkl = hkl_input.value
        size = size_input.value
        layer = layer_input.value
        cmd = subprocess.run(["cmmdepre.py","-j","surface","-hkl","{}".format(hkl),"-s","{}".format(size),"-n","{}".format(layer),"-i","POSCAR"])
        os.system("mv slab_{}.xyz cmmd.xyz".format(hkl))
        os.system("mv slab_{}.vasp slab.vasp".format(hkl))
        cmd = subprocess.run(["echo","({})-surface construction done!".format(hkl)],capture_output=True,text=True)
        TextArea.value = TextArea.value + "\n" + cmd.stdout
        # xyzview = py3Dmol.view()
        
        
        # with open('slab_{}.xyz'.format(hkl),'r') as f:
        #     xyz = f.read()
        # xyzview.addModel(xyz,'xyz')
        # xyzview.setStyle('stick')
        # xyzviewer.object = xyzview
    slabbuilder_btn.on_click(slab_builder)
# Solution Builder
    filename_solute = pn.widgets.TextInput(title="Solute name")
    filename_solvent = pn.widgets.TextInput(title="Solvent name")
    solute_upload = pn.widgets.FileInput(accept=".xyz",multiple=True)
    solvent_upload = pn.widgets.FileInput(accept=".xyz")
    def SolutionBuilder(event):
        solute = filename_solute.value.split(",")
        for i in solute:
            solute_folder = workdir + "/" + id_input.value + "/" + "{}".format(i)
            os.makedirs(i)
            os.chdir(i)
            solute_upload.save(solute_upload.filename)
        solvent = filename_solvent.value
        solvent_folder = workdir + "/" + id_input.value + "/" + solvent

        
        
        os.chdir("../{}".format(solvent_folder))
        solvent_upload.save(solvent_upload.filename)
        os.chdir("../")
        cmd = subprocess.run(["cmmde.py","-s","gromacs","-mt","{}".format(solute_folder)])

        solvent = filename_solvent.value
        solvent_upload.save(solvent_upload.filename)




# Visualize the results
    xyzview = py3Dmol.view()
    
    xyzviewer = Py3DMol(xyzview, height=400, sizing_mode="stretch_width",name="CMMDE viewer")
    def visualize(event):
        xyzview = py3Dmol.view()
        
        if Material_input.value == "":
            Folder = workdir + "/" + id_input.value + "/" + Molecule_input.value + "/" + job[job_main.value] 
        else:
            Folder = workdir + "/" + id_input.value + "/" + Material_input.value
        os.chdir(Folder)
        with open('cmmd.xyz','r') as f:
            xyz = f.read()
        xyzview.addModel(xyz,'xyz')
        xyzview.setStyle({'stick':{},'sphere':{'scale':.30}},viewer=(0,0))
        xyzview.zoomTo()
        xyzviewer.object = xyzview
        

    visual_btn = pn.widgets.Button(name="Visualize!", button_type='primary')
    visual_btn.on_click(visualize)

    def set_background(color='0xeeeeee'):
        xyzview.setBackgroundColor(color)
        xyzviewer.param.trigger("object")
    set_background("#e6f6ff")

    accent = "#0072B5"

    # background = pn.widgets.ColorPicker(value="#e6f6ff", name="Background")
    # pn.bind(set_background, color=background, watch=True)     
    # def set_style(style="stick"):
    #     xyzview.setStyle({style: {}})
    #     xyzview.zoomTo()
    #     xyzviewer.param.trigger("object")
    
    # set_style("stick")
    # style=pn.widgets.RadioButtonGroup(value="stick", options=["stick", "sphere"], name="Style", button_type="success")
    # set_style=pn.bind(set_style, style=style, watch=True) 
    
       
#### Wrap them all together ##########       
    
    return pn.template.MaterialTemplate(
        sidebar_width=410,
        site="Computational Molecular and Material Design Environment",
        title="CMMDE",
        main=[TextArea, editor, xyzviewer],
        sidebar=[pn.Card(id_input,title="User Information",collapsed=True),pn.Card(Molecule_input,charge,mult,checkdir_btn,TextWarning,pn.Card(FileInput,title="Upload molecule",collapsed=True),title="Molecule Information",collapsed=True),pn.Card(Material_input, materialdir_btn,TextWarning, pn.Card(Material_upload,title="Unit cell",collapsed=True),hkl_input,size_input,layer_input,pn.Row(slabbuilder_btn,visual_btn),pn.Row(download_slab_xyz,download_slab),title="Surface Builder",collapsed=True),pn.Card(software_main,job_main,method_btn, dispersion_cor,solvent, pn.Row(Run_btn,Progress_btn),pn.Row(Checkcalc_btn_main,visual_btn),download_opt,title="Main Calculation",collapsed=True),pn.Card(post_software_main,post_btn,post_method_btn,dispersion_cor,solvent,pn.Row(runpost_btn,Progress_btn),pn.Row(Checkcalc_btn_post),pn.Row(download_spec_raw,download_spec_plot),pn.Row(download_opt_raw,download_opt_plot),title="Post-Calculation",collapsed=True)],
        header_background=accent, accent_base_color=accent
    )    


if __name__.startswith("bokeh"):
    cmmde_gui().servable()


