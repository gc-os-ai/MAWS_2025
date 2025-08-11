import subprocess
import os
from openmm import app

## makeLib
# parametrizes a pdb or mol2 file to generate a .lib library file for tleap / nab
def makeLib(file_path, residue_name, connect0=None, connect1=None, charges='bcc', atom_type='gaff', force_field_aptamer='leaprc.RNA.OL3', force_field_ligand='leaprc.protein.ff19SB', parameterized=False, conda_env="maws_p3"):
	name, extension = file_path.split('/')[-1].split(".")
	lib_path = '/'.join(file_path.split('/')[:-1]+[residue_name])
	tleap_input = """
	source %s
	source %s
	%s = loadmol2 %s.mol2
	loadamberparams %s.frcmod"""%(force_field_aptamer,
		 force_field_ligand,
		 residue_name, name, 
		 lib_path) 
	if connect0 and connect1:
		tleap_input += """
	set %s head %s.1.%s
	set %s tail %s.1.%s
	set %s.1 connect0 %s.head
	set %s.1 connect1 %s.tail"""%(residue_name, residue_name, connect0, 
					residue_name, residue_name, connect1, 
					residue_name, residue_name, 
					residue_name, residue_name)
	tleap_input += """
	check %s
        saveoff %s %s.lib
	savepdb %s %s_tmp.pdb
	quit
	"""%(residue_name, residue_name, lib_path, residue_name, lib_path)
	if parameterized:
		tleap_input = """
		source %s
		source %s
		%s = loadpdb %s.pdb
		saveoff %s %s.lib
		savepdb %s %s_tmp.pdb
		quit
		"""%(force_field_aptamer, force_field_ligand, residue_name, name, residue_name, lib_path,
			residue_name, lib_path)
	#Write LEaP infile
	with open("%s.in"%name,"w") as fil:
		fil.write(tleap_input)
	if not parameterized:
		#Execute
		subprocess.call("conda run -n %s && antechamber -i %s -fi %s -o %s.mol2 -fo mol2 -c %s -rn %s -at %s"%(conda_env,
													 								   file_path,
																					   extension, 
																					   name, 
																					   charges, 
																					   residue_name, 
																					   atom_type), 
						shell=True)
		subprocess.call("conda run -n %s && parmchk2 -i %s.mol2 -f mol2 -o %s.frcmod"%(conda_env, name, lib_path), shell=True)
	subprocess.call("conda run -n %s && tleap -f %s.in"%(conda_env, name), shell=True)
	PDB = app.PDBFile(lib_path + "_tmp.pdb")
	length = sum([1 for atom in PDB.topology.atoms()])
	#Cleanup
	if not parameterized:
		os.remove(name+".mol2")
	#os.remove(name+".frcmod")
	os.remove("%s.in"%name)
	os.remove("%s_tmp.pdb"%(lib_path)) # lib_path = PDB
	os.remove("leap.log")
	return length

## toggleHydrogens
# Toggles presence/abscence of Hydrogens as indicated by 
# @param boolean
# true: hydrogens or false: hydrogens
def toggleHydrogens(path, boolean=True, conda_env="maws_p3"):
	subprocess.call("conda run -n %s && reduce -Trim %s > %s"%(conda_env, path, path), shell=True)
	if boolean:
		subprocess.call("conda run -n %s && reduce -Build %s > %s"%(conda_env, path, path), shell=True)


