#!/usr/bin/env python3

'''moldl: Molecule structure downloader
Written by Vikram Kashyap, Autumn 2020'''

#Imports
import pubchempy as pcp		#PubChem API
import mendeleev as mdv		#Mendeleev Periodic Table Library
import argparse
import sys
import os
import re
import json
import glob


#Convert MOL formated structure to XYZ
def convert_mol_to_xyz(mol_data, compound):
	lines = mol_data.split(b'\n')

	istart = 0
	for l in lines:
		if b'V2000' in l.split() or b'V3000' in l.split():
			break
		else:
			istart = istart + 1

	atomre = re.compile(b'\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\w+)\s+')
	atoms = []
	for l in lines[istart + 1:]:
		m = atomre.match(l)
		if m is not None:
			atoms.append(m.groups())
		else:
			break

	xyz_data = str(len(atoms)).encode() + b'\n' + str.encode(compound.synonyms[0]) + b',' \
		+ str.encode(compound.synonyms[1]) + b',' \
		+ str.encode(compound.synonyms[2]) + b'\n'
	atomcoords = [(a[-1] + b'\t' + a[0] + b'\t' + a[1] + b'\t' + a[2]) for a in atoms]
	xyz_data += b'\n'.join(atomcoords) + b'\n'

	return xyz_data


#Convert SDF formatted structure to MOL
def convert_sdf_to_mol(sdf_data, compound):
	linenum = 1
	i = 0
	while (linenum < 3):
		if (sdf_data[i] == ord('\n')):
			linenum+=1
		i+=1
	mol_data = sdf_data[:i] + str.encode(compound.synonyms[0]) + b',' \
		+ str.encode(compound.synonyms[1]) + b',' \
		+ str.encode(compound.synonyms[2]) + sdf_data[i:sdf_data.index(b'M  END')+len('M  END\n')]
	
	return mol_data


#Retrieve files
def retrieve_files(args):
	#Check element
	try:
		if (args.element.isnumeric()):
			esymb = mdv.element(int(args.element)).symbol
		else:
			esymb = mdv.element(args.element.capitalize()).symbol
	except:
		sys.stderr.write(sys.argv[0] + ': Unable to identify target element\n')
		sys.exit()

	#Check number of atoms
	if (args.numatoms == None):
		num_specified = False
	elif (args.numatoms < 1):
		sys.stderr.write(sys.argv[0]
			+ ": Invalid number of atoms of target element, numatoms must be positive integer")
		sys.exit()
	else:
		num_specified = True

	#Check max results
	if (args.maxresults < 1):
		sys.stderr.write(sys.argv[0] + \
			": Invalid number of results requested, maxresults must be positive integer")
		sys.exit()

	#Load database log
	try:
		with open(args.path+'/dblog.json', 'r') as f:
			dblog = json.load(f)
	except:
		print('Local database log not found. Initializing directory as new database.')
		dblog = {}

	#Retrieve search results from PubChem 
	search_results = pcp.get_compounds(esymb, 'smiles', searchtype='substructure', \
		MaxRecords=args.maxresults, record_type='3d')

	#Filter results
	refined_results = []
	for compound in search_results:
		if (num_specified and compound.elements.count(esymb) == args.numatoms
			and compound.charge == 0):
			refined_results.append(compound)

	#Convert to correct format and save based on format
	dlcount = 0
	for compound in refined_results:
		#Generate filepath of structure
		filepath = args.path+'/'+str(compound.cid)
		if (args.format.lower() == 'mol'):
			filepath += '.mol'
		elif (args.format.lower() == 'xyz'):
			filepath += '.xyz'
		else:
			filepath += '.sdf'
	
		#Check if already in local db (don't redownload)
		if os.path.isfile(filepath):
			continue

		#Retrieve SDF from PubChem
		file_data = pcp.get(compound.cid, output='SDF', record_type='3d')
		
		if file_data is None:
			continue
	
		#Apply any requested format conversions
		if (args.format.lower() == 'mol'):
			file_data = convert_sdf_to_mol(file_data, compound)
		elif (args.format.lower() == 'xyz'):
			file_data = convert_mol_to_xyz(convert_sdf_to_mol(file_data, compound), compound)

		#Write file
		with open(filepath, 'wb') as f:
			f.write(file_data)

		#Add 3 most common names with CID to database log
		dblog[compound.synonyms[0].lower()] = compound.cid
		dblog[compound.synonyms[1].lower()] = compound.cid
		dblog[compound.synonyms[2].lower()] = compound.cid

		dlcount+=1

	#Save log file
	with open(args.path+'/dblog.json', 'w') as f:
		json.dump(dblog, f)

	
#Lookup CID based on name, first searching local database, then online if not found
def lookup_cid(args):
	#Load database log
	try:
		with open(args.path+'/dblog.json', 'r') as f:
			dblog = json.load(f)
		if args.lookup.lower() in dblog:
			cid = dblog[args.lookup.lower()]
			print('CID ' + str(cid) + ': ' + ','.join(glob.glob(args.path+'/'+str(cid)+'.*')))
			return
		else:
			print('"'+args.lookup+'" not found in local database. Searching online...')
	except Exception as e:
		print(e)
		print('Local database log not found.')

	print('Name could correspond to following CIDs and local files (if present):')
	cids = pcp.get_cids(args.lookup, 'name', 'substance', list_return='flat')
	if len(cids) == 0: print('None')
	for cid in cids:
		 print(str(cid) + ': ' + ','.join(glob.glob(args.path+'/'+str(cid)+'.*')))
			
		

#Main function, parse and handle command line requests
def main():

	#Parse arguments
	argparser = argparse.ArgumentParser(description='moldl: Download molecular structure files', \
		formatter_class = argparse.RawDescriptionHelpFormatter,
		epilog='Example:\n\t moldl -e P -n 1 -m 5 --path molecules --format mol \n\n \
		Peace, and happy sciencing.')
	commandaction = argparser.add_mutually_exclusive_group(required=True)
	commandaction.add_argument('-e', '--element',
				metavar = 'Element',
				type=str,
				help='element to search for (name, symbol, or atomic number)')
	commandaction.add_argument('-l', '--lookup',
				metavar = 'Name',
				type=str,
				help='Compound name to search for in local database or online')
	argparser.add_argument('-n', '--numatoms',
				metavar = 'NumAtoms',
				type=int,
				help='number of atom of the target element to be present in molecule (default any)')
	argparser.add_argument('-m', '--maxresults',
				metavar = 'MaxResults',
				type=int,
				help='maximum number of results to download (will likely not achieve max)',
				default = 10)
	argparser.add_argument('-f', '--format',
				metavar = "Format",
				type=str,
				help='file format in which to save structure file (SDF, MOL, XYZ)',
				choices = ['sdf', 'SDF', 'mol', 'MOL', 'xyz', 'XYZ'],
				default = 'sdf')
	argparser.add_argument('-p', '--path',
				metavar = "Path",
				type=str,
				help='path to directory in which to download files (default current dir)',
				default = os.path.dirname(os.path.realpath(__file__)))

	args = argparser.parse_args()

	
	#Actions
	#Retrieve files
	if args.element is not None:
		retrieve_files(args)

	#Lookup molecules locally
	elif args.lookup is not None:
		lookup_cid(args)


#Run main if called as command
if __name__ == "__main__":
	main()
