#!/usr/bin/env python3

'''moldl: Molecular Structure Downloader
Written by Vikram Kashyap, starting Autumn 2020'''

#Imports
import pubchempy as pcp		#PubChem API
import mendeleev as mdv		#Mendeleev Periodic Table Library
import argparse
import sys
import os
import re
import json
import glob
import random
from pathlib import Path

##########################################
# Utility Functions
##########################################

def getatomicnum(element):
	'''
	Get atomic number of element

		Parameters:
			element (int or str): atomic number or symbol of element

		Returns:
			atomic number of element (int)
	'''
	TABLE =
		['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg',
		'al', 'si', 'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr',
		'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
		'rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd',
		'in', 'sn', 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd',
		'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf',
		'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po',
		'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm',
		'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs',
		'mt', 'ds ', 'rg ', 'cn ', 'nh', 'fl', 'mc', 'lv', 'ts', 'og']
	# If argument is already int, assume it is atomic num and pass back
	if type(element) is int:
		return element
	# If string, find in list of element symbols
	elif type(element) is str:
		try:
			return TABLE.index(str.lower())
		except:
			raise ValueError("Invalid element symbol string")
	# Perhaps in future allow identification by other types?
	else:
		raise TypeError( \
			"Elements can only be identified by atomic number of symbol")

##########################################
# Core Classes
##########################################

class Atom:
	'''
	Represents an individual atom in a molecule.
	This is a mostly hollow object because all actual data is stored in the related Molecule.
	Atom objects should only be created by a Molecule object

		Attributes:
			aid (int): atom ID unique within molecule
			element (int): atomic number of atom
			bonds ([Bond]): list of Bonds that Atom is part of
	'''

	def __init__(self, molecule, aid):
		'''
		Define atom by the molecule it is contained in and its atom ID within molecule

			Parameters:
				molecule (Molecule): Molecule object that contains atom
				aid (int): Atom ID within molecule

			Returns:
				Object representing atom in molecule (Atom)
		'''
		self.molecule = molecule
		self.aid = aid
		self.element = molecule.getelement(self.aid)

	def __getattr__(self, key):
		'''
		Get attributes of Atom that should be generated on demand (bonds)

			Paramters:
				key (str): name of attribute

			Returns:
				attribute corresponding to key
		'''
		if key == 'bonds':
			self.bonds = self._genBonds()
			return self.bonds
		else:
			raise AttributeError('Invalid Bond object attribute')

	def _genBonds(self):
		'''
		Get a list of Bond objects representing the bonds this atom has
		'''

		bids = self.molecule.getBonds(self.molecule.getbondsbyatom(aid=self.aid))


class Bond:
	'''
	Represents a bond between Atoms
	This is a mostly hollow object because all actual data is stored in the related Molecule
	Bond objects should only be created by a Molecule or Atom object
	'''
	def __init__(self, molecule, aid, bid):
		self.molecule = molecule
		self.bid = bid
		self.record = molecule.struct['bonds'][bid]
		self.order = molecule.struct['bonds'][bid]['order']


	def partnerof(self, atom):
		if self.record['atoms'][0] == atom.aid:
			return Atom(self.molecule, self.record['atoms'][1])
		elif self.record['atoms'][1] == atom.aid:
			return Atom(self.molecule, self.record['atoms'][0])
		else:
			return ValueError('Passed Atom is not part of Bond')

	@property
	def order(self):
		return self.record['order']

	@order.setter
	def order(self, value):
		self.record

	def __getattr__(self, key):
		if key == 'order':


class Molecule:
	'''
	Represents a certain molecule and stores its properties as attributes

		Attributes:
			cid (int): PubChem CID of molecule
			id ({}): dictionary of identifiers
			struct ({}): dictionary of structural information
			chem ({}): dictionary of chemical information
	'''

	def __init__(self, cid, properties=None, download=True):
		'''
		Define molecule by CID

			Parameters:
				cid (int): PubChem CID
				properties ({}): properly JMF formatted dictionary of properties
				download (bool): whether to download from PubChem
						if no properties passed (only set False if you know
						what you're doing)

			Returns:
				Molecule object
		'''
		self.cid = cid
		if properties is not None:
			self.__dict__ = properties
		elif download:
			self.download()

	def download(self):
		'''
		Load properties of molecule from PubChem online database

			Parameters:
				None

			Returns:
				True if properly loaded, else False
		'''
		# Fetch record from PubChem in JSON format, return False if record not found
		record = pcp.get_json(self.cid, record_type='3d')['PC_Compound'][0]
		if record is None: return False
		self.pcrecord = record

		# Extract properties from PubChem record
		# Hardcoded since PubChem record system is a mess. Hopefully they don't change their format

		# Top level layout
		self.id = {}	# Identifying information
		self.struct = {}	# Structure information
		self.chem = {}	# Chemical information

		# Import ID/Aliases of molecules
		self.id['cid'] = self.cid
		#TODO other ids

		# Import structure
		coords = record['coords'][0]['conformers'][0]
		self.struct['atoms'] = []
		for e, x, y, z in zip(record['atoms']['element'], coords['x'], coords['y'], coords['z']):
			atom = {'element': e,
				'coord': (x, y, z)}
			self.struct['atoms'].append(atom)
		self.struct['bonds'] = []
		for a1, a2, order in zip(record['bonds']['aid1'], record['bonds']['aid2'], record['bonds']['order']):
			bond = {'atoms': (a1, a2),
				'order': order}
			self.struct['bonds'].append(bond)

		# Derive structure properties
		self._derivestructproperties()

		# TODO Import chemical properties

		# Set attributes
		self.atoms = list(range(len(self.struct['atoms'])))
		self.bonds = list(range(len(self.struct['bonds'])))
		self.elements = self.struct['deriv']['elements']

	def getAtom(self, aid):
		return Atom(self, aid)


	def getAtoms(self, element=None):
		return [self.getAtom(self, aid) for aid in self.atomsbyelement(element)]


	def getBond(self, bid, aid):
		return Bond(self, aid, bid)


	def getBonds(self, group):
		if isinstance(group, Atom):


	def getelement(self, aid):
		'''
		Get atomic number of atom by Atom ID

			Parameters:
				aid (int): Atom ID within molecule

			Returns:
				atomic number of atom (int)
		'''
		return self.struct['atoms'][aid]['element']


	def getbondorder(self, bid):
		'''
		Get bond order of bond by Bond ID

			Parameters:
				bid (int): Bond ID within molecule

			Returns:
				order of bond (int)
		'''
		return self.struct['bonds'][bid]['order']


	def getbondpartner(self, bid, aid):
		b = _getbond(bid)
		if b[0] == aid: return b[1]
		if b[1] == aid: return b[0]
		else: return None


	def atomsbyelement(self, element):
		'''
		Get list of atom IDs by element

			Parameters:
				element (int or str): atomic number or symbol of element

			Returns:
				list of atom IDs for atoms in molecules of element
		'''
		return self.struct['deriv']['aidsbyelement'][getatomicnum(element)]

	def bondsbyatom(self, aid):
		'''
		Get list of bond IDs by atom ID

			Parameters:
				aid (int): atom ID to get bonds for

			Returns:
				list of bond IDs for bonds that involve given aid ([int])
		'''
		return self.struct['deriv']['bidsbyaid']


	def _getbond(self, bid):
		'''
		Get a reference to the dictionary for a bond

			Parameters:
				bid (int): Bond ID within molecule

			Returns:
				reference to dictionary within molecule that represents bond
		'''
		return self.struct['bonds'][bid]


	def _derivestructproperties(self):
		'''
		Calculate properties from core information
		Populates Molecule.struct:'deriv'

			Parameters:
				None
			Returns:
				None (int)
		'''

		# Make derived properties dictionary
		self.struct['deriv'] = {}

		# Sort atoms by element and generate list of elements in molecule
		atomsbyelement = {}
		elements = set()
		for aid in self.atoms():
			e = self.elementof(aid)
			elements.add(e)
			if e in aidsbyelement:
				aidsbyelement[e].append(aid)
			else:
				aidsbyelement[e] = [aid]
		self.struct['deriv']['aidsbyelement'] = aidsbyelement
		self.struct['deriv']['elements'] = elements

		# Sort bonds by atom
		bidsbyaid = [[]]*len(self.getatoms())
		for bid in self.getbonds():
			aids = self._getbond(bid)['atoms']
			bidsbyaid[aids[0]].append(bid)
			bidsbyaid[aids[1]].append(bid)
		self.struct['deriv']['bidsbyaid'] = bidsbyaid


class MolDB:
	'''
	Represents a database of molecules
	'''

	def __init__(self, path):
		'''
		Define database

			Parameters:
				path (Path or string): path of database folder

			Returns:
				nothing
		'''
		self.path = Path(path)
		self.cids = self.path.glob('*.json')

	def getcids(filters):
		'''
		Return list of cids present in database that pass given filters

			Parameters:
				filters ([Filter]): list of filters to apply

			Returns:
				list of cids ([int])
		'''
		cids = self.cids
		for molfilter in filters:
			cids = filter(molfilter.check, cids)
		return list(cids)


	def filtermolecs:



	def getmolecules(filters):
		'''
		Get list of initialized Molecule objects from database that pass given filters

			Parameters:
				filters ([Filter]): list of filters to apply

			Returns:
				list of molecules ([Molecule])
		'''
		cids = self.getcids(filters)
		molecs = []
		for cid in cids:
			try:


	def add(*molecules):
		for cid in cids:
			m = Molecule(cid)
			m.load()
			m.

	def getmolecule(cid):
		pass

#####################################
# Filters
# Gets passed specs list and compound, returns whether compound passes the check
# Check pass by default if check is not being used

class MoleculeFilter:
	'''
	Abstract superclass of all filters (Do not instantiate)
	'''
	def check(mol):
		'''
		Check if molecule passes filter
		'''
		return False

class AtomCountFilter(MoleculeFilter):
	'''
	Filter for molecules based on number of atoms/atoms of a certain element
	'''

	def __init__(locount, hicount, element=None){
		self.locount = locount
		self.hicount = hicount
		self.element = getatomicnum(element)
	}

	def check(mol):
		count = len(mol.getAtoms(element=self.element))
		return count >= locount and count <= hicount

class BondedElementFilter(MoleculeFilter):
	'''
	Filter for molecules based on bonded elements
	'''

	def __init__(element, allowed):
		self.element = element
		self.allowed = [getatomicnum(e) for e in allowed]

	def check(mol):
		for atom in mol.getAtoms(element=self.element):
			for bond in atom.bonds:
				if bond.partnerof(atom).element not in self.allowed:
					return False
		return True


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

	#Generate names
	namesrecorded = 0
	names = b''
	while (namesrecorded < 3 and namesrecorded < len(compound.synonyms)):
		names += (str.encode(compound.synonyms[namesrecorded]) + b',')
		namesrecorded += 1
	names = names[:-1]

	xyz_data = str(len(atoms)).encode() + b'\n' + names + b'\n'
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

	#Generate names
	namesrecorded = 0
	names = b''
	while (namesrecorded < 3 and namesrecorded < len(compound.synonyms)):
		names += (str.encode(compound.synonyms[namesrecorded]) + b',')
		namesrecorded += 1
	names = names[:-1]

	mol_data = sdf_data[:i] + names + sdf_data[i:sdf_data.index(b'M  END')+len('M  END\n')]

	return mol_data


#Retrieve files in batches from a list of cids. Lowers risk of PUG REST Timeout
def retrieve_by_cids(cids):
	compounds = []
	for cid in cids:
		print(cid)
		compounds.append(pcp.Compound.from_cid(cid))
	return compounds

#Retrieve files
def retrieve_files(specs, dbpath, filterlist, fileformat='xyz', maxresults=1000):
	#Load database log
	try:
		with open(dbpath/'dblog.json', 'r') as f:
			dblog = json.load(f)
	except:
		print('Local database log not found. Initializing directory as new database.')
		dblog = {}

	searchterm = specs[0][1]    # Assumed to be first line of specs

	#Retrieve search results from PubChem
	search_result_cids = pcp.get_cids(searchterm , namespace='smiles', searchtype='substructure', \
		MaxRecords=100000, record_type='3d')

	if search_result_cids is None:
		print('Error retrieving search results from PubChem')
		sys.exit()

	numdownloaded = 0
	random.shuffle(search_result_cids)
	shuffledcids = iter(search_result_cids)

	#Filter results
	while numdownloaded < maxresults:
		try:
			cid = shuffledcids.__next__()
		except:
			print('Did not achieve maximum requested results')
			break

		#Check that file not already downloaded
		#Generate filepath of structure
		filepath = dbpath/(str(cid)+'.'+fileformat)
		#Check if already in local db (don't redownload)
		if filepath.exists():
			continue

		compound = pcp.Compound.from_cid(cid)

		# Do checks
		passcheck = True
		for filterfunc in filterlist:
			if not filterfunc(specs, compound):
				passcheck = False
				break

		if passcheck:
			#Convert to correct format and save based on format
			print('Saving ' + str(compound.cid))
			#Retrieve SDF from PubChem
			try:
				file_data = pcp.get(compound.cid, output='SDF', record_type='3d')
			except Exception as e:
				print(e)
				continue

			if file_data is None:
				print(f"Could not download SDF for {compound.cid}")
				continue

			#Apply any requested format conversions
			if (fileformat.lower() == 'mol'):
				file_data = convert_sdf_to_mol(file_data, compound)
			elif (fileformat.lower() == 'xyz'):
				file_data = convert_mol_to_xyz(convert_sdf_to_mol(file_data, compound), compound)

			#Write file
			with open(filepath, 'wb') as f:
				f.write(file_data)

			#Add 3 most common names with CID to database log
			numnames = 0
			while (numnames < 3 and numnames < len(compound.synonyms)):
				dblog[compound.synonyms[numnames].lower()] = compound.cid
				numnames += 1

			numdownloaded +=1


	#Save log file
	with open(dbpath/'dblog.json', 'w') as f:
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
	except:
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
	commandaction.add_argument('-l', '--lookup',
				metavar = 'Name',
				type=str,
				help='Compound name to search for in local database or online')
	commandaction.add_argument('-d', '--download',
				metavar = 'SpecsFile',
				type=str,
				help='Download molecules using specs from file')
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
	if args.download is not None:
		filterlist = [check_charge,
			check_max_size,
			check_atom_counts,
			check_bonded_elements,
			check_num_bonds_between_elements]	#Should put this specs file

		specs = get_specs(Path(args.download))
		print(specs)
		retrieve_files(specs, Path(args.path), filterlist,
			maxresults=args.maxresults, fileformat=args.format.lower())

	#Lookup molecules locally
	elif args.lookup is not None:
		lookup_cid(args)

#Run main if called as command
if __name__ == "__main__":
	main()
