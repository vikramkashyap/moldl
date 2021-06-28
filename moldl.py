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
from copy import copy

##########################################
# Utility Functions
##########################################

PERIODIC_TABLE = \
	['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
	'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V',
	'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
	'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',
	'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
	'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
	'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
	'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am',
	'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh',
	'Hs', 'Mt', 'Ds ', 'Rg ', 'Cn ', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

def getatomicnum(element):
	'''
	Get atomic number of element

		Parameters:
			element (int or str): atomic number or symbol of element

		Returns:
			atomic number of element (int)
	'''

	# If argument is already int, assume it is atomic num and pass back
	if type(element) is int:
		if element > 0 and element < len(PERIODIC_TABLE):
			return element
		else:
			raise ValueError("Invalid atomic number")
	# If string, find in list of element symbols
	elif type(element) is str:
		try:
			return PERIODIC_TABLE.index(element.capitalize())+1
		except:
			raise ValueError("Invalid element symbol")
	# Perhaps in future allow identification by other types?
	else:
		raise TypeError( \
			"Elements can only be identified by atomic number or symbol")


def getelementsymbol(element):
	if type(element) is int:
		try:
			return PERIODIC_TABLE[element-1]
		except:
			raise ValueError("Invalid atomic number")
	elif type(element) is str:
		if element.capitalize() in PERIODIC_TABLE:
			return element.capitalize()
		else:
			raise ValueError("Invalid element symbol")
	else:
		raise TypeError( \
			"Elements can only be identified by atomic number or symbol")

##########################################
# Core Classes
##########################################


class PropertyBlob:
	'''
	Contains properties and ID
	Superclass of Atom, Bond, and Molecule
	'''
	def __init__(self, **kwargs):
		self._setattribute('_properties', {})
		for key in kwargs.keys():
			self.__setattr__(key, kwargs[key])


	@property
	def properties(self):
		return copy(self._properties)


	def __getattr__(self, key):
		'''
		Get attributes of Atom
		either directly from properties dictionary or aliased

			Paramters:
				key (str): name of attribute

			Returns:
				attribute corresponding to key
		'''
		return self.getProperty(key)

	def __setattr__(self, key, val):
		'''
		Set attributes of Atom
		either directly through properties dictionary or aliased

			Paramters:
				key (str): name of attribute
				val: value of attribute

			Returns:
				None
		'''
		if key == '_properties': super().__setattr__(key, val)
		else: self.setProperty(key, val)


	def __eq__(self, other):
		if self.id == None or other.id == None: return self is other
		return self.id == other.id


	def getProperty(self, key):
		if key in self._properties: return self._properties[key]
		else: return None


	def setProperty(self, key, val):
		self._properties[key] = val


	def _setattribute(self, key, val):
		object.__setattr__(self, key, val)

class Atom(PropertyBlob):
	'''
	Represents an individual atom in a molecule

		Attributes:
			element (int): atomic number of atom
			bonds ([Bond]): list of Bonds that Atom is part of
	'''

	def __init__(self, element=None, **kwargs):
		'''
		Define atom by the molecule it is contained in and its atom ID within molecule

			Parameters:

			Returns:
				Object representing atom in molecule (Atom)
		'''
		super().__init__(**kwargs)
		if element is not None:
			self.element = getatomicnum(element)


class Bond(PropertyBlob):
	def __init__(self, atom1, atom2, **kwargs):
		'''
		Represents a bond between Atoms

			Parameters:
				atom1 (Atom): one Atom in bond
				atom2 (Atom): other Atom in bond

			Returns:
				Object representing bond between atom1 and atom2
		'''
		super().__init__(**kwargs)
		self._setattribute('_atoms', (atom1, atom2))

	@property
	def atoms(self):
		return copy(self._atoms)

	def partnerof(self, atom):
		if self._atoms[0] is atom:
			return self._atoms[1]
		elif self._atoms[1] is atom:
			return self._atoms[0]
		else:
			return ValueError('Given Atom is not part of Bond')




class Molecule(PropertyBlob):
	'''
	Represents a certain molecule and stores its properties as attributes

		Attributes:
	'''

	def __init__(self, id, bonds=[], **kwargs):
		'''
			Define Molecule by name
			Parameters:
				properties ({}): properly JMF formatted dictionary of properties
				download (bool): whether to download from PubChem
						if no properties passed (only set False if you know
						what you're doing)

			Returns:
				Molecule object
		'''
		super().__init__(**kwargs)
		self.id = str(id)
		self._setattribute('_bonds', [])
		self._setattribute('_atoms', [])
		self._setattribute('_bidcounter', 0)
		self._setattribute('_aidcounter', 0)
		for bond in bonds: self.addBond(bond)


	@staticmethod
	def from_jmf(f):
		jmf = json.load(f)
		atoms = {}
		for a in jmf['atoms']:
			atom = Atom()
			for property in a.keys():
				atom.setProperty(property, a[property])
			atoms[atom.id] = atom
		bonds = []
		for b in jmf['bonds']:
			bond = Bond(atoms[b['aids'][0]], atoms[b['aids'][1]])
			for property in b.keys():
				bond.setProperty(property, b[property])
			bonds.append(bond)
		mol = Molecule(jmf['id'], bonds)
		jmf.pop('atoms')
		jmf.pop('bonds')
		for property in jmf.keys():
			mol.setProperty(property, jmf[property])
		return mol

	@staticmethod
	def from_cid(cid):
		'''
		Create Molecule by retrieving information from PubChem by CID
		'''
		try:
			pcpcompound = pcp.Compound.from_cid(cid, record_type='3d')
		except Exception:
			return None
		# Fail if 3d compound not available
		# TODO: put option to include 2d only structures
		if pcpcompound.atoms[0].z is None:
			return None
		atoms = {}
		for pcpatom in pcpcompound.atoms:
			a = Atom(pcpatom.number)
			a.coords = [pcpatom.x, pcpatom.y, pcpatom.z]
			a.formalcharge = pcpatom.charge
			a.coordinate_type = pcpatom.coordinate_type
			atoms[pcpatom.aid] = a
		bonds = []
		for pcpbond in pcpcompound.bonds:
			b = Bond(atoms[pcpbond.aid1], atoms[pcpbond.aid2])
			b.order = pcpbond.order
			bonds.append(b)
		# TODO: Submit fix to pubchempy bug that crashes to_dict() on 3d structures
		properties = pcp.Compound.from_cid(cid, record_type='2d').to_dict()
		properties.pop('atoms')
		properties.pop('bonds')
		mol = Molecule(cid, bonds, **properties)
		return mol

	@property
	def atoms(self):
		return copy(self._atoms)


	@property
	def bonds(self):
		return copy(self._bonds)


	def addBond(self, bond):
		if bond in self.bonds: return
		bond.id = self._bidcounter
		self._setattribute('_bidcounter', self._bidcounter+1)
		self._bonds.append(bond)
		self._addAtom(bond.atoms[0])
		self._addAtom(bond.atoms[1])
		#self._bondpropertymaps['Atom'][bond.atoms[0]].add(bond)
		#self._bondpropertymaps['Atom'][bond.atoms[1]].add(bond)


	def _addAtom(self, atom):
		if atom not in self.atoms:
			atom.id = self._aidcounter
			self._setattribute('_aidcounter', self._aidcounter+1)
			self._atoms.append(atom)
			#self._atompropertymaps['element'][atom.element].add(atom)


	def getAtoms(self, **kwargs):
		'''
		Get a list of atoms in molecule filtered to match a set of properties
		'''
		atoms = []
		for atom in self.atoms:
			selected = True
			for key in kwargs.keys():
				val = kwargs[key]
				if atom.getProperty(key) != val: selected = False
			if selected: atoms.append(atom)
		return atoms


	def getBonds(self, atom=None, **kwargs):
		bonds = []
		for bond in self.bonds:
			selected = True
			if atom is None or atom in bond.atoms:
				for key in kwargs.keys():
					val = kwargs[key]
					if bond.key != val: selected = False
			else: selected = False
			if selected: bonds.append(bond)
		return bonds


	def to_json(self):
		moldict = self.properties
		moldict['atoms'] = []
		for atom in self.atoms:
			moldict['atoms'].append(atom.properties)
		moldict['bonds'] = []
		for bond in self.bonds:
			bonddict = bond.properties
			bonddict['aids'] = (bond.atoms[0].id, bond.atoms[1].id)
			moldict['bonds'].append(bonddict)
		return json.dumps(moldict, indent=4)


	def to_xyz(self):
		return Molecule.jmf_to_xyz(self.to_json())

	@classmethod
	def jmf_to_xyz(cls, jmf):
		jmf = json.loads(jmf)
		xyz = str(len(jmf['atoms']))
		xyz += "\n\n"
		for atom in jmf['atoms']:
			xyz += f"{getelementsymbol(atom['element'])}\t" + \
			f"{atom['coords'][0]}\t{atom['coords'][1]}\t{atom['coords'][2]}\n"
		return xyz

class MolDB:
	'''
	Represents a database of molecules
	'''

	def __init__(self, path, molecules=[]):
		'''
		Define database

			Parameters:
				path (Path or string): path of database folder

			Returns:
				nothing
		'''
		self.path = Path(path)
		if not self.path.exists():
			self.path.mkdir()
		self.molecules = [str(file.name).rstrip('.jmf') \
		 	for file in self.path.glob('*.jmf')]
		self._loadedmols = {}
		for mol in molecules:
			self.addMolecule(mol)

	def getMolecule(self, id):
		if id in self._loadedmols:
			return self._loadedmols[id]
		else:
			mol = Molecule.from_jmf((self.path/(f"{id}.jmf")).open('r'))
			self._loadedmols[mol.id] = mol
			return mol

	def addMolecule(self, mol):
		if mol.name in self.molecules:
			raise ValueError('Molecule already exists in database')
		else:
			self._loadedmols.append(mol)
			mol.save(self.path/(mol.id+'.jmf'))
			self.molecules.append(mol.id)


	def filter(self, filters):
		'''
		Return list of cids present in database that pass given filters

			Parameters:
				filters ([Filter]): list of filters to apply

			Returns:
				list of cids ([int])
		'''
		filteredids = []
		for id in self.molecules:
			mol = self.getMolecule(id)
			passed = True
			for molfilter in filters:
				if molfilter.check(mol) is False:
					passed = False
					break
			if passed:
				filteredids.append(id)
		return filteredids


	def save(self, mol):
		data = mol.to_json()
		path = self.path/f"{mol.cid}.jmf"
		with path.open('w') as f:
			f.write(data)

	def extract(self, ids, path, format='xyz'):
		'''
		Save moleculuar structures in various formats
		'''

		path = Path(path)
		if not path.exists():
			path.mkdir()
		for id in ids:
			data = self.getMolecule(id).to_xyz()
			with (path/f"{id}.xyz").open('w') as f:
				f.write(data)


	def searchPubChem(self, searchterm='', filters=[], numresults=10,\
		randomized=True, save=True):
		'''
		Get list of initialized Molecule cids from database that pass given filters

			Parameters:
				filters ([Filter]): list of filters to apply

			Returns:
				list of cids ([int])
		'''
		#Retrieve search results from PubChem
		searchcids = \
			pcp.get_cids(searchterm , namespace='smiles', \
			searchtype='substructure', MaxRecords=10000, record_type='3d')
		if randomized: random.shuffle(searchcids)

		results = []
		for cid in searchcids:
			if len(results) == numresults: break
			if str(cid) in self.molecules: continue
			print(f"fetching molecule {cid}")
			mol = Molecule.from_cid(cid)
			if mol is None: continue
			print("saving")
			self.save(mol)
			passed = True
			for molfilter in filters:
				if not molfilter.check(mol):
					passed = False
					break
			if passed:
				print('passed')
				results.append(mol.id)

		return results


#####################################
# Filters
# Gets passed specs list and compound, returns whether compound passes the check
# Check pass by default if check is not being used

class MoleculeFilter:
	'''
	Abstract superclass of all filters (Do not instantiate)
	'''
	def check(self, mol):
		'''
		Check if molecule passes filter
		'''
		return False


class AtomCountFilter(MoleculeFilter):
	'''
	Filter for molecules based on number of atoms/atoms of a certain element
	'''

	def __init__(self, locount, hicount, element=None):
		self.locount = locount
		self.hicount = hicount
		if element is not None:
			self.element = getatomicnum(element)
		else: self.element = None


	def check(self, mol):
		if self.element == None: count = len(mol.getAtoms())
		else: count = len(mol.getAtoms(element=self.element))
		return count >= self.locount and count <= self.hicount


class BondedElementFilter(MoleculeFilter):
	'''
	Filter for molecules based on bonded elements
	'''

	def __init__(self, element, allowed):
		self.element = getatomicnum(element)
		self.allowed = [getatomicnum(e) for e in allowed]

	def check(self, mol):
		for atom in mol.getAtoms(element=self.element):
			for bond in mol.getBonds(atom):
				if bond.partnerof(atom).element not in self.allowed:
					return False
		return True


class BondCountFilter(MoleculeFilter):
	'''
	Filter molecules based on number of certain type of bond
	'''

	def __init__(self, element1, element2):
		pass
