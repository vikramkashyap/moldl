# moldl Documentation

`moldl` (MOLecule DownLoader) is a python module that helps in collecting and managing datasets of molecular structure files. Version 1.0 of moldl (the current version) is centered around downloading molecular structure information from the [PubChem database](https://pubchem.ncbi.nlm.nih.gov/). Interaction with the PubChem API is handled using [`pubchempy`](https://pypi.org/project/PubChemPy/).

Check out the pages linked on the side to learn how to get started with `moldl`.

*Full API documentation coming soon!*

# Installation

Currently, in order to use `moldl` you must download "moldl.py" from [https://raw.githubusercontent.com/vikramkashyap/moldl/master/moldl.py](https://raw.githubusercontent.com/vikramkashyap/moldl/master/moldl.py). Include this file in your python path or the same folder as your python script.

The module can then be imported in your python script using

	import moldl

# Molecules

Molecules are represented using `moldl.Molecule` objects. Molecules contain a list of atoms, represented by `moldl.Atom` objects, and a list of bonds, represented by `moldl.Bond` objects.

Import the `Molecule` class using:
	
	from moldl import Molecule

## Loading and Saving Molecules

`Molecule` objects can be loaded from the online PubChem database using the PubChem CID (compound ID) number. This example loads the structure for Nitrous Oxide (CID 948):

	mol = Molecule.from_cid(948)

moldl stores molecule data in JSON files that we call "JMF" files ("JSON Molecule Format"). A molecule can be converted to JMF format like so:

	jmfdata = mol.to_json()

This JMF data could then be saved to a file like any text data.

To load from a JMF file, we pass a corresponding file object to `.from_jmf()`:

	mol = Molecule.from_jmf(jmffile)

## Exporting Molecules

Currently, molecules can only be exported as .XYZ, though support for .MOL and .SDF is coming soon. To get the .XYZ representation of a molecule as string, use

	xyzstring = mol.to_xyz()

## Atoms and Bonds

### `Atom` Objects

Atoms in a molecule are represented by `Atom` objects. `Atom`s have a `element` attribute that stores the atom's atomic number. They also have a `bonds` attribute that is a list of `Bond` objects representing the bonds the atom is part of.

`Atom`s in `Molecule`s loaded from PubChem will have additional attributes:

* `coords`: 3D coordinates of atom within molecule (list [x,y,z]). If coordinates are 2D, z coordinate is set to 0.
* `formalcharge`: Formal charge on atom
* `coordinate_type`: "3d" or "2d" depending on type of coordinates

### `Bond` Objects

Bonds in a molecule are represented by `Bond` objects. `Bond`s have an attribute `atoms` that is a 2-tuple of the `Atoms` that make up that bond.

`Bond`s loaded from PubChem will have additional attributes:

* `order`: order of the bond

### Custom Properties

Custom properties of `Atom`s and `Bond`s can be set freely using attributes. For example, to set property `foo` to `1` of some `Atom a`, simply write `a.foo = 1`. The same can be done for any `Bond`. If an attribute is accessed that has not already been set, `None` is returned.


## Inspecting `Molecule`s

To get a list of `Atom` objects in the `Molecule`, we use
	
	atoms = mol.atoms

To get only atoms with certain properties, we use `getAtoms` and pass keyword arguments according to what properties we want to select for. For example, to select only oxygen atoms, we use

	oxygens = mol.getAtoms(element=8)

To get a list of `Bond` objects in the `Molecule`, we use

	bonds = mol.bonds

To get the bonds that include a certain `Atom` `a`, we use:

	bonds_with_a = mol.getBonds(atom=a)

To get only bonds with certain properties, we use `getBonds` and pass keyword arguments according to what properties we want to select for. For example, to select double bonds, we use

	doublebonds = mol.getBonds(order=2)

# Databases

`moldl` handles sets of molecular structures as "databases". `moldl` stores the molecular structure data of each molecule in a JSON file &mdash; called a "JMF" file for "JSON Molecule Format" &mdash; and in reality a database is simply a directory containing JMF files. Within `moldl` however, databases are represented using `MolDB` objects.

## Opening a Database

To open a database, we simply pass the path of the database directory to the `MolDB` constructor. Any directory can be treated as a database, so making a new database simply means passing an empty directory to the `MolDB` constructor.

	from moldl import MolDB
	db = MolDB('/path/to/database/directory')


## Filtering Molecules

The molecules in the database can be filtered using `MoleculeFilter` objects. For example, to get a list of the ids of molecules in the database that contain a single phosphorus atom, have less than 30 atoms total, and are net neutral, we use:

	filters = [
		AtomCountFilter(1, 1, 'P')
		AtomCountFilter(0, 30)
		ChargeFilter([0])
	]
	ids = db.filter(filters)

See the ["Filters"](#filters) section for more information on what `MoleculeFilters` are provided by `moldl` and how custom filters can be implemented. 

## Downloading Structures

To populate the database, we can search PubChem for molecular structures matching a certain search string. Currently, only substructure searches are supported. This means that to collect 100 structures containing phosphorus, we use:

	db.searchPubChem('P', numresults=100)

If we want to only select certain molecules, for example the same molecules we filter for in the [previous section](#filtering-molecules), we use:

	db.searchPubChem('P', numresults=100, filters=filters)

This will continue checking molecules until it has placed 100 molecules that match the filter setting in the database.

By default, any structures that are downloaded that do not pass the filters are still saved in the database. To not save these, pass `save=false`.

PubChem substructure searches will return results in order of "popularity". To avoid only getting commonly known molecules, pass `randomized=True` to randomize the search results.

# Filters

To learn how to use filters with databases, see [the database section](#filtering-molecules).

Filters are used to check if a molecule passes a certain test. They are implemented as subclasses of `MoleculeFilter` and are required to implement a `check` function that takes the molecule to be checked as the argument. The `check` function returns true if the molecule passes the filter and false if it fails. The filters that are provided with `moldl` are as follows:

* `AtomCountFilter(locount, hicount, element=None)`
	* If `element=None`, returns `True` if number of elements in molecule is between `locount` and `hicount` inclusive
	* If `element` is set to an atomic number (ex. `15`) or element symbol (ex. `P`), returns true if number of atoms of that element is between `locount` and `highcount` inclusive

* `ChargeFilter(acceptedrange)`
	* Returns `True` if `charge` attribute of molecule is within `acceptedrange`, which can be either a `range` or `list` (anything that implements `in`), else returns `False`
	* Note that the `charge` attribute is not guaranteed to be set for every `Molecule`, but will be set for all `Molecules` loaded from PubChem

* `BondedElementFilter(element, allowed)`
	* Returns `True` if all atoms of `element` (specified using atomic number or symbol) are only bonded to elements specified in the `allowed` list

# Example

This is an example of how to use `moldl`.

	#Import moldl classes
	from moldl import *

	#Open database (could be an empty directory)
	db = MolDB('/path/to/dir')

	#Populate the database with 100 ranodm molecules containing phosphorus
	db.searchPubChem(searchterm='P', numresults=100)

	#Create a set of filters
	filters = [
		AtomCountFilter(0, 30)		#Allow between 0 and 30 atoms
		AtomCountFilter(1, 1, 'P')	#Allow only 1 phosphorus atom
		BondedElementFilter('P', ('H', 'N', 'O', 'C'))	#Allow phosphorus to only be bonded to H, N, O, or C
	]

	#Apply filters to database to get list of ids (PubChem CIDs in this case)
	filtered_ids = db.filter(filters)

	#Save the filtered molecules to a different folder in XYZ format
	db.extract(filtered_ids, '/path/to/output/dir', format='xyz')
