# moldl: A tool for downloading molecular structure files from PubChem
==========
## Installation
`moldl` uses the PubChem API through ```pubchempy``` and gets info about the elements using ```mendeleev```. Install them using pip:

```bash
pip install pubchempy
```

```bash
pip install mendeleev
```

----------
## Use

### Downloading files
`moldl` can be instructed to bulk download structure files (in SDF, MOL, or XYZ) that contain a specified element. Further paramenters can also be specified, such as the number of target element atoms to be present in the molecule. All specifications for what molecules to download are passed in an input file (see below). moldl will also save a JSON file mapping the 3 most common names for the compounds in the local downloaded database to the compounds' CIDS.
#### Input files
moldl downloads structures based on the specifications given to it in an input file. Each line in the input file references a certain "filter" and has the form 

```python
[filtername] [arg1] [arg2] [arg3] ...
```

All currently implemented filters are described below. `SEARCH_TERM` is the only required filter and must be the first line.

```python
SEARCH_TERM [SMILES formatted substructure]
CHARGE [lowest allowed charge] [highest allowed charge]
MAX_SIZE [number of atoms]
ATOM_COUNT [element] [lowest allowed number] [highest allowed number]
BONDED_ELEMENTS [element] [allowed bonded elements]
NUM_BONDS_BETWEEN_ELEMENTS [element] [bonded element] [lowest num bonds] [highest num bonds]
```

Here is an example input file:

```python
# Consider any molecule containing phosphorus
SEARCH_TERM P
# Only take neutral charged molecules (0 charge)
CHARGE 0 0
# Maximum number of atoms allowed in molecule is 30
MAX_SIZE 30
# Only allow one phophorus per molecule
ATOM_COUNT P 1 1
# Only allow P to be bonded to H,N,O,C, or S
BONDED_ELEMENTS P H N O C S
# Only let P be bonded to up to 2 Os
NUM_BONDS_BETWEEN_ELEMENTS P O 0 2
```

### Searching names
moldl can also be passed a name and will search for it in the local database to see if it is downloaded. If it cannot be found, the name will be looked up on PubChem and the CID will be retrieved.

## Help

For information on how to use the `moldl` tool, run 

```
moldl -h
```