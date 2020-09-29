# moldl: A tool for downloading molecular structure files from PubChem
==========
## Installation
moldl uses the PubChem API through ```pubchempy``` and gets info about the elements using ```mendeleev```. Install them using pip:

```pip install pubchempy```

```pip install mendeleev```

----------
## Use
For information on how to use the moldl tool, run ```moldl -h```
### Downloading files
moldl can be instructed to bulk download structure files (in SDF, MOL, or XYZ) that contain a specified element. Further paramenters can also be specified, such as the number of target element atoms to be present in the molecule. moldl will also save a JSON file mapping the 3 most common names for the compounds in the local downloaded database to the compounds' CIDS.
### Searching names
moldl can also be passed a name and will search for it in the local database to see if it is downloaded. If it cannot be found, the name will be looked up on PubChem and the CID will be retrieved.
