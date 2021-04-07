# First line is search key, which is used as a SMILES-formatted substructure search
SEARCH_TERM P
# Only take neutral charged molecules
CHARGE 0 0
# Maximum number of atoms allowed in molecule
MAX_SIZE 30
# Maximum number of atoms allowed for P
ATOM_COUNT P 1 1
# Only allow P to be bonded to H,N,O,C, or S
BONDED_ELEMENTS P H N O C S
# Only let P be bonded to up to 2 Os
NUM_BONDS_BETWEEN_ELEMENTS P O 0 2
