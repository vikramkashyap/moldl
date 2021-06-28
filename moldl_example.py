import moldl
from moldl import *

db = MolDB('/path/to/dir')

filters = [
    AtomCountFilter(0, 30),
    AtomCountFilter(1, 1, 'P'),
    BondedElementFilter('P', ('H', 'N', 'O', 'C'))
]

ids = db.searchPubChem(searchterm='P', numresults=1000)

filtered_ids = db.filter(filters)

db.extract(ids, '/path/to/output/dir', format='xyz')
