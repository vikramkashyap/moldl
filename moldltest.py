import moldl as mdl

database = mdl.MolDB('./database')

filters = [mdl.AtomCountFilter(5, 15), mdl.AtomCountFilter(1, 'P')]

database.fetch('P', filters, 100, randomized=True)

newfilters = [mdl.BondedElementFilter('P', ['H', 'N', 'O', 'C'])]

cids = database.getcids(newfilters)

databse.extractfiles(cids, 'XYZ', './xyzfiles')