def elements(self):
  '''
  Get list of elements present in molecule

      Parameters:
        None

      Returns:
        list of atomic numbers of elements present in molecule, without duplicates ([int])
  '''
  return self.struct['deriv']['elements']

def atoms(self, element=None):
	'''
	Get list of atom IDs

		Parameters:
			element (int or string): atomic number or symbol of element to limit list to

		Returns:
			list of atom IDs for all atoms in molecule (default)
			or only atoms of specified element
	'''
	if element is None:
		return list(range(len(self.struct['atoms'])))
	else:
		e = getatomicnum(element)
		if e in self.elements():
			return self.struct['deriv']['aidsbyelement'][getatomicnum(element)]