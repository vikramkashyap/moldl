class T:
	def __init__(self, partner):
		self.p = partner

	def __getattr__(self, key):
		if key == 'goodstuff':
			print('heres summa da good stuff')
t = T(1)
print(t.p)
t.goodstuff
t.a