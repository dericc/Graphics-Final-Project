f = open('levels/generated.scn', 'w')

f.write('material \n 1 1 1 \n 1 1 1 \n 0 0 0 \n0 0 0 \n 0 0 0 \n 1 1 stonewall.jpg \n\n')

f.write('material \n 1 1 1 \n 1 1 1 \n 0 0 0 \n 0 0 0 \n 0 0 0 \n 1 1 stonewallRotate.jpg \n\n')

f.write('material 1.000000 1.000000 1.000000 .90000 .700000 0.200000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000 1.000000 coin.jpg \n\n')

f.write('material \n 1 1 1 \n 1 1 1 \n 0 0 0 \n 0 0 0 \n 0 0 0 \n 1 1 CCube.jpg \n\n')

f.write('player  \n 3 \n 0 0 0 \n 1 1 1 \n 5 \n 1 \n\n')

f.write('box \n 1 \n -5 -1 -2 \n 20 0 2 \n\n')

for i in range(0, 14):
	for j in range(0, 7):
		f.write('coin\n2\n' + str(i) + ' ' + str(j + 1) + ' 0 \n\n')
