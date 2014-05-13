f = open('levels/generated.scn', 'w')

f.write('material \n 1 1 1 \n 1 1 1 \n 0 0 0 \n0 0 0 \n 0 0 0 \n 1 1 stonewall.jpg \n\n')

f.write('material \n 1 1 1 \n 1 1 1 \n 0 0 0 \n 0 0 0 \n 0 0 0 \n 1 1 stonewallRotate.jpg \n\n')

f.write('material 1.000000 1.000000 1.000000 .90000 .700000 0.200000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000 1.000000 coin.jpg \n\n')

f.write('material \n 1 1 1 \n 1 1 1 \n 0 0 0 \n 0 0 0 \n 0 0 0 \n 1 1 CCube.jpg \n\n')

f.write('coin_material 2\n\n')

f.write('player  \n 3 \n 0 0 0 \n 1 1 1 \n 5 \n 1 \n\n')

f.write('goal 1 \n 28 0 0 \n 29 1 1\n\n')

f.write('box \n 1 \n -5 -1 -2 \n 30 0 2 \n\n')

f.write('skybox black.jpg\n\n')

f.write('soundtrack /../sounds/secret_level.wav\n\n')

for i in range(1, 12):
	for j in range(0, 4):
		f.write('coin\n' + str(i * 2) + ' ' + str(j * 2 + 1) + ' 0 \n\n')
