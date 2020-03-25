import sys, random

lines = open(sys.argv[1]).readlines()
random.shuffle(lines)
open(sys.argv[2], 'w').write(''.join(lines))
