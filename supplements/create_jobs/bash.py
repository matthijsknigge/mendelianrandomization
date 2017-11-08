#!/usr/bin/env python3.5

import random
import string


def main(args):
	# aquire traits
	trait_list = list()
	
	file = open("traits.txt", "r")
	for line in file:
		# append all traits to list
		line = line.strip()
		trait_list.append(line)
	
	# devide all traits in 1000 jobs
	chunks = lambda l, n: [l[x: x+n] for x in range(0, len(l), n)]
	chunks_trait = (chunks(trait_list, 50))
	
	for chunk in chunks_trait:
		traits = ",".join(chunk)
		name = "".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
		bash = "#!/bin/bash\n#SBATCH --job-name={0}\n#SBATCH --output={0}.out\n#SBATCH --err={0}.err\n#SBATCH --time=02:00:00\n#SBATCH --cpus-per-task=1\n#SBATCH --export=NONE\n#SBATCH --get-user-env=L\n#SBATCH --mem=5gb\n\nmodule load R\n\nRscript mrbase_2010.R -f {1}".format(name, traits)
		name = "{}.sh".format(name)
		output = open(name, "w")
		output.write(bash)
		output.close()


if __name__ == '__main__':
	import sys
	
	
	sys.exit(main(sys.argv))






