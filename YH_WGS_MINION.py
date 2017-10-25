#!/usr/bin/env python


##now write something about your progam 

import argparse
import warnings

#### now add all the things we want to have parsed form the command line

parser = argparse.ArgumentParser(
		description="""
		Add here the description of your program
		""")

parser.add_argument('INDIR', help="This is the basesfolder for data analysis\
					it should contain the Albacore output", type=str)

parser.add_argument('--barcodes', help="These are the barcodes that do NOT\
		get combined when using this flag. All others will be combined.\
		Input should be a comma separated list like '01,02,03'", type=str)

args = parser.parse_args()

#get the variables into your script

INDIR= args.INDIR

sep_barcodes = args.barcodes

print(args)

print(barcodes)
exit()
print("All good")