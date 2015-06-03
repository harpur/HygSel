"""
Description: The GFF_Window.py file works with GFF3_Gamma.txt file and a window file to determine the genes that fall within a specific window in a scaffold. The script also calculates gamma within the scaffold window. 

	Input: The format of the GFF3_Gamma file is as follows...
		Groupchrom	Gene	Start	End		Strand		Genes          Gamma 
		Group1.1	gene	1	100000 		+/-		GB123 		0.420
		
		The format of the window file is as follows...
			Group		Start	End
			Group1.1	1	100000

	Output: The name of the output file is user-defined and has the following format...
		Scaffold	Start	End	Number of Genes	Gamma	Gene List
		1.1		1	1000	23		0.420	GB123, GB456

"""

from  __future__ import division
import sys
import os.path

def snpIterator():
	CPU_NUM = 4
	
	if len(sys.argv) != 2:
        	print >> sys.stderr, "Usage Incorrect!\nCorrect Usage: ./Phased_Script.py Inputfile"
       		sys.exit()
	
	inputfile = sys.argv[1] #window_100kb.txt
	outputfile = inputfile + ".out"
	
	my_file = open(inputfile, "r")
	lines = my_file.readlines() #skips header, the following line also skips the header but has to be used in conjuction with readlines()
	my_file.close()
	
	#snpFiles_in_dir = glob.glob(snpFolderArg) 	

	time.sleep(1)
	procs = []
	    
	while True:
		for proc in procs[:]:
       			if not proc.is_alive(): # Check if process terminated
            			procs.remove(proc)
				
		if not len(lines) and not len(procs):
        		break	
		
		while len(procs) < CPU_NUM and len(lines): #if lines remain in file containing windows and processes are less then CPU_NUM, run the loop
			nuline = lines.pop()
			p = multiprocessing.Process(target=task, args=(nuline, outputfile, ))
			procs.append(p)
			p.start()

	for p in procs:
		p.join()	
	
	# uncomment if you want a txt file with groups sorted
	#sortByGroup(outputFileArg

def task(nuline, outputfile):
	
	chrom = nuline.strip().split()[0]
	chrom_pos = nuline.strip().split()[1]
	pos = nuline.strip().split()[2]
	ref = nuline.strip().split()[3]
	alt = nuline.strip().split()[4]
	
	line_split = nuline.strip().split()
	
	array = [chrom, chrom_pos, pos, ref, alt]
	
	for i in xrange(5, 87) :
		
		position = line_split[i]
		
		if position == "0":
			position = ref
		elif position == "1":
			position = alt
			
		array.append(position)	
			
	outputFileWriter = open(outputfile, 'a')
	outputFileHeader = "\t".join(array) + "\n"
	outputFileWriter.write(outputFileHeader)
	outputFileWriter.close()			
			

if __name__ == "__main__":
	import collections
	from collections import deque
	from collections import defaultdict
	import csv
	import os, os.path
	import time
	import multiprocessing
	import glob	
	snpIterator()
