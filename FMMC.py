#!/usr/bin/python

# things to do
# non-random? output is always the same?
# stereochemistry check
# replace checksame by RMS comparison in RDKit - the getTorsion is not required
# get rid of writing classes and include in main loop
# fix realtime

# Standard Python Libraries #
import sys, os, random, math, datetime, time
import numpy as np

# Non-Standard Python Libraries #
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import ForceField

# The time elapsed between two specified Y/M/D 24H/M/S format #
def RealTime(time1, time2):
	FMT = "%H:%M:%S"
	#t_delta = datetime.strptime(str(time2), FMT) - datetime.strptime(str(time1), FMT)
	#return [t_delta.days, t_delta.hours, t_delta.minutes,t_delta.seconds]
	return [0,0,0,0]

#Pythagoras and Simple Trig #
def calcdist(atoma,atomb,coords):
	a = np.array([coords[atoma][0], coords[atoma][1],coords[atoma][2]])
	b = np.array([coords[atomb][0], coords[atomb][1],coords[atomb][2]])
	return np.linalg.norm(numpy.array(a)-numpy.array(b))

#Angle between two vectors (degrees)#
def calcangle(atoma,atomb,atomc,coords):
	u = np.array([coords[atoma][0]-coords[atomb][0], coords[atoma][1]-coords[atomb][1],coords[atoma][2]-coords[atomb][2]])
	v = np.array([coords[atomc][0]-coords[atomb][0], coords[atomc][1]-coords[atomb][1],coords[atomc][2]-coords[atomb][2]])
	c = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v)
	return 180.0/math.pi * np.arccos(clip(c, -1, 1))

#Dihedral angle between two bonds (degrees) #
def calcdihedral(atoma,atomb,atomc,atomd,coords):
	ab = np.array([coords[atomb][0]-coords[atoma][0], coords[atomb][1]-coords[atoma][1],coords[atomb][2]-coords[atoma][2]])
	bc = np.array([coords[atomc][0]-coords[atomb][0], coords[atomc][1]-coords[atomb][1],coords[atomc][2]-coords[atomb][2]])
	cd = np.array([coords[atomd][0]-coords[atomc][0], coords[atomd][1]-coords[atomc][1],coords[atomd][2]-coords[atomc][2]])
	vec1 = np.cross(ab,bc)/np.linalg.norm(ab)/np.linalg.norm(bc)
	vec2 = np.cross(bc,cd)/np.linalg.norm(bc)/np.linalg.norm(cd)
	torsion=180.0/math.pi*np.arccos(np.dot(vec1,vec2)/np.linalg.norm(vec1)/np.linalg.norm(vec2))
	sign=180.0/math.pi*np.arccos(np.dot(vec1,ab)/np.linalg.norm(vec1)/np.linalg.norm(ab))
	if sign<90.0: torsion=torsion*-1.0
	return torsion

# Returns a matrix of dihedral angles (with sign) given connecitivty and coordinates in numerical order. Only between heavy atoms and NH,OH and SH protons
# to be removed when checkSame is revised
def getTorsion(MolSpec):
	torval=[]
	for atoma in range(0,len(MolSpec.CARTESIANS)):
		for partner1 in MolSpec.CONNECTIVITY[atoma]:
			atomb = int(partner1.split("__")[0])-1
			for partner2 in MolSpec.CONNECTIVITY[atomb]:
				atomc = int(partner2.split("__")[0])-1
				if atomc!=atoma:
					for partner3 in MolSpec.CONNECTIVITY[atomc]:
						atomd = int(partner3.split("__")[0])-1
						if atomd>atoma and atomd!=atomb:
							if MolSpec.ATOMTYPES[atoma]=="H" and MolSpec.ATOMTYPES[atomb]=="C": ignore=1
							elif MolSpec.ATOMTYPES[atomd]=="H" and MolSpec.ATOMTYPES[atomc]=="C": ignore=1
							else:
								endA=""
								endD=""
								for endAatom in MolSpec.CONNECTIVITY[atomb]:
									if (int(endAatom.split("__")[0])-1)!=atomc:
										endA = endA+MolSpec.ATOMTYPES[int(endAatom.split("__")[0])-1]
								for endDatom in MolSpec.CONNECTIVITY[atomc]:
									if (int(endDatom.split("__")[0])-1)!=atomb:endD = endD+MolSpec.ATOMTYPES[int(endDatom.split("__")[0])-1]
								if endA!="HHH" and endD!="HHH":
									torsion=calcdihedral(atoma,atomb,atomc,atomd,MolSpec.CARTESIANS)
									#print (MolSpec.ATOMTYPES[atoma],atoma, MolSpec.ATOMTYPES[atomb],atomb, MolSpec.ATOMTYPES[atomc],atomc, MolSpec.ATOMTYPES[atomd],atomd, torsion)
									torval.append(torsion)
	return torval

# Filter post optimization - checks whether two conformers are identical on the basis of non-bonded distances and energy. Needs to consider equivalent coordinate descriptions
# also do enantiomers here
def checkSame(ConfSpec, CSearch, SearchParams, savedconf):
	
	if not hasattr(ConfSpec, "CARTESIANS"): return 1
	
	tordiff=0.0; besttordiff=180.0; sameval=0
	
	if len(SearchParams.EQUI)==0:
		#print "No equivalent coordinate descriptions"
		torval1=getTorsion(ConfSpec)
		#print ConfSpec.NAME
		#print ConfSpec.CARTESIANS
		#print torval1
		#print CSearch.NAME[savedconf]
		#print CSearch.CARTESIANS[savedconf]
		#print CSearch.TORVAL[savedconf]
		for x in range(0,len(torval1)):
			difftor=math.sqrt((torval1[x]-CSearch.TORVAL[savedconf][x])*(torval1[x]-CSearch.TORVAL[savedconf][x]))
			if difftor>180.0:
				difftor=360.0-difftor
			tordiff=tordiff + difftor*difftor
		
		# print torval1[x], CSearch.TORVAL[savedconf][x], difftor
		if len(torval1)!=0:
			besttordiff=math.sqrt(tordiff/len(torval1))
	#print "----------"
	#print besttordiff
	#print "----------"

	else:
		tempcart=[]
		for i in range(0,len(ConfSpec.CARTESIANS)): tempcart.append(ConfSpec.CARTESIANS[i])
		allposscoords=[tempcart]
		#print len(allposscoords)
		#print SearchParams.EQUI
		equilist=[]
		for i in range(0,len(SearchParams.EQUI)): equilist.append(SearchParams.EQUI[i])
		for equivcoords in equilist:
			equivstring=equivcoords.split(" ")
			if len(equivstring)==3:
				equilist.append(equivstring[0]+" "+equivstring[1])
				equilist.append(equivstring[0]+" "+equivstring[2])
				equilist.append(equivstring[1]+" "+equivstring[2])
		#print equilist
		for equivcoords in equilist:
			#print equivcoords
			equivstring=equivcoords.split(" ")
			nequiv=len(equivstring)
			if nequiv == 2:
				equivloop=[]
				for item in equivstring:
					equivloop.append(int(item)-1)
				for item in equivstring:
					equivloop.append(int(item)-1)
				
				for l in range(0,len(allposscoords)):
					for i in range(1, nequiv):
						orig=[]
						swap=[]
						for j in range(0,nequiv):
							orig.append(equivloop[j])
							swap.append(equivloop[i+j])
						
						swappedcoords=[0]*len(ConfSpec.CARTESIANS)
						for j in range(0,len(ConfSpec.CARTESIANS)):
							swappedcoords[j]=allposscoords[l][j]
							for k in range(0,nequiv):
								if j == orig[k]:
									#print "Interchanging coordinates",j, swap[k]
									swappedcoords[j]=allposscoords[l][swap[k]]
						#print swappedcoords[j]
						
						#PossSpec = ConfSpec
						#PossSpec.CARTESIANS = swappedcoords
						#torval1=getTorsion(PossSpec)
						#print torval1
						allposscoords.append(swappedcoords)
		#print "appending"

		#print swappedcoords
		#print len(allposscoords)
		originaltorval=getTorsion(ConfSpec)
		originalsum = 0.0
		alteredsum = 0.0
		for x in range(0,len(originaltorval)):
			originalsum = originalsum + math.pow(originaltorval[x],2.0)
		for poss in allposscoords:
			#print "poss"
			#for i in range(0,len(poss)):
			#	print ConfSpec.ATOMTYPES[i], poss[i][0], poss[i][1], poss[i][2]
			PossSpec = ConfSpec
			PossSpec.CARTESIANS = poss
			torval1=getTorsion(PossSpec)
			tordiff=0.0
			alteredsum = 0.0
			for y in range(0,len(torval1)):
				alteredsum = alteredsum + math.pow(torval1[y],2.0)

			if math.pow((originalsum-alteredsum),2.0) < 0.1:
				#print originalsum, alteredsum
				#print "comparing torsions"
				for z in range(0,len(torval1)):
					#print savedconf, "-", x, torval1[x], CSearch.TORVAL[savedconf][x]
					difftor=math.sqrt((torval1[z]-CSearch.TORVAL[savedconf][z])*(torval1[z]-CSearch.TORVAL[savedconf][z]))
					if difftor>180.0:
						difftor=360.0-difftor
					tordiff=tordiff + difftor*difftor
					#print torval1[z], CSearch.TORVAL[savedconf][z], difftor
				if len(torval1)!=0:
					tordiff=math.sqrt(tordiff/len(torval1))
					#print tordiff
					if tordiff<besttordiff:
						besttordiff=tordiff


		#this is a horrible hack which returns the cartesians back to before equivalent coordinate systems were considered....
		ConfSpec.CARTESIANS = tempcart

	if besttordiff<SearchParams.COMP: sameval=sameval+1
	return sameval


# Find how many separate molecules there are
def howmanyMol(bondmatrix,startatom):
	molecule1=[]
	count = 1
	stop=0
	currentatom=[]
	nextlot=[]
	currentatom.append([-1])
	currentatom.append([startatom])
	while count<100 and stop==0:
		nextlot=[]
		for onecurrentatom in currentatom[count]:
			for partners in bondmatrix[onecurrentatom]:
					inf = partners.split("__")
					for n in range(0,len(inf)/2):
							noback=0
							for onepreviousatom in currentatom[count-1]:
								if (int(inf[2*n])-1)==onepreviousatom: # Can't go back - make sure atom isn't in previous currentatom
									noback=noback+1
								if noback==0:
									nextlot.append(int(inf[2*n])-1)
		count=count+1
		if len(nextlot)==0: stop=stop+1
		currentatom.append(nextlot)

	for i in range(0,len(bondmatrix)):
		for j in range(0,len(currentatom)):
			for atom in currentatom[j]:
				if atom==i:
					same = 0
					for alreadyfound in molecule1:
						if atom == alreadyfound: same = same + 1
					if same ==0: molecule1.append(atom)
	return molecule1

# Looks for rotatable single bonds. Requires connectivity information. Uninteresting torsions (e.g. methyl groups) are excluded
class Assign_Variables:
	def __init__(self, MolSpec, PARAMS,log):
		self.MOLATOMS = [range(0,MolSpec.NATOMS)]
		self.NMOLS = 1
		# Find rotatable bonds
		self.ETOZ = []
		#if len(PARAMS.ETOZ) > 0:
		#	self.ETOZ.append([int(Params.ETOZ[0].split()[0]), int(Params.ETOZ[0].split()[1])])
		self.TORSION = []; ring = []
		for i in range(0, MolSpec.NATOMS):
			for partner in MolSpec.CONNECTIVITY[i]:
				nextatom = int(partner.split("__")[0])-1
				bondorder = (partner.split("__")[1])
				
				if nextatom>i: # Avoid duplication
					
					fixed=0
					
					for fix in PARAMS.FIXT:
						fixA=int(fix.split(" ")[0])
						fixB=int(fix.split(" ")[1])
						if fixA == i+1 and fixB == nextatom+1: fixed=fixed+1
						if fixB == i+1 and fixA == nextatom+1: fixed=fixed+1
					
					# Must be single bonds not specified as fixed
					if bondorder == "SINGLE" and fixed == 0:
						terminal = 0
						CX3 = 0
						nextCX3 = 0
						for member in [i, nextatom]:
							nextatomstring = ""
							nextnextatomstring = ""
							nextnumh = 0
							templist1 = []
							
							#print member, MolSpec.CONNECTIVITY[member]
							# Each atom must have other atoms attached
							if len(MolSpec.CONNECTIVITY[member])>1:
								terminal = terminal + 1
								for nextpartner in MolSpec.CONNECTIVITY[member]:
									for z in range(0,len(nextpartner.split("__"))/2):
										nextnextatom = int(nextpartner.split("__")[2*z])-1
										if nextnextatom != i and nextnextatom != nextatom:
											nextatomstring=nextatomstring+MolSpec.ATOMTYPES[nextnextatom]
											templist1.append(nextnextatom)
											for nextnextpartner in MolSpec.CONNECTIVITY[nextnextatom]:
												for v in range(0,len(nextnextpartner.split("__"))/2):
													nextnextnextatom = int(nextnextpartner.split("__")[2*v])-1
													if nextnextnextatom!=member:
														nextnextatomstring=nextnextatomstring+MolSpec.ATOMTYPES[nextnextnextatom]
									
										#Checks if either of two connected atoms are CX3 groups
										if nextatomstring.find("HHH") > -1 or nextatomstring.find("FFF") > -1 or nextatomstring.find("ClClCl") >- 1 or nextatomstring.find("III")>-1: CX3 = CX3 + 1
										#Checks if either of two connected atoms are bonded to three identical CX3 or NX3 or PX3 or SiX3 groups
										if nextatomstring.find("CCC")>-1 or nextatomstring.find("NNN")>-1 or nextatomstring.find("PPP")>-1 or nextatomstring.find("SiSiSi")>-1:
											if nextnextatomstring.find("HHHHHHHHH")>-1 or nextnextatomstring.find("FFFFFFFFF")>-1 or nextnextatomstring.find("ClClClClClClClClCl")>-1 or nextnextatomstring.find("IIIIIIIII")>-1: nextCX3=nextCX3+1
			
						# Given the above criteria are satisfied, we also need to make sure the torsion isn't part of a ring, as they must be dealt with differently
						if CX3 == 0 and nextCX3 == 0 and terminal == 2:
							count = 1
							mem = 0
							currentatom=[[i],[nextatom]]
							nextlot=[]
							while count<1000 and mem==0:
								nextlot=[]
								for onecurrentatom in currentatom[count]:
									for partners in MolSpec.CONNECTIVITY[onecurrentatom]:
										inf = partners.split("__")
										for n in range(0,len(inf)/2):
											noback=0
											for onepreviousatom in currentatom[count-1]:
												if (int(inf[2*n])-1) == onepreviousatom: noback = noback + 1
												if (int(inf[2*n])-1) == nextatom: noback = noback + 1
											
											if noback == 0:
												if (int(inf[2*n])-1) == i: mem = count+1
												nextlot.append(int(inf[2*n])-1)
								count=count+1
								currentatom.append(nextlot)
									
							if mem == 0: self.TORSION.append([i,nextatom])


		self.MCNV = len(self.TORSION)
		
		#A random number of changes are made to generate each new structure. Chang, Guida and Still found between 2 and ntors-1 to work well
		self.MCNVmin = 2
		self.MCNVmax = self.MCNV - 1
		
		#However if there is only one or two rotatable bonds, adjust the upper limit to ntors
		if self.MCNVmin > self.MCNV: self.MCNVmin = self.MCNV
		if self.MCNVmax < self.MCNVmin: self.MCNVmax = self.MCNVmin
		
		#If there is nothing to vary then exit
		if self.MCNV == 0 and self.MCRI == 0 and self.NMOLS == 1: print ("\nFATAL ERROR: Found zero rotatable torsions and only one molecule in %s  \n"%MolSpec.NAME); sys.exit()


#Rotates specified single bond through a specified angle, returning the modified coordinates
def AtomRot(MolSpec, torsion, geometry):
	
	atomA = geometry[torsion[0]-1]
	atomB = geometry[torsion[1]-1]
	theta = float(torsion[2])/180.0*math.pi
	cyclic = 0
	newcoord=[]
	
	if cyclic == 0:
		vectorAB = [float(atomB[0])-float(atomA[0]),float(atomB[1])-float(atomA[1]),float(atomB[2])-float(atomA[2])]
		distAB = math.sqrt(vectorAB[0]*vectorAB[0]+vectorAB[1]*vectorAB[1]+vectorAB[2]*vectorAB[2])
		unitAB = [vectorAB[0]/distAB,vectorAB[1]/distAB,vectorAB[2]/distAB]
		
		count = 1
		stop=0
		currentatom=[]
		nextlot=[]
		currentatom.append([torsion[0]-1])
		currentatom.append([torsion[1]-1])
		
		while count<100 and stop==0:
			nextlot=[]
			for onecurrentatom in currentatom[count]:
				for partners in MolSpec.CONNECTIVITY[onecurrentatom]:
					inf = partners.split("__")
					for n in range(0,len(inf)/2):
						noback=0
						for onepreviousatom in currentatom[count-1]:
							if (int(inf[2*n])-1)==onepreviousatom: noback=noback+1
						for onepreviousatom in currentatom[count]:
							if (int(inf[2*n])-1)==onepreviousatom: noback=noback+1
						if noback==0: nextlot.append(int(inf[2*n])-1)
			count=count+1
			#print count
			if len(nextlot) == 0:stop=stop+1
			currentatom.append(nextlot)
		
		#print "current atom", currentatom
		for i in range(0,len(geometry)):
			#print i,geometry[i]
			newcoord.append(geometry[i])
		for i in range(2,len(currentatom)-1):
			for atom in currentatom[i]:
				#print "Rotating", (atom+1), "about", torsion[0]
				dotproduct = unitAB[0]*(float(geometry[atom][0]) - float(atomA[0])) + unitAB[1]*(float(geometry[atom][1]) - float(atomA[1])) + unitAB[2]*(float(geometry[atom][2]) - float(atomA[2]))
				centre = [float(atomA[0]) + dotproduct*unitAB[0], float(atomA[1]) + dotproduct*unitAB[1], float(atomA[2]) + dotproduct*unitAB[2]]
				v = [float(geometry[atom][0]) - centre[0], float(geometry[atom][1]) - centre[1], float(geometry[atom][2]) - centre[2]]
				d = math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
				px = v[0]*math.cos(theta) + v[1]*math.sin(theta)*unitAB[2] - v[2]*math.sin(theta)*unitAB[1]
				py = v[1]*math.cos(theta) + v[2]*math.sin(theta)*unitAB[0] - v[0]*math.sin(theta)*unitAB[2]
				pz = v[2]*math.cos(theta) + v[0]*math.sin(theta)*unitAB[1] - v[1]*math.sin(theta)*unitAB[0]
				newv = [px + centre[0], py + centre[1], pz + centre[2]]
				newdist = math.sqrt(px*px + py*py + pz*pz)
				newcoord[int(atom)]=newv

	if len(newcoord) !=0 : return newcoord
	else:
		print "didn't do anything!!!"
		for i in range(0,len(geometry)): newcoord.append([0.0,0.0,0.0])
		return newcoord


class OrderConfs:
	def __init__(self, CSEARCH, PARAMS, start, log):
		#Order the low energy conformers by energy
		self.CARTESIANS = []
		self.NAME = []
		self.TIMESFOUND = []
		self.USED = []
		self.TORVAL =[]
		self.MATCHEDALREADY = []
		for j in range(0, CSEARCH.NSAVED): self.MATCHEDALREADY.append(0)
		
		self.ENERGY = sorted(CSEARCH.ENERGY)
		for j in range(0, CSEARCH.NSAVED):
			for i in range(0, CSEARCH.NSAVED):
				if CSEARCH.ENERGY[i] == self.ENERGY[j] and self.MATCHEDALREADY[i] == 0:
					match = i
			
			self.MATCHEDALREADY[match] = 1
			self.CARTESIANS.append(CSEARCH.CARTESIANS[match])
			self.NAME.append(CSEARCH.NAME[match])
			self.TIMESFOUND.append(CSEARCH.TIMESFOUND[match])
			self.USED.append(CSEARCH.USED[match])
			self.TORVAL.append(CSEARCH.TORVAL[match])
		
		CSEARCH.CARTESIANS = self.CARTESIANS
		CSEARCH.NAME = self.NAME
		CSEARCH.ENERGY = self.ENERGY
		CSEARCH.TIMESFOUND = self.TIMESFOUND
		CSEARCH.USED = self.USED
		CSEARCH.TORVAL = self.TORVAL

class AddConformer:
	def __init__(self, CSEARCH, CONFSPEC):
		CSEARCH.NAME.append(CONFSPEC.NAME)
		CSEARCH.ENERGY.append(CONFSPEC.ENERGY)
		CSEARCH.CARTESIANS.append(CONFSPEC.CARTESIANS)
		CSEARCH.CONNECTIVITY.append(CONFSPEC.CONNECTIVITY)
		CSEARCH.TIMESFOUND.append(1)
		CSEARCH.USED.append(0)
		CSEARCH.TORVAL.append(getTorsion(CONFSPEC))
		CSEARCH.NSAVED = CSEARCH.NSAVED + 1


class RemoveConformer:
	def __init__(self, CSEARCH, todel):
		j=0
		for i in range(0,len(CSEARCH.NAME)): print CSEARCH.NAME[i], CSEARCH.ENERGY[i]
		print todel, len(todel),
		cutoff = (len(CSEARCH.NAME)-len(todel))
		print cutoff
		newtodel=[]
		for i in range(len(todel)-1, -1, -1): newtodel.append(todel[i])
		print CSEARCH.NAME[cutoff:]
		print CSEARCH.NAME[:cutoff]
		
		for i in range(0,len(todel)):
			#print i, todel[i], CSEARCH.TIMESFOUND[todel[i]]
			CSEARCH.NREJECT = CSEARCH.NREJECT + CSEARCH.TIMESFOUND[todel[i]]
		
		del CSEARCH.NAME[cutoff:]
		del CSEARCH.ENERGY[cutoff:]
		del CSEARCH.CARTESIANS[cutoff:]
		del CSEARCH.CONNECTIVITY[cutoff:]
		del CSEARCH.USED[cutoff:]
		del CSEARCH.TIMESFOUND[cutoff:]
		del CSEARCH.TORVAL[cutoff:]
		print "AFTER REMOVAL"
		for i in range(0,len(CSEARCH.NAME)): print CSEARCH.NAME[i], CSEARCH.ENERGY[i]
		CSEARCH.NSAVED = len(CSEARCH.NAME)


#       Formatted output to command line and log file         #
class FMLog:
	# Designated initializer
	def __init__(self,filein,suffix,append):
		# Create the log file at the input path
		self.log = open(filein+"_"+append+"."+suffix, 'w' )
	
	# Write a message to the log
	def Write(self, message):
		# Print the message
		print message
		
		# Write to log
		self.log.write(message + "\n")
	
	# Write a message only to the log and not to the terminal
	def Writeonlyfile(self, message):
		# Write to log
		self.log.write(message)
	
	# Write a fatal error, finalize and terminate the program
	def Fatal(self, message):
		# Print the message
		print message+"\n"
		
		# Write to log
		self.log.write(message + "\n")
		
		# Finalize the log
		self.Finalize()
		
		# End the program
		sys.exit(1)
	
	# Finalize the log file
	def Finalize(self):
		self.log.close()


class Writeintro:
	# Formatted text printed to terminal and log file at the beginning of a new search
	def __init__(self, MolSpec, Params, Variables, time, log):
		
		strucname = MolSpec.NAME.split("_step_0")[0]
		torstring=""
		for torsion in Variables.TORSION: torstring = torstring+"{"+str(MolSpec.ATOMTYPES[torsion[0]])+str(torsion[0]+1)+"-"+str(MolSpec.ATOMTYPES[torsion[1]])+str((torsion[1]+1))+"} "
		
		fixtstring=""
		for fixed in Params.FIXEDATOMS: fixtstring = fixtstring+"{"+str(MolSpec.ATOMTYPES[fixed[0]-1])+str(fixed[0])+"-"+str(MolSpec.ATOMTYPES[fixed[1]-1])+str((fixed[1]))+"} "
		
		molarray=[]
		for mol in Variables.MOLATOMS:
			molstring ="{ "
			for atom in mol: molstring =molstring+str(MolSpec.ATOMTYPES[int(atom)])+str(int(atom)+1)+" "
			molstring =molstring+"} "
			molarray.append(molstring)
		
		equistring=""
		for equi in Params.EQUI:
			equistring =equistring+"{ "
			for atom in equi.split(" "): equistring = equistring+str(MolSpec.ATOMTYPES[int(atom)-1])+str(int(atom))+" "
			equistring =equistring+"} "
		
		log.Write(dashedline+"\n   |    "+("FULL_MONTE search on "+strucname).ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("COMP: "+str(Params.COMP)+" degrees").ljust(leftcol)+("|").rjust(rightcol))
		if len(Params.FIXT) > 0:
			log.Write("   | o  "+("FIXT: Manually constrained "+str(len(Params.FIXT))+" torsional variables").ljust(leftcol)+("|").rjust(rightcol)); log.Write("   |    "+fixtstring.ljust(leftcol)+("|").rjust(rightcol))

		if Variables.NMOLS > 1:
			log.Write("   | o  "+("Detected "+str(Variables.NMOLS)+" separate molecules - this adds additional search coordinates").ljust(leftcol)+("|").rjust(rightcol));
			for i in range(0,len(molarray)):
				chunks, chunk_size = len(molarray[i]), leftcol
				for j in range(0, chunks, chunk_size):
					log.Write("   |    "+(molarray[i][j:j+chunk_size]).ljust(leftcol)+("|").rjust(rightcol))
		
		log.Write("   | o  "+("LEVL: "+str(Params.LEVL)+" force field").ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("DEMX: "+str(Params.DEMX)+" kcal/mol").ljust(leftcol)+("|").rjust(rightcol))
		if len(Params.EQUI) > 0:
			log.Write("   | o  "+("EQUI: The following sets of atoms are equivalent ").ljust(leftcol)+("|").rjust(rightcol))
			log.Write("   |    "+equistring.ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("EWIN: "+str(Params.EWIN)+" kcal/mol").ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("MCNV: "+str(Variables.MCNV)).ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   |    "+torstring.ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("STEP: "+str(Params.MAXSTEP)).ljust(leftcol)+("|").rjust(rightcol))
		log.Write(dashedline+"\n")


class WriteSummary:
	# Formatted text printed to terminal and log file at the end of each search step
	def __init__(self, CSearch, SearchParams, start, log):
		now = time.strftime("%H:%M:%S", time.localtime())
		runningtime = RealTime(start, now)
		
		if CSearch.COMPLETE == 0: log.Write("\no  "+("STEP "+str(CSearch.STEP)+" COMPLETE: "+str(CSearch.NSAVED)+" unique conformations. Global minimum energy = "+str(round(CSearch.GLOBMIN,5))).ljust(leftcol)+("").rjust(rightcol))
		if CSearch.COMPLETE == 1: log.Write("\no  "+("FULL MONTE SEARCH COMPLETE: "+str(CSearch.NSAVED)+" unique conformations. Global minimum energy = "+str(round(CSearch.GLOBMIN,5))).ljust(leftcol)+("").rjust(rightcol))
		
		log.Write(("\n     Conformer Name             Absolute Energy       Erel (kcal/mol)        Times found    Times used   ").ljust(leftcol)+"\n"+dashedline)
		
		for i in range(0, CSearch.NSAVED):
			absenergy = str(round(float(CSearch.ENERGY[i]),5))
			if len(absenergy.split(".")[1])!=5: absenergy = absenergy+"0"
			relenergy = str(round(float((CSearch.ENERGY[i]-CSearch.GLOBMIN)),2))
			if len(relenergy.split(".")[1])!=2: relenergy = relenergy+"0"
			log.Write("     "+CSearch.NAME[i].ljust(30)+(absenergy).ljust(20)+(relenergy).rjust(10)+ (str(CSearch.TIMESFOUND[i])).rjust(15)+(str(CSearch.USED[i])).rjust(15)+("").rjust(2))
		log.Write(dashedline+"\n     o  "+("Execute time: "+str(runningtime[0])+"d "+str(runningtime[1])+"h "+str(runningtime[2])+"m "+str(runningtime[3])+"s").ljust(leftcol)+("").rjust(rightcol))
		if CSearch.COMPLETE == 0: log.Write("     o  "+("SE = "+str(round(float(CSearch.AERATE),1))+"    DMIN = "+str(CSearch.DMIN)+"    NOPT = "+str(CSearch.STEP)+"    NFAIL = "+str(CSearch.NFAILED)).ljust(leftcol)+("").rjust(rightcol)+"\n"+ dashedline)
		if CSearch.COMPLETE == 1: log.Write("     o  "+("SE = "+str(round(float(CSearch.AERATE),1))+"    DMIN = "+str(CSearch.DMIN)+"    NOPT = "+str((CSearch.STEP-1))+"    NFAIL = "+str(CSearch.NFAILED)).ljust(leftcol)+("").rjust(rightcol)+"\n"+ dashedline)


class makePDBformat:
	#Write a PDF file for viewing that contains the low energy conformations in ascending order of energy
	def __init__(self, filein, MolSpec, CSearch,append):
		pdbfile = open(filein+"_"+append+".pdb", 'w' )
		if CSearch.NSAVED > 0:
			for i in range(0, CSearch.NSAVED):
				Erel = (CSearch.ENERGY[i]-CSearch.GLOBMIN)
				pdbfile.write("COMPND    "+CSearch.NAME[i]+"   E = "+str(Erel)+"\nAUTHOR    FULL MONTE SEARCH - paton.chem.ox.ac.uk")
				for j in range(0, MolSpec.NATOMS):
					x = "%.3f" % CSearch.CARTESIANS[i][j][0]
					y = "%.3f" % CSearch.CARTESIANS[i][j][1]
					z = "%.3f" % CSearch.CARTESIANS[i][j][2]
					pdbfile.write("\nHETATM"+str(j+1).rjust(5)+MolSpec.ATOMTYPES[j].rjust(3)+"   LIG     1".rjust(10)+x.rjust(12)+y.rjust(8)+z.rjust(8))
				pdbfile.write("\nEND    \n")
			pdbfile.close()


class makeSDFformat:
	#Write a PDF file for viewing that contains the low energy conformations in ascending order of energy
	def __init__(self, filein, MolSpec, CSearch,append):
		sdffile = open(filein+"_"+append+".sdf", 'w' )
		if CSearch.NSAVED > 0:
			for i in range(0, CSearch.NSAVED):
				Erel = (CSearch.ENERGY[i]-CSearch.GLOBMIN)
				sdffile.write(CSearch.NAME[i]+"\n")
				sdffile.write("     E="+str(Erel)+"\n\n")
				sdffile.write(str(MolSpec.NATOMS).rjust(3)+str(MolSpec.NBONDS).rjust(3)+"  0  0  0  0  0  0  0  0  0999 V2000")
				for j in range(0, MolSpec.NATOMS):
					x = "%.4f" % CSearch.CARTESIANS[i][j][0]
					y = "%.4f" % CSearch.CARTESIANS[i][j][1]
					z = "%.4f" % CSearch.CARTESIANS[i][j][2]
					sdffile.write("\n"+x.rjust(10)+y.rjust(10)+z.rjust(10)+MolSpec.ATOMTYPES[j].rjust(2)+"   0  0  0  0  0  0  0  0  0  0  0  0")
				for atomi in range(0,MolSpec.GetNumAtoms()):
						for atomj in range(atomi,MolSpec.GetNumAtoms()):
							if MolSpec.GetBondBetweenAtoms(atomi,atomj):
								sdffile.write("\n"+str(atomi+1).rjust(3)+str(atomj+1).rjust(3)+str(int(MolSpec.GetBondBetweenAtoms(atomi,atomj).GetBondTypeAsDouble())).rjust(2)+" 0")
				sdffile.write("\nM  END\n$$$$    \n")
			sdffile.close()


# Formatting
dashedline = "   ------------------------------------------------------------------------------------------------------------------"
emptyline = "   |                                                                                                                     |"
normaltermination = "\n   -----------------       N   O   R   M   A   L      T   E   R   M   I   N   A   T   I   O   N      ----------------\n"
leftcol=97
rightcol=12

def asciiArt(now):
	print "     ___       ___                                    ___          ___          ___                   ___ "
	print "    /  /\\     /__/\\                                  /__/\\        /  /\\        /__/\\         ___     /  /\\"
	print "   /  /:/_    \\  \\:\\                                |  |::\\      /  /::\\       \\  \\:\\       /  /\\   /  /:/_"
	print "  /  /:/ /\\    \\  \\:\\   ___     ___  ___     ___    |  |:|:\\    /  /:/\\:\\       \\  \\:\\     /  /:/  /  /:/ /\\"
	print " /  /:/ /:/___  \\  \\:\\ /__/\\   /  /\\/__/\\   /  /\\ __|__|:|\\:\\  /  /:/  \\:\\  _____\\__\\:\\   /  /:/  /  /:/ /:/_"
	print "/__/:/ /://__/\\  \\__\\:\\\\  \\:\\ /  /:/\\  \\:\\ /  /://__/::::| \\:\\/__/:/ \\__\\:\\/__/::::::::\\ /  /::\\ /__/:/ /:/ /\\"
	print "\\  \\:\\/:/ \\  \\:\\ /  /:/ \\  \\:\\  /:/  \\  \\:\\  /:/ \\  \\:\\~~\\__\\/\\  \\:\\ /  /:/\\  \\:\\~~\\~~\\//__/:/\\:\\\\  \\:\\/:/ /:/"
	print " \\  \\::/   \\  \\:\\  /:/   \\  \\:\\/:/    \\  \\:\\/:/   \\  \\:\\       \\  \\:\\  /:/  \\  \\:\\  ~~~ \\__\\/  \\:\\\\  \\::/ /:/"
	print "  \\  \\:\\    \\  \\:\\/:/     \\  \\::/      \\  \\::/     \\  \\:\\       \\  \\:\\/:/    \\  \\:\\          \\  \\:\\\\  \\:\\/:/"
	print "   \\  \\:\\    \\  \\::/       \\__\\/        \\__\\/       \\  \\:\\       \\  \\::/      \\  \\:\\          \\__\\/ \\  \\::/"
	print "    \\__\\/     \\__\\/                                  \\__\\/        \\__\\/        \\__\\/                 \\__\\/ "
	print "  ", now


class PARAMS: pass
PARAMS.MAXSTEP = 100
PARAMS.LEVL = "UFF"
PARAMS.COMP = 10
PARAMS.FIXT = []
PARAMS.FIXEDATOMS = []
PARAMS.EQUI = []
PARAMS.EWIN=20.0
PARAMS.DEMX=41.84

# Define conformational search statistics
class CSEARCH: pass
# Names from Chang, Guida and Still's definitions
CSEARCH.NREJECT = 0
CSEARCH.NFAILED = 0
CSEARCH.AERATE = 0
CSEARCH.DMIN = 1
CSEARCH.STEP = 20
CSEARCH.NAME = []
CSEARCH.CARTESIANS = []
CSEARCH.CONNECTIVITY = []
CSEARCH.USED = [0]
CSEARCH.TIMESFOUND = [1]
CSEARCH.NSAVED = 1
CSEARCH.COMPLETE = 0

if __name__ == "__main__":
	# An input file must be specified - format must be MOL #
	if len(sys.argv)>1: 
		filein = sys.argv[1].split(".")[0]
		if len(sys.argv[1].split(".mol"))>1: filetype = sys.argv[1].split(".")[1]
		else: print "MOL file name required"; sys.exit()
		for i in range(1,len(sys.argv)):
			if sys.argv[i] == "-step": PARAMS.MAXSTEP = int(sys.argv[i+1])
			elif sys.argv[i] == "-levl": PARAMS.LEVL = (sys.argv[i+1])
	else: print "\nWrong number of arguments used. Correct format: FullMonte molecule.mol \n"; sys.exit()


	# Initialize the logfile for all text output #
	if os.path.exists(filein+"_fm.dat"):
		var = raw_input("\no  Log file already exists! OK to overwrite this file ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": print "   Overwriting ..."
		else: print "\nExiting\n";  sys.exit(1)
	log = FMLog(filein,"dat", "fm")

	# Open the structure file #
	log.Write("\no  Extracting structure from "+filein+"."+filetype+" ...")
	if filetype == "mol": MOLSPEC = Chem.MolFromMolFile('./'+filein+'.mol', removeHs=False)
	MOLSPEC.NAME = filein

	# Model Chemistry to be used
	for level in ["UFF", "MMFF"]:
		if PARAMS.LEVL.upper() == level: JOBTYPE = level
	log.Write("\no  Using "+JOBTYPE+" force field ... ")

	# Perform an optimization of the starting geometry #
	if JOBTYPE == "MMFF" or JOBTYPE == "UFF":
		AllChem.EmbedMolecule(MOLSPEC)
		if JOBTYPE == "UFF":
			if AllChem.UFFHasAllMoleculeParams(MOLSPEC):
				ff = AllChem.UFFGetMoleculeForceField(MOLSPEC)
				AllChem.UFFOptimizeMolecule(MOLSPEC)
		if JOBTYPE == "MMFF":
			if AllChem.MMFFHasAllMoleculeParams(MOLSPEC):
				ff = AllChem.MMFFGetMoleculeForceField(MOLSPEC,AllChem.MMFFGetMoleculeProperties(MOLSPEC))
				AllChem.MMFFOptimizeMolecule(MOLSPEC)
		MOLSPEC.ENERGY = ff.CalcEnergy()
	else: log.Fatal("\nFATAL ERROR"%file)

	MOLSPEC.ATOMTYPES = []
	MOLSPEC.CONNECTIVITY = []
	MOLSPEC.CARTESIANS = []
	MOLSPEC.CHARGE = Chem.GetFormalCharge(MOLSPEC)
	MOLSPEC.NATOMS = MOLSPEC.GetNumAtoms()
	MOLSPEC.NBONDS = MOLSPEC.GetNumBonds()
	for atom in MOLSPEC.GetAtoms(): MOLSPEC.ATOMTYPES.append(atom.GetSymbol())
	for atom in range(0,MOLSPEC.NATOMS):
		pos = MOLSPEC.GetConformer().GetAtomPosition(atom)
		print pos.x, pos.y, pos.z
		MOLSPEC.CARTESIANS.append([pos.x, pos.y, pos.z])

	for atomi in range(0,MOLSPEC.GetNumAtoms()):
		MOLSPEC.CONNECTIVITY.append([])
		for atomj in range(0,MOLSPEC.GetNumAtoms()):
			if MOLSPEC.GetBondBetweenAtoms(atomi,atomj): MOLSPEC.CONNECTIVITY[atomi].append(str(atomj+1)+"__"+str(MOLSPEC.GetBondBetweenAtoms(atomi,atomj).GetBondType()))
	#print Chem.FindMolChiralCenters(MOLSPEC,force=True,includeUnassigned=True)


	# Assign variable torsions, number of separate molecules and ##
	# If number of steps is not assigned use 3^rotatable torsions #
	FMVAR = Assign_Variables(MOLSPEC, PARAMS, log)

	if not PARAMS.MAXSTEP: PARAMS.MAXSTEP = int(math.pow(3,FMVAR.MCNV))
	start = time.strftime("%H:%M:%S", time.localtime())
	asciiArt(start); Writeintro(MOLSPEC, PARAMS, FMVAR, start, log)

	CONFSPEC = MOLSPEC
	CSEARCH.NAME.append(MOLSPEC.NAME)
	CSEARCH.CARTESIANS.append(MOLSPEC.CARTESIANS)
	CSEARCH.TORVAL = [getTorsion(MOLSPEC)]
	CSEARCH.CONNECTIVITY.append(MOLSPEC.CONNECTIVITY)
	CSEARCH.ENERGY = [MOLSPEC.ENERGY]
	CSEARCH.GLOBMIN = MOLSPEC.ENERGY
	CSEARCH.LASTFOUND = 0
	CSEARCH.NSAVED = 1
	CSEARCH.STEP = 0


	# Stop once number of steps exceeded or no new conformers found	
	while CSEARCH.STEP <= PARAMS.MAXSTEP:
		
	# Setting the geometry that will be altered to generate new conformers
		for i in range(0, CSEARCH.NSAVED):
			if CSEARCH.ENERGY[i] - CSEARCH.GLOBMIN == 0.0: startgeom = i

		# Generate new geometries
		CONFSPEC.NAME = filein+"_step_"+str(CSEARCH.STEP)
		
		for j in range(0, CSEARCH.NSAVED):
			if (CSEARCH.ENERGY[j] - CSEARCH.GLOBMIN) < PARAMS.EWIN:
				if CSEARCH.USED[j] < CSEARCH.USED[startgeom]: startgeom = j
				if CSEARCH.USED[j] == CSEARCH.USED[startgeom] and CSEARCH.ENERGY[j] < CSEARCH.ENERGY[startgeom]: startgeom = j
		CSEARCH.USED[startgeom] = CSEARCH.USED[startgeom] + 1

		CONFSPEC.CARTESIANS = []
		# The coordinates of the lowest energy, least used structure will be altered
		log.Write("o  STEP "+str(CSEARCH.STEP)+": Generating structure from "+CSEARCH.NAME[startgeom]+" ...")
		for i in range (0,len(CSEARCH.CARTESIANS[startgeom])):
			CONFSPEC.CARTESIANS.append([])
			for cart in (CSEARCH.CARTESIANS[startgeom][i]): CONFSPEC.CARTESIANS[i].append(cart)
		CONFSPEC.CONNECTIVITY = CSEARCH.CONNECTIVITY[startgeom]
		CONFSPEC.ATOMTYPES = MOLSPEC.ATOMTYPES
		if FMVAR.MCNVmin < FMVAR.MCNVmax: nrandom = random.randint(FMVAR.MCNVmin, FMVAR.MCNVmax)
		else: nrandom = FMVAR.MCNVmax

		if FMVAR.MCNV != 0:
			FMVAR.ADJUST = []
			for dihedral in random.sample(FMVAR.TORSION, nrandom): 
				FMVAR.ADJUST.append([int(dihedral[0])+1, int(dihedral[1])+1, random.randint(0,360)])
			
			if len(FMVAR.ETOZ) > 0:
				ezisomerize = random.choice([0,1])
				for dihedral in random.sample(FMVAR.ETOZ,ezisomerize):
					FMVAR.ADJUST.append([int(dihedral[0]), int(dihedral[1]), 180])
		
			# Take input geometry and apply specified torsional changes
			if hasattr(FMVAR, "ADJUST"):
				print FMVAR.ADJUST
				for torsion in FMVAR.ADJUST: CONFSPEC.CARTESIANS = AtomRot(MOLSPEC, torsion, CONFSPEC.CARTESIANS)


		conf = CONFSPEC.GetConformer(id=-1)
		for atomi in range(0,MOLSPEC.GetNumAtoms()):
			conf.SetAtomPosition(atomi,CONFSPEC.CARTESIANS[atomi])
			print conf.GetAtomPosition(atomi).x, conf.GetAtomPosition(atomi).y, conf.GetAtomPosition(atomi).z


		#for nconf in range(0,CONFSPEC.GetNumConformers()):
		#	for atomi in range(0,CONFSPEC.GetNumAtoms()):
		#		print CONFSPEC.GetConformer(id=nconf-1).GetAtomPosition(atomi).x, CONFSPEC.GetConformer(id=nconf-1).GetAtomPosition(atomi).y, CONFSPEC.GetConformer(id=nconf-1).GetAtomPosition(atomi).z

		# Perform an optimization
		if JOBTYPE == "MMFF" or JOBTYPE == "UFF":
			AllChem.EmbedMolecule(CONFSPEC)
			print Chem.MolToMolBlock(CONFSPEC,confId=-1)

			if JOBTYPE == "UFF":
				if AllChem.UFFHasAllMoleculeParams(CONFSPEC):
					ff = AllChem.UFFGetMoleculeForceField(CONFSPEC)
					AllChem.UFFOptimizeMolecule(CONFSPEC)
			if JOBTYPE == "MMFF":
				if AllChem.MMFFHasAllMoleculeParams(CONFSPEC):
					ff = AllChem.MMFFGetMoleculeForceField(CONFSPEC,AllChem.MMFFGetMoleculeProperties(CONFSPEC))
					AllChem.MMFFOptimizeMolecule(CONFSPEC)
			CONFSPEC.ENERGY = ff.CalcEnergy()
		else: log.Fatal("\nFATAL ERROR"%file)


		CONFSPEC.CARTESIANS = []
		for atom in range(0,MOLSPEC.NATOMS):
			pos = CONFSPEC.GetConformer().GetAtomPosition(atom)
			print pos.x, pos.y, pos.z
			CONFSPEC.CARTESIANS.append([pos.x, pos.y, pos.z])

		#Check whether the molecule has high energy
		if ((CONFSPEC.ENERGY-CSEARCH.GLOBMIN)) < PARAMS.DEMX:
			samecheck = 0
			# also check whether a duplicate conformation has been found
			for j in range(0, CSEARCH.NSAVED):
				if checkSame(CONFSPEC, CSEARCH, PARAMS, j) > 0 : #or checkSame(makemirror(CONFSPEC), CSEARCH, PARAMS, j) > 0:
					log.Write("   "+(CONFSPEC.NAME+" is a duplicate of conformer "+CSEARCH.NAME[j]+" ... ").ljust(50))
					CSEARCH.TIMESFOUND[j] = CSEARCH.TIMESFOUND[j] + 1
					CSEARCH.NREJECT = CSEARCH.NREJECT + 1
					samecheck = samecheck + 1
					break
			
			# Unique conformation with low energy! #
			if samecheck == 0:
				if CONFSPEC.ENERGY < CSEARCH.GLOBMIN:
					CSEARCH.GLOBMIN = CONFSPEC.ENERGY
					log.Write("   "+(CONFSPEC.NAME+" is a new Global Minimum!").ljust(80)+("E = "+str(CSEARCH.GLOBMIN)).rjust(rightcol))
				else : log.Write("   "+(CONFSPEC.NAME+" is saved").ljust(80)+("E = "+str(CONFSPEC.ENERGY)).rjust(rightcol))
				AddConformer(CSEARCH, CONFSPEC)
				if (CONFSPEC.ENERGY-CSEARCH.GLOBMIN) < PARAMS.EWIN: CSEARCH.LASTFOUND = CSEARCH.STEP


		# Rejection - discard #
		else: log.Write("   "+(CONFSPEC.NAME+" is rejected due to high energy ... ").ljust(50)); CSEARCH.NREJECT = CSEARCH.NREJECT + 1
		OrderConfs(CSEARCH, PARAMS, start, log)


		# End of step - update the search statistics #
		if (CSEARCH.STEP-CSEARCH.NFAILED) != 0: CSEARCH.AERATE = float(CSEARCH.STEP-CSEARCH.NREJECT-CSEARCH.NFAILED)/float(CSEARCH.STEP-CSEARCH.NFAILED)*100.0
		else: CSEARCH.AERATE = 0.0
		if len(CSEARCH.TIMESFOUND) > 0:
			for dup in CSEARCH.TIMESFOUND:
				if dup < CSEARCH.DMIN: CSEARCH.DMIN = dup
		else: CSEARCH.DMIN = 0

		#Tidy up the results - if the lowest energy has dropped then it may be necessary to remove some previously saved conformers
		if CSEARCH.STEP % 100 == 0 and CSEARCH.STEP > 0:
			todel=[]
			for i in range(0,len(CSEARCH.NAME)):
				if ((CSEARCH.ENERGY[i] - CSEARCH.GLOBMIN)) > PARAMS.DEMX or (i > 199): todel.append(i)
			if len(todel) !=0: RemoveConformer(CSEARCH, todel)
			WriteSummary(CSEARCH, PARAMS, start, log)	
					
		#End of step - update step number
		CSEARCH.STEP = CSEARCH.STEP + 1

	#Summary of completed Full Monte search #######################
	CSEARCH.COMPLETE = 1
	WriteSummary(CSEARCH, PARAMS, start, log)
	makeSDFformat(filein, MOLSPEC, CSEARCH, "fm")
	end = time.strftime("%H:%M:%S", time.localtime())
	asciiArt(end); log.Write(normaltermination); log.Finalize()	

