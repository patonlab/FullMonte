#!/usr/bin/python

###                   ###                     ###          ###      ###
###                   ###                     ###          ###      ###
    #####b.   ####b.  ###### .d##b.  #####b.  ###  ####b.  #####b.  ###
### ### "##b     "##b ###   d##""##b ### "##b ###     "##b ### "##b ###
### ###  ### .d###### ###   ###  ### ###  ### ### .d###### ###  ###
### ### d##P ###  ### Y##b. Y##..##P ###  ### ### ###  ### ### d##P ###
### #####P"  "Y######  "Y### "Y##P"  ###  ### ### "Y###### #####P"  ###
    ###
    ###

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# robert.paton@chem.ox.ac.uk


###############################################################
#                        FMTools.py                           #
#       Libraries and methods for Full Monte Carlo            #
###############################################################
#######  Written by:  Rob Paton ###############################
#######  Last modified:  Mar 20, 2013 #########################
###############################################################

# Python Libraries ############################################
import subprocess, sys, os, commands, math, time, tarfile, random
###############################################################
	
# EXECECTUBALE ################################################
G09_EXEC = 'g09sub'
MOPAC_EXEC = '/u/rsp/rpaton/mopac/MOPAC2012.exe'
###############################################################

# The time elapsed between two specified Y/M/D 24H/M/S format #
def RealTime(time1, time2):
	timeTuple1 = time.strptime(time1, "%Y/%m/%d %H:%M:%S")
	timeTuple2 = time.strptime(time2, "%Y/%m/%d %H:%M:%S")
	time_difference = time.mktime(timeTuple2) - time.mktime(timeTuple1)
	realdays=int(time_difference/(60.0*60*24))
	realhours=int(time_difference/(60.0*60))-realdays*24
	realmins=int(time_difference/60.0)-realdays*24-realhours*60
	realsecs=int(time_difference)-realdays*24*60*60-realhours*60*60-realmins*60
	timediff=[realdays,realhours,realmins,realsecs]
	return timediff
###############################################################

# Tidies up times in [D,H,M,S] format #########################
def CPUTime(CSearch):
	totalcpu = [0,0,0,0]	
	for cpu in CSearch.ALLCPU: 
		for i in range(0,4): totalcpu[i] = totalcpu[i] + cpu[i]
	totalcpu[2] = totalcpu[2]+int(totalcpu[3]/60.0)
	totalcpu[3] = totalcpu[3]-60*int(totalcpu[3]/60.0)
	totalcpu[1] = totalcpu[1] + int(totalcpu[2]/60.0)
	totalcpu[2] = totalcpu[2] - 60*int(totalcpu[2]/60.0)
	totalcpu[0] = totalcpu[0] + int(totalcpu[1]/24.0)
	totalcpu[1] = totalcpu[1] - 24*int(totalcpu[1]/24.0)
	return totalcpu	
###############################################################

# Define Job type #############################################
class JobSpec: 		
	def __init__(self, software):

		self.PROGRAM = software
		self.CONSTRAINED = []
		if software == "Mopac":
			self.EXEC = MOPAC_EXEC
			self.INPUT = ".mop"
			self.ARGS = "&"
		if software == "Gaussian":
			self.EXEC = G09_EXEC
			self.INPUT = " "
			self.ARGS = " "
	
###############################################################
		
# Submits a computational chemistry job #######################
def submitJob(JobSpec,MolSpec,log):			
	if not os.path.exists(MolSpec.NAME+".log"):command = JobSpec.EXEC+" "+MolSpec.NAME+JobSpec.INPUT+" "+JobSpec.ARGS+" > /dev/null"  		
	else: command = ""
	try:
		#print "deactivated submission"
		retcode = subprocess.call(command, shell=True)
		if retcode != 0:
			print >>sys.stderr, log.Write("\nERROR: Submission of "+MolSpec.NAME+" failed")
			return -1
		else: return 1
	except OSError, e:
		print >>sys.stderr, log.Write("\nERROR: Submission of "+MolSpec.NAME+" failed")
		return -1
###############################################################	

# Check that a computational chemistry job has finished #######
def isJobFinished(JobSpec, MolSpec): 
	if JobSpec.PROGRAM == "Mopac":
		if not os.path.exists(MolSpec.NAME+".out"): return 0
		else: 
			outfile = open(MolSpec.NAME+".out","r") 
			jobdone=0; normal=0
			for line in outfile.readlines():
				if JobSpec.PROGRAM == "Mopac":
					if line.find("== MOPAC DONE ==") > -1:
						jobdone = jobdone+1
						normal=normal+1
					if line.find("EXCESS NUMBER OF OPTIMIZATION CYCLES") > -1:
						jobdone = jobdone+1
						normal=normal-1
				outfile.close()

	if JobSpec.PROGRAM == "Gaussian":
		if os.path.exists(MolSpec.NAME+".out"):
                        outfile = open(MolSpec.NAME+".out","r")
                if os.path.exists(MolSpec.NAME+".log"):
                        outfile = open(MolSpec.NAME+".log","r");
		jobdone=0; normal=0
		for line in outfile.readlines():
			if line.find("Normal termination") > -1:
				jobdone = jobdone+1
				normal=normal+1
		outfile.close()

	if jobdone>0 and normal>0: return 1
	if jobdone>0 and normal==0: return 2
	else:
		if JobSpec.PROGRAM == "Mopac": modtime=commands.getoutput("ls -l -t "+MolSpec.NAME+".out")
		if JobSpec.PROGRAM == "Gaussian": 
			if os.path.exists(MolSpec.NAME+".log"):modtime=commands.getoutput("ls -l -t "+MolSpec.NAME+".log")
			if os.path.exists(MolSpec.NAME+".out"):modtime=commands.getoutput("ls -l -t "+MolSpec.NAME+".out")
		#print modtime
		for mod in modtime.split(): 
			if mod.find(":") > -1: timeofday = mod		
		Elapsed = RealTime(time.strftime("%Y/%m/%d" , time.localtime())+" "+modtime.split()[7]+":00", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()))
		ElapsedMins = Elapsed[0]*24*60+Elapsed[1]*60+Elapsed[2]
		if JobSpec.PROGRAM == "Mopac":
			if ElapsedMins < 5: return 0
			else: return -1
###############################################################

# Filter prior to optimization - if there are any very close nonbonded contacts, a non-zero value is returned
def checkDists(MolSpec, SearchParams):
	checkval = 0
	for i in range(0,len(MolSpec.CARTESIANS)): 
		bondedatomlist = []
		for partners in MolSpec.CONNECTIVITY[i]: bondedatomlist.append(int(partners.split("__")[0])-1)
		for j in range(i+1,len(MolSpec.CARTESIANS)):
			bond = 0
			for bondedatom in bondedatomlist:
				if j == bondedatom: bond = bond + bond + 1
			if bond == 0:	
				totdist = abs(calcdist(i, j, MolSpec.CARTESIANS))
				bump = SearchParams.RJCT*(bondiRadius(atomicnumber(MolSpec.ATOMTYPES[i]))+bondiRadius(atomicnumber(MolSpec.ATOMTYPES[j])))
				# If heteroatom - hydrogen bonds are not specified as fixed... 
				if SearchParams.HSWAP != 0:
					if MolSpec.ATOMTYPES[i]=="N" or MolSpec.ATOMTYPES[i]=="O" or MolSpec.ATOMTYPES[i]=="S":
						if MolSpec.ATOMTYPES[j]=="H": bump=0.75*bump
					if MolSpec.ATOMTYPES[j]=="N" or MolSpec.ATOMTYPES[j]=="O" or MolSpec.ATOMTYPES[j]=="S":
						if MolSpec.ATOMTYPES[i]=="H": bump=0.75*bump
				if totdist<bump: 
					checkval = checkval+1
					#print "   PREOPT: Rejecting structure!",MolSpec.ATOMTYPES[i],(i+1),MolSpec.ATOMTYPES[j],(j+1),"distance =", totdist,"Ang"
	return checkval		
###############################################################

# Some useful arrays ##########################################
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
				 "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
				 "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]
atomicmass = [0.0,1.008, 4.003, 6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07, 35.45, 39.95, 39.10, 40.08, 44.96, 47.87, 50.94, 52.00, 54.94, 55.84, 58.93, 58.69, 
			  63.55, 65.39, 69.72, 72.61, 74.92, 78.96, 79.90, 83.80, 85.47, 87.62, 88.91, 91.22, 92.91, 95.94, 99.0, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90, 131.29]
calendar=["","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]

def digitalMonth(month):
	digital = 0
	for i in range(0,len(calendar)):
		if calendar[i] in month.lower(): digital = i
	return digital

def elementID(massno):
	if massno < len(periodictable): return periodictable[massno]	
	else: return "XX"

def atomicnumber(element):
	atomicno = 0
	for i in range(0,len(periodictable)):
		if element == periodictable[i]: atomicno = i
	return atomicno		

def bondiRadius(massno):
	#Bondi van der Waals radii for all atoms from: Bondi, A. J. Phys. Chem. 1964, 68, 441-452, except hydrogen, which is taken from Rowland, R. S.; Taylor, R. J. Phys. Chem. 1996, 100, 7384-7391
	#Radii that are not available in either of these publications have RvdW = 2.00 Angstrom
	
	bondi = [0.0,1.09, 1.40, 1.82,2.00,2.00,1.70,1.55,1.52,1.47,1.54,2.27,1.73,2.00,2.10,1.80,1.80,1.75,1.88,2.75,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.40,1.39,1.87,2.00,1.85,1.90,
			 1.85,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.72,1.58,1.93,2.17,2.00,2.06,1.98,2.16,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.72,1.66,1.55,1.96,2.02,2.00,2.00,2.00,
			 2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.86]
	if massno<len(bondi): radius = bondi[massno]
	else: radius = 2.0 
	return radius
###############################################################

#Geometric calculations #######################################
def calcdist(atoma,atomb,coords):
	x1=coords[atoma][0]
	y1=coords[atoma][1]
	z1=coords[atoma][2]
	x2=coords[atomb][0]
	y2=coords[atomb][1]
	z2=coords[atomb][2]
	ba = [x1-x2, y1-y2, z1-z2]
	dist = math.sqrt(ba[0]*ba[0]+ba[1]*ba[1]+ba[2]*ba[2])
	return dist

def calcangle(atoma,atomb,atomc,coords):
	x1=coords[atoma][0]
	y1=coords[atoma][1]
	z1=coords[atoma][2]
	x2=coords[atomb][0]
	y2=coords[atomb][1]
	z2=coords[atomb][2]
	x3=coords[atomc][0]
	y3=coords[atomc][1]
	z3=coords[atomc][2]
	ba = [x1-x2, y1-y2, z1-z2]
	bc = [x3-x2, y3-y2, z3-z2]
	angle = 180.0/math.pi*math.acos((ba[0]*bc[0]+ba[1]*bc[1]+ba[2]*bc[2])/(math.sqrt(ba[0]*ba[0]+ba[1]*ba[1]+ba[2]*ba[2])*math.sqrt(bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2]))) 
	return angle

def calcdihedral(atoma,atomb,atomc,atomd,coords):
	x1=coords[atoma][0]
	y1=coords[atoma][1]
	z1=coords[atoma][2]
	x2=coords[atomb][0]
	y2=coords[atomb][1]
	z2=coords[atomb][2]
	x3=coords[atomc][0]
	y3=coords[atomc][1]
	z3=coords[atomc][2]
	x4=coords[atomd][0]
	y4=coords[atomd][1]
	z4=coords[atomd][2]
	ax= (y2-y1)*(z2-z3)-(z2-z1)*(y2-y3)
	ay= (z2-z1)*(x2-x3)-(x2-x1)*(z2-z3)
	az= (x2-x1)*(y2-y3)-(y2-y1)*(x2-x3)
	bx= (y3-y2)*(z3-z4)-(z3-z2)*(y3-y4)
	by= (z3-z2)*(x3-x4)-(x3-x2)*(z3-z4)
	bz= (x3-x2)*(y3-y4)-(y3-y2)*(x3-x4)
	nbx= (y2-y3)*(z4-z3)-(z2-z3)*(y4-y3)
	nby= (z2-z3)*(x4-x3)-(x2-x3)*(z4-z3)
	nbz= (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3)
	torsion=180.0/math.pi*math.acos((ax*bx+ay*by+az*bz)/(math.sqrt(ax*ax+ay*ay+az*az)*math.sqrt(bx*bx+by*by+bz*bz)))
	sign=180.0/math.pi*math.acos((nbx*(x2-x1)+nby*(y2-y1)+nbz*(z2-z1))/(math.sqrt(nbx*nbx+nby*nby+nbz*nbz)*math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))))
	if sign<90.0: torsion=torsion*-1.0
	return torsion
###############################################################


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
				print originalsum, alteredsum
				print "comparing torsions"
				for z in range(0,len(torval1)):
					#print savedconf, "-", x, torval1[x], CSearch.TORVAL[savedconf][x]
					difftor=math.sqrt((torval1[z]-CSearch.TORVAL[savedconf][z])*(torval1[z]-CSearch.TORVAL[savedconf][z]))
					if difftor>180.0:
						difftor=360.0-difftor
					tordiff=tordiff + difftor*difftor
					#print torval1[z], CSearch.TORVAL[savedconf][z], difftor
				if len(torval1)!=0:
					tordiff=math.sqrt(tordiff/len(torval1))	
					print tordiff
					if tordiff<besttordiff:
						besttordiff=tordiff
					
		
		#this is a horrible hack which returns the cartesians back to before equivalent coordinate systems were considered....
		ConfSpec.CARTESIANS = tempcart
		
	#print "   ----------"
	#print "   "+str(besttordiff)
	#print "   ----------"
	#print ConfSpec.NAME, CSearch.NAME[savedconf], besttordiff
	if besttordiff<SearchParams.COMP: sameval=sameval+1
	return sameval
	

#Check the stereochemistry has not been changed
def checkchir(ConfSpec, MolSpec, CSearch, SearchParams):
	epimerized = 0; epimatom = 0
	for i in range(0,len(ConfSpec.CARTESIANS)):
		if MolSpec.ATOMTYPES[i] == "C":
			if len(MolSpec.CONNECTIVITY[i]) == 4:
				#print "\nATOM", MolSpec.ATOMTYPES[i], (i),
				abcd = []; types = []
				for partners in MolSpec.CONNECTIVITY[i]: 
					abcd.append(int(partners.split("__")[0])-1)
					types.append(MolSpec.ATOMTYPES[int(partners.split("__")[0])-1])
				numh = 0
				for type in types:
					if type == "H": numh = numh + 1
				
				if numh <= 1:
					#print "   Computing dihedral angle between", abcd 
					if math.fabs(calcdihedral(abcd[0],abcd[1],abcd[2],abcd[3],ConfSpec.CARTESIANS) - calcdihedral(abcd[0],abcd[1],abcd[2],abcd[3],MolSpec.CARTESIANS)) > 10.0:
						print "   POSSIBLY EPIMERIZED!",i, abcd, types
						print calcdihedral(abcd[0],abcd[1],abcd[2],abcd[3],ConfSpec.CARTESIANS), calcdihedral(abcd[0],abcd[1],abcd[2],abcd[3],MolSpec.CARTESIANS)
						epimerized = 1; epimatom = (i+1)
	return [epimerized,epimatom]

#Check the connectivity and compare to the starting structure
def checkconn(ConfSpec, MolSpec, CSearch, SearchParams):
	checkval=0
	a = "X"; b = 0; c = "X"; d = 0
	for i in range(0,len(ConfSpec.CARTESIANS)): 
		bondedatomlist=[]
		for partners in MolSpec.CONNECTIVITY[i]: bondedatomlist.append(int(partners.split("__")[0])-1)
					
		nonbondedlist=[]
		for j in range(i+1,len(ConfSpec.CARTESIANS)):
			bond=0
			nonbondedlist.append(j)
			for bondedatom in bondedatomlist:
							
				#This deals with breaking existing bonds
				if j==bondedatom:
								
					#print (i+1), (j+1)
					nonbondedlist.pop()
					xdist=float(ConfSpec.CARTESIANS[i][0])-float(ConfSpec.CARTESIANS[j][0])
					ydist=float(ConfSpec.CARTESIANS[i][1])-float(ConfSpec.CARTESIANS[j][1])
					zdist=float(ConfSpec.CARTESIANS[i][2])-float(ConfSpec.CARTESIANS[j][2])
					totdist=math.sqrt(xdist*xdist+ydist*ydist+zdist*zdist)
					#print totdist
					origxdist=float(CSearch.CARTESIANS[0][i][0])-float(CSearch.CARTESIANS[0][j][0])
					origydist=float(CSearch.CARTESIANS[0][i][1])-float(CSearch.CARTESIANS[0][j][1])
					origzdist=float(CSearch.CARTESIANS[0][i][2])-float(CSearch.CARTESIANS[0][j][2])
					origdist=math.sqrt(origxdist*origxdist+origydist*origydist+origzdist*origzdist)
					#print origdist
					if (totdist-origdist)>0.5*origdist:
						#print "Looks like ",MolSpec.ATOMTYPES[i],(i+1)," has broken from ",MolSpec.ATOMTYPES[j],(j+1)
						a = MolSpec.ATOMTYPES[i]; b = (i+1); c = MolSpec.ATOMTYPES[j]; d = (j+1)
						checkval=checkval+1                     
										#if HSWAP!=0:
											#		if MolSpec.ATOMTYPES[i]=="N" or MolSpec.ATOMTYPES[i]=="O" or MolSpec.ATOMTYPES[i]=="S":
											#				if MolSpec.ATOMTYPES[j]=="H":
																	#Has an acidic proton swapped positions?
																	#for k in range(0,len(ConfSpec.CARTESIANS)):
																			#if k!=i and k!=j:
																					#hxdist=float(ConfSpec.CARTESIANS[k][0])-float(ConfSpec.CARTESIANS[j][0])
																					#hydist=float(ConfSpec.CARTESIANS[k][1])-float(ConfSpec.CARTESIANS[j][1])
																					#hzdist=float(ConfSpec.CARTESIANS[k][2])-float(ConfSpec.CARTESIANS[j][2])
																					#htotdist=math.sqrt(hxdist*hxdist+hydist*hydist+hzdist*hzdist)
																					#print htotdist
																					#if htotdist<0.5*(bondiradius(atomicnumber(MolSpec.ATOMTYPES[i]))+bondiradius(atomicnumber(MolSpec.ATOMTYPES[j]))):
																							#print "Looks like ",MolSpec.ATOMTYPES[i],(i+1)," has broken from ",MolSpec.ATOMTYPES[j],(j+1)
											#												breakbond.append([i+1,j+1])
											#												checkval=checkval-1     
											#		if MolSpec.ATOMTYPES[j]=="N" or MolSpec.ATOMTYPES[j]=="O" or MolSpec.ATOMTYPES[j]=="S":
											#				if MolSpec.ATOMTYPES[i]=="H":
																	#Has an acidic proton swapped positions?
																	#for k in range(0,len(ConfSpec.CARTESIANS)):
																			#if k!=j and k!=i:
																					#hxdist=float(ConfSpec.CARTESIANS[k][0])-float(ConfSpec.CARTESIANS[i][0])
																					#hydist=float(ConfSpec.CARTESIANS[k][1])-float(ConfSpec.CARTESIANS[i][1])
																					#hzdist=float(ConfSpec.CARTESIANS[k][2])-float(ConfSpec.CARTESIANS[i][2])
																					#htotdist=math.sqrt(hxdist*hxdist+hydist*hydist+hzdist*hzdist)
																					#if htotdist<0.5*(bondiradius(atomicnumber(MolSpec.ATOMTYPES[i]))+bondiradius(atomicnumber(MolSpec.ATOMTYPES[k]))):
																							#print "Looks like ",MolSpec.ATOMTYPES[i],(i+1)," has broken from ",MolSpec.ATOMTYPES[j],(j+1)
																							#breakbond.append([i+1,j+1])
																							#checkval=checkval-1     
		#This deals with forming new bonds
		for j in nonbondedlist:
					
			xdist=float(ConfSpec.CARTESIANS[i][0])-float(ConfSpec.CARTESIANS[j][0])
			ydist=float(ConfSpec.CARTESIANS[i][1])-float(ConfSpec.CARTESIANS[j][1])
			zdist=float(ConfSpec.CARTESIANS[i][2])-float(ConfSpec.CARTESIANS[j][2])
			totdist=math.sqrt(xdist*xdist+ydist*ydist+zdist*zdist)

			if totdist<0.5*(bondiRadius(atomicnumber(MolSpec.ATOMTYPES[i]))+bondiRadius(atomicnumber(MolSpec.ATOMTYPES[j]))):
				print "Looks like ",MolSpec.ATOMTYPES[i],(i+1),":",MolSpec.ATOMTYPES[j],(j+1),"have formed a new bond"
				checkval=checkval+1                     
					
				if SearchParams.NNBO == 0:
					if MolSpec.ATOMTYPES[i]=="N" or MolSpec.ATOMTYPES[i]=="O" or MolSpec.ATOMTYPES[i]=="S":
						if MolSpec.ATOMTYPES[j]=="H": checkval=checkval-1
							
					if MolSpec.ATOMTYPES[j]=="N" or MolSpec.ATOMTYPES[j]=="O" or MolSpec.ATOMTYPES[j]=="S":
						if MolSpec.ATOMTYPES[i]=="H": checkval=checkval-1

	                 
	return [checkval, str(a), str(b),str(c),str(d)]               

	
# Returns a matrix of dihedral angles (with sign) given connecitivty and coordinates in numerical order. Only between heavy atoms and NH,OH and SH protons
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
		if len(nextlot)==0:
			stop=stop+1
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
	#Repeat for all starting atoms. How may different molecules are there?
		
# Looks for rotatable single bonds. Requires connectivity information. Uninteresting torsions (e.g. methyl groups) are excluded			
class Assign_Variables:
		
	def __init__(self, MolSpec, Params, log):
		
		# Find out if there are separate molecules
		self.MOLATOMS = [howmanyMol(MolSpec.CONNECTIVITY,0)]
		for i in range(1,MolSpec.NATOMS):
				k=0
				for j in range(0,len(self.MOLATOMS)):
					if howmanyMol(MolSpec.CONNECTIVITY,i)[0] == self.MOLATOMS[j][0]:
						k=k+1
				if k == 0: self.MOLATOMS.append(howmanyMol(MolSpec.CONNECTIVITY,i))	
		self.NMOLS = len(self.MOLATOMS)
		
		
		# Find rotatable bonds
		self.ETOZ = []
		if len(Params.ETOZ) > 0:
			self.ETOZ.append([int(Params.ETOZ[0].split()[0]), int(Params.ETOZ[0].split()[1])])
		self.TORSION = []
		for i in range(0, MolSpec.NATOMS): 
			for partner in MolSpec.CONNECTIVITY[i]:
				nextatom = int(partner.split("__")[0])-1
				bondorder = float(partner.split("__")[1])
				if nextatom>i: # Avoid duplication	
					fixed=0
					for fix in Params.FIXT:
						fixA=int(fix.split(" ")[0])
						fixB=int(fix.split(" ")[1])
						if fixA == i+1 and fixB == nextatom+1: fixed=fixed+1
						if fixB == i+1 and fixA == nextatom+1: fixed=fixed+1
					
					# Must be single bonds not specified as fixed
					if bondorder == 1.0 and fixed == 0:
						terminal = 0
						CX3 = 0
						nextCX3 = 0
						for member in [i, nextatom]:
							nextatomstring = ""
							nextnextatomstring = ""
							nextnumh = 0
							templist1 = []
							
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
							
							#if mem>5: #Only larger than 5mem rings are interesting conformationally!
									#ring.append([x,nextatom])		
		
		
		self.MCNV = len(self.TORSION) 
		
		#A random number of changes are made to generate each new structure. Chang, Guida and Still found between 2 and ntors-1 to work well
		self.MCNVmin = 2 
		self.MCNVmax = self.MCNV - 1
		
		#However if there is only one or two rotatable bonds, adjust the upper limit to ntors 
		if self.MCNVmin > self.MCNV: self.MCNVmin = self.MCNV
		if self.MCNVmax < self.MCNVmin: self.MCNVmax = self.MCNVmin

		#Ring variables - not coded yet
		self.RING = []
		self.MCRI = 0

		#If there is nothing to vary then exit
		if self.MCNV == 0 and self.MCRI == 0 and self.NMOLS == 1: print ("\nFATAL ERROR: Found zero rotatable torsions and only one molecule in %s  \n"%MolSpec.NAME); sys.exit()			
			    
#Translate a random number of molecules by a given vector
def translateMol(FMVAR, ConfSpec):
	newcoord=[]
	for i in range(0,len(ConfSpec.CARTESIANS)): newcoord.append(ConfSpec.CARTESIANS[i])
	
	
	#Random intermolecular vectors are altered
	nrand = random.randint(1, FMVAR.NMOLS-1)
	movemol = random.sample(xrange(0,FMVAR.NMOLS), nrand)


	for mol in movemol:
		trans = [random.uniform(-1.0,1.0), random.uniform(-1.0,1.0), random.uniform(-1.0,1.0)]
		
		for atom in FMVAR.MOLATOMS[mol]:
			newcoord[atom]=[newcoord[atom][0]+trans[0],newcoord[atom][1]+trans[1],newcoord[atom][2]+trans[2]]
	
	return newcoord		
						
								
#Rotate a molecule about its centre of mass
def rotateMol(FMVAR, ConfSpec):
	newcoord=[]
	for i in range(0,len(ConfSpec.CARTESIANS)): newcoord.append(ConfSpec.CARTESIANS[i])
	
	#Random molecules are spun about their centre of mass
	nrand = random.randint(1, FMVAR.NMOLS-1)
	movemol = random.sample(xrange(0,FMVAR.NMOLS), nrand)
	
	
	for mol in movemol:
		rot = [random.randint(90,180), random.randint(90,180), random.randint(90,180)]
		
		coords1 = []
		types1 = []
		mass1 = 0.0
		
		#print "Spinning molecule",(k+1),"in",savedname[startgeom],"about its centre of mass by",xrot,yrot,zrot
		for atom in FMVAR.MOLATOMS[mol]:
			coords1.append(newcoord[atom])
			types1.append(ConfSpec.ATOMTYPES[atom])
			#print ConfSpec.ATOMTYPES[atom]
			#print atomicnumber(ConfSpec.ATOMTYPES[atom])
			#print atomicmass[atomicnumber(ConfSpec.ATOMTYPES[atom])]
			mass1 = mass1 + atomicmass[atomicnumber(ConfSpec.ATOMTYPES[atom])]
		
		
		#print mass1	
		com1x=0.0
		com1y=0.0
		com1z=0.0		
		for i in range(0,len(coords1)):
			com1x=com1x+coords1[i][0]*atomicmass[atomicnumber(types1[i])]
			com1y=com1y+coords1[i][1]*atomicmass[atomicnumber(types1[i])]
			com1z=com1z+coords1[i][2]*atomicmass[atomicnumber(types1[i])]
		c_o_mass = [com1x/mass1, com1y/mass1, com1z/mass1]
		
		
		xvector=[c_o_mass[0]+1.0, c_o_mass[1], c_o_mass[2], rot[0]]
		yvector=[c_o_mass[0], c_o_mass[1]+1.0, c_o_mass[2], rot[1]]
		zvector=[c_o_mass[0], c_o_mass[1], c_o_mass[2]+1.0, rot[2]]
		rotvector=[xvector,yvector,zvector]
		
		
		for vector in rotvector:
			magvector=math.sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2])
			unitvector=[vector[0]/magvector,vector[1]/magvector,vector[2]/magvector]
			theta=vector[3]/180.0*math.pi
			
			for atom in FMVAR.MOLATOMS[mol]:
				
				dotproduct=unitvector[0]*(float(newcoord[atom][0])-float(c_o_mass[0]))+unitvector[1]*(float(newcoord[atom][1])-float(c_o_mass[1]))+unitvector[2]*(float(newcoord[atom][2])-float(c_o_mass[2]))
				centre=[float(c_o_mass[0])+dotproduct*unitvector[0],float(c_o_mass[1])+dotproduct*unitvector[1],float(c_o_mass[2])+dotproduct*unitvector[2]]
				v=[float(newcoord[atom][0])-centre[0],float(newcoord[atom][1])-centre[1],float(newcoord[atom][2])-centre[2]]
				d=math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
				px=v[0]*math.cos(theta)+v[1]*math.sin(theta)*unitvector[2]-v[2]*math.sin(theta)*unitvector[1]
				py=v[1]*math.cos(theta)+v[2]*math.sin(theta)*unitvector[0]-v[0]*math.sin(theta)*unitvector[2]
				pz=v[2]*math.cos(theta)+v[0]*math.sin(theta)*unitvector[1]-v[1]*math.sin(theta)*unitvector[0]
				newv=[px+centre[0],py+centre[1],pz+centre[2]]
				newdist=math.sqrt(px*px+py*py+pz*pz)	
				newcoord[atom]=newv


	return newcoord
	

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
		
		#print currentatom
		for i in range(0,len(geometry)):
			#print i,geometry[i]
			newcoord.append(geometry[i])
		for i in range(2,len(currentatom)-1):
			for atom in currentatom[i]:
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
					
class makemirror:
	
	def __init__(self, CONFSPEC):
		
		self.CARTESIANS = []
		for cart in CONFSPEC.CARTESIANS:
			self.CARTESIANS.append([-1.0*cart[0], cart[1], cart[2]])
		self.ATOMTYPES = CONFSPEC.ATOMTYPES
		self.CONNECTIVITY = CONFSPEC.CONNECTIVITY
		self.NAME = CONFSPEC.NAME			
	
					
class OrderConfs:
	def __init__(self, CSEARCH, SEARCHPARAMS, start, log):
		#Order the low energy conformers by energy
		self.CARTESIANS = []
		self.NAME = []
		self.TIMESFOUND = []
		self.USED = []
		self.CPU = []
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
			self.CPU.append(CSEARCH.CPU[match])
			self.TORVAL.append(CSEARCH.TORVAL[match])
		
		CSEARCH.CARTESIANS = self.CARTESIANS
		CSEARCH.NAME = self.NAME
		CSEARCH.ENERGY = self.ENERGY
		CSEARCH.TIMESFOUND = self.TIMESFOUND
		CSEARCH.USED = self.USED
		CSEARCH.CPU = self.CPU
		CSEARCH.TORVAL = self.TORVAL

class AddConformer:
	def __init__(self, CSEARCH, CONFSPEC):
		CSEARCH.NAME.append(CONFSPEC.NAME)
		CSEARCH.ENERGY.append(CONFSPEC.ENERGY)
		CSEARCH.CARTESIANS.append(CONFSPEC.CARTESIANS)						
		CSEARCH.CONNECTIVITY.append(CONFSPEC.CONNECTIVITY)
		CSEARCH.TIMESFOUND.append(1)
		CSEARCH.USED.append(0)
		CSEARCH.CPU.append(CONFSPEC.CPU)
		CSEARCH.TORVAL.append(getTorsion(CONFSPEC))
		CSEARCH.NSAVED = CSEARCH.NSAVED + 1


class RemoveConformer:
	def __init__(self, CSEARCH, todel):
		j=0
		for i in range(0,len(CSEARCH.NAME)): print CSEARCH.NAME[i], CSEARCH.ENERGY[i]
		print todel, len(todel), 
		cutoff = (len(CSEARCH.NAME)-len(todel))
		print cutoff
		#print len(CSEARCH.NAME)-len(todel)
		newtodel=[]
		for i in range(len(todel)-1, -1, -1): newtodel.append(todel[i])
		#print newtodel
		#print len(CSEARCH.NAME), CSEARCH.NAME
		print CSEARCH.NAME[cutoff:]
		print CSEARCH.NAME[:cutoff]

		for i in range(0,len(todel)):
			#print i, todel[i], CSEARCH.TIMESFOUND[todel[i]]
			CSEARCH.NREJECT = CSEARCH.NREJECT + CSEARCH.TIMESFOUND[todel[i]]

		del CSEARCH.NAME[cutoff:]
		#print len(CSEARCH.NAME), CSEARCH.NAME
		del CSEARCH.ENERGY[cutoff:]
		#print len(CSEARCH.ENERGY), CSEARCH.ENERGY
		del CSEARCH.CARTESIANS[cutoff:]
		del CSEARCH.CONNECTIVITY[cutoff:]
		del CSEARCH.USED[cutoff:]
		del CSEARCH.CPU[cutoff:]
		del CSEARCH.TIMESFOUND[cutoff:]
		del CSEARCH.TORVAL[cutoff:]
		print "AFTER REMOVAL"
		for i in range(0,len(CSEARCH.NAME)): print CSEARCH.NAME[i], CSEARCH.ENERGY[i]
	
		#print len(CSEARCH.CONNECTIVITY),len(CSEARCH.USED), len(CSEARCH.CPU), len(CSEARCH.TIMESFOUND), len(CSEARCH.TORVAL)
		CSEARCH.NSAVED = len(CSEARCH.NAME)
		

class CleanAfterJob:	
	def __init__(self, Job, Confspec, samecheck, toohigh, isomerize):					
		try:
			for suffix in [".com", ".mop", ".arc", ".temp", ".end", ".chk", ".joblog", ".csh", ".errlog"]: 
				if os.path.exists(Confspec.NAME+suffix): os.remove(Confspec.NAME+suffix)

			# If discarded remove the outfile as well		
			if isJobFinished(Job, Confspec) != 1 or samecheck > 0 or toohigh == 1 or isomerize == 1:		
				if os.path.exists(Confspec.NAME+".out"): os.remove(Confspec.NAME+".out")
		except: pass

class CleanUp:
	def __init__(self, CSearch, SearchParams, filein, log):
		# Create a tarball of low energy outfiles
		try:
			# First open previous incarnation of the tarfile
			if os.path.isfile(filein+"_fm.tgz") == 1: 
				tar = tarfile.open(filein+"_fm.tgz", "r:gz")			
				prevfiles = tar.getnames()
				for prev in prevfiles: tar.extract(prev, path="")
				os.remove(filein+"_fm.tgz")
				
				
			# Now create and write
			tar = tarfile.open(filein+"_fm.tgz", "w|gz")			
			for saved in CSEARCH.NAME: 
				#print "TARRING", saved+".out"
				if os.path.isfile(saved+".out") == 1: tar.add(saved+".out")
				#print "found and zipping", saved+".out"
				#print commands.getoutput("ls -l -t "+filein+"_fm.tgz")
				if os.path.isfile(saved+".log") == 1: tar.add(saved+".log")
				#print "found and zipping", saved+".out"
				#print commands.getoutput("ls -l -t "+filein+"_fm.tgz")
				time.sleep(0.1)
			tar.close()	
			
			
			for i in range(0,CSearch.STEP*SearchParams.POOL+1):
				if os.path.isfile(filein+"_step_"+str(i)+".out") == 1: os.remove(filein+"_step_"+str(i)+".out")
					#print "found and deleting", filein+"_step_"+str(i)+".out"
				if os.path.isfile(filein+"_step_"+str(i)+".log") == 1: os.remove(filein+"_step_"+str(i)+".log")
	

		except: print "ERROR IN ZIPPING!!!"

def reName(dir, oldname, newname):
	files=commands.getstatusoutput("ls "+dir+"/*"+oldname+"*")
	if files[0] == 0:
		for file in string.split(files[1],'\n'):
			newfile = file.replace(oldname, newname)
			print "Renaming",file,"to",newfile
			commands.getoutput("mv "+file+" "+newfile) 
	else:
		print "No files Found"

# Set an executable
class SETUPEXE:
	def __init__(self, executable):
		modify = "false"
		if os.path.exists(executable):
			var = raw_input("\no  Path to Mopac executable seem to be ok, do you really want to change this? ... ")
			if var.lower() == "y" or var.lower() == "": modify = "true"
			else: print "\nExiting";  sys.exit(1)
		else: modify = "true"
		
		if modify == "true":
			executable = ''
			enter_values = "true"
			while enter_values == "true":
				fmlocation = raw_input("\no  Enter location of FullMonte ($FULL_MONTE_DIR), including full path ... ")
				
				if os.path.exists(fmlocation):
					enter_values = "false"
				
				else: print "\no  Wrong path, please repeat the procedure"
			enter_values = "true"
			while enter_values == "true":
				executable = raw_input("\no  Enter filename of mopac executable, including full path ... ")
				
				if os.path.exists(executable):
					enter_values = "false"
				
				else: print "\no  Wrong path, please repeat the procedure"
		
		# Modify the value of MOPAC_EXEC in FMTools and then quit
		print "\no  Modifying path to Mopac exectubale ...\n"
		thisfile = open(fmlocation+"/FMTools.py","r") 
		theselines = thisfile.readlines()
		theselines[38] = "MOPAC_EXEC = \'"+executable+"\'\n"
		replacethisfile = open(fmlocation+"/FMTools.py","w")
		for line in theselines: replacethisfile.write(line)
		sys.exit()

# Define conformational search statistics 
class FMRESULTS: pass

# Names from Chang, Guida and Still's definitions
CSEARCH = FMRESULTS()
CSEARCH.NREJECT = 0
CSEARCH.NFAILED = 0
CSEARCH.AERATE = 0
CSEARCH.DMIN = 1
CSEARCH.STEP = 1	
CSEARCH.NAME = []
CSEARCH.CARTESIANS = []
CSEARCH.CONNECTIVITY = []
CSEARCH.USED = [0]
CSEARCH.TIMESFOUND = [1]
CSEARCH.NSAVED = 1
CSEARCH.TODO=[]
CSEARCH.NJOBS = 0
CSEARCH.DONE = 0
CSEARCH.RUNNING = 0
CSEARCH.DONELIST = []
CSEARCH.RUNNINGLIST = []
CSEARCH.DEAD = 0
CSEARCH.DEADLIST = []
CSEARCH.COMPLETE = 0


###############################################################
#                          FMOutput                           #
#       Formatted output to command line and log file         #
###############################################################
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
		
		# ringstring=""
        # for ringatom in Variables.RING: ringstring = ringstring+"{"+str(MolSpec.ATOMTYPES[ringatom[0]])+str(ringatom[0]+1)+"-"+str(MolSpec.ATOMTYPES[ringatom[1]])+str((ringatom[1]+1))+"} "
		
		log.Write(dashedline+"\n   |    "+("FULL_MONTE "+Params.CSEARCH+" search on "+strucname).ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("COMP: "+str(Params.COMP)+" degrees").ljust(leftcol)+("|").rjust(rightcol))
		if len(Params.FIXT) > 0:
			log.Write("   | o  "+("FIXT: Manually constrained "+str(len(Params.FIXT))+" torsional variables").ljust(leftcol)+("|").rjust(rightcol)); log.Write("   |    "+fixtstring.ljust(leftcol)+("|").rjust(rightcol))
		# if MCRI!=0:
		# print "|  ",("MCRI: There are "+str(MCRI)+" rotatable ring bonds").ljust(leftcol),("|").rjust(rightcol); print "|  ",ringstring.ljust(leftcol),("|").rjust(rightcol),"\n",emptyline
		if Variables.NMOLS > 1:
			log.Write("   | o  "+("Detected "+str(Variables.NMOLS)+" separate molecules - this adds additional search coordinates").ljust(leftcol)+("|").rjust(rightcol));
			for i in range(0,len(molarray)): 
				chunks, chunk_size = len(molarray[i]), leftcol
				for j in range(0, chunks, chunk_size):
					log.Write("   |    "+(molarray[i][j:j+chunk_size]).ljust(leftcol)+("|").rjust(rightcol))
				#log.Write("   |  "+molarray[i].ljust(leftcol)+("|").rjust(rightcol-1))
		log.Write("   | o  "+("LEVL: "+str(Params.LEVL)+" model chemistry").ljust(leftcol)+("|").rjust(rightcol))
		if Params.CSEARCH == "MCMM":
			log.Write("   | o  "+("DEMX: "+str(Params.DEMX)+" kJ/mol").ljust(leftcol)+("|").rjust(rightcol))
			if len(Params.EQUI) > 0:
				log.Write("   | o  "+("EQUI: The following sets of atoms are equivalent ").ljust(leftcol)+("|").rjust(rightcol))
				log.Write("   |    "+equistring.ljust(leftcol)+("|").rjust(rightcol))
			log.Write("   | o  "+("EWIN: "+str(Params.EWIN)+" kJ/mol").ljust(leftcol)+("|").rjust(rightcol))
			log.Write("   | o  "+("MCSS: "+Params.MCSS+" within "+str(Params.EWIN)+" kJ/mol").ljust(leftcol)+("|").rjust(rightcol))
			
		if Params.CSEARCH == "SUMM":
			log.Write("   | o  "+("ITVL: "+str(Params.ITVL)+" degrees").ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("MCNV: "+str(Variables.MCNV)).ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   |    "+torstring.ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("PROC: "+str(Params.POOL)+" processors will be used").ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("RJCT: "+str(Params.RJCT)+"x sum of Bondi radii").ljust(leftcol)+("|").rjust(rightcol))
		log.Write("   | o  "+("STEP: "+str(Params.STEP)).ljust(leftcol)+("|").rjust(rightcol))
		#log.Write("|  "+("POOL: pool size for simultaneous optimizations: "+str(Params.POOL)).ljust(leftcol)+("|").rjust(rightcol)+"\n"+emptyline)
		if Params.HSWAP>0:
			log.Write("   | "+("HSWAP: "+str(Params.HSWAP)+" meaning NH,OH and SH protons are allowed to wander...").ljust(leftcol)+("|").rjust(rightcol)+"\n"+emptyline)
		log.Write(dashedline+"\n")


class WriteSummary:
	# Formatted text printed to terminal and log file at the end of each search step
	def __init__(self, CSearch, SearchParams, start, log):
		
		now = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
		totalCPU = CPUTime(CSearch)
		runningtime = RealTime(start, now)
		
		if CSearch.COMPLETE == 0: log.Write("\no  "+("STEP "+str(CSearch.STEP)+" COMPLETE: "+str(CSearch.NSAVED)+" unique conformations. Global minimum energy = "+str(round(CSearch.GLOBMIN,5))).ljust(leftcol)+("").rjust(rightcol))
		if CSearch.COMPLETE == 1: log.Write("\no  "+("FULL MONTE SEARCH COMPLETE: "+str(CSearch.NSAVED)+" unique conformations. Global minimum energy = "+str(round(CSearch.GLOBMIN,5))).ljust(leftcol)+("").rjust(rightcol))
		
		log.Write(("\n     Conformer Name             Absolute Energy       Erel (kJ/mol)        Times found    Times used         CPU").ljust(leftcol)+"\n"+dashedline)
		
		for i in range(0, CSearch.NSAVED):
			absenergy = str(round(float(CSearch.ENERGY[i]),5))
			if len(absenergy.split(".")[1])!=5: absenergy = absenergy+"0"
			relenergy = str(round(float((CSearch.ENERGY[i]-CSearch.GLOBMIN)*2625.5),2))
			if len(relenergy.split(".")[1])!=2: relenergy = relenergy+"0"
			log.Write("   | "+CSearch.NAME[i].ljust(30)+(absenergy).ljust(20)+(relenergy).rjust(10)+ (str(CSearch.TIMESFOUND[i])).rjust(15)+(str(CSearch.USED[i])).rjust(15)+(str(CSearch.CPU[i][0])+"d "+str(CSearch.CPU[i][1])+"h "+str(CSearch.CPU[i][2])+"m "+str(CSearch.CPU[i][3])+"s").rjust(20)+("|").rjust(2))
		log.Write(dashedline+"\n   | o  "+("CPU time: "+str(totalCPU[0])+"d "+str(totalCPU[1])+"h "+str(totalCPU[2])+"m "+str(totalCPU[3])+"s   Execute time: "+str(runningtime[0])+"d "+str(runningtime[1])+"h "+str(runningtime[2])+"m "+str(runningtime[3])+"s").ljust(leftcol)+("|").rjust(rightcol))
		if CSearch.COMPLETE == 0: log.Write("   | o  "+("SE = "+str(round(float(CSearch.AERATE),1))+"    DMIN = "+str(CSearch.DMIN)+"    NOPT = "+str(CSearch.STEP*SearchParams.POOL)+"    NFAIL = "+str(CSearch.NFAILED)).ljust(leftcol)+("|").rjust(rightcol)+"\n"+ dashedline)
		if CSearch.COMPLETE == 1: log.Write("   | o  "+("SE = "+str(round(float(CSearch.AERATE),1))+"    DMIN = "+str(CSearch.DMIN)+"    NOPT = "+str((CSearch.STEP-1)*SearchParams.POOL)+"    NFAIL = "+str(CSearch.NFAILED)).ljust(leftcol)+("|").rjust(rightcol)+"\n"+ dashedline)


class WriteQMSummary:
	
	# Formatted text printed to terminal and log file at the end of each search step
	def __init__(self, CSearch, OldSearch, SearchParams, start, log, qmsuffix):
		
		now = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
		runningtime = RealTime(start, now)
		totalCPU = CPUTime(CSearch)
		
		if CSearch.COMPLETE == 0: log.Write("\no  "+(str(CSearch.DONE)+" QM JOBS COMPLETED: "+str(CSearch.NSAVED)+" unique conformations. Global minimum energy = "+str(round(CSearch.GLOBMIN,5))).ljust(leftcol)+("").rjust(rightcol))
		if CSearch.COMPLETE == 1: log.Write("\no  "+("FULL MONTE SEARCH COMPLETE: "+str(CSearch.NSAVED)+" unique conformations. Global minimum energy = "+str(round(CSearch.GLOBMIN,5))).ljust(leftcol)+("").rjust(rightcol))
		
		log.Write("\n     "+"Conformer Name".ljust(50)+"Eabs (H)".ljust(12)+"Erel (kJ/mol)".ljust(16)+"X found".ljust(14)+"CPU".ljust(10)+"S.E Erel".ljust(28)+(" ")+"\n"+dashedline)
		
		for i in range(0, CSearch.NSAVED):
			absenergy = str(round(float(CSearch.ENERGY[i]),5))
			if len(absenergy.split(".")[1])!=5: absenergy = absenergy+"0"
			relenergy = str(round(float((CSearch.ENERGY[i]-CSearch.GLOBMIN)*2625.5),2))
			if len(relenergy.split(".")[1])!=2: relenergy = relenergy+"0"
			for j in range(0,len(OldSearch.NAME)):
				if (OldSearch.NAME[j]+"_"+qmsuffix).find(CSearch.NAME[i]) > -1:
					olderel = str(round(float((OldSearch.ENERGY[j]-OldSearch.GLOBMIN)*2625.5),2))
			
			log.Write("   | "+CSearch.NAME[i].ljust(50)+(absenergy).ljust(10)+(relenergy).rjust(10)+ (str(CSearch.TIMESFOUND[i])).rjust(10)+(str(CSearch.CPU[i][0])+"d "+str(CSearch.CPU[i][1])+"h "+str(CSearch.CPU[i][2])+"m "+str(CSearch.CPU[i][3])+"s").rjust(20)+olderel.rjust(10)+("|").rjust(1))
		log.Write(dashedline)
		log.Write("   | o  "+("CPU time: "+str(totalCPU[0])+"d "+str(totalCPU[1])+"h "+str(totalCPU[2])+"m "+str(totalCPU[3])+"s   Execute time: "+str(runningtime[0])+"d "+str(runningtime[1])+"h "+str(runningtime[2])+"m "+str(runningtime[3])+"s").ljust(leftcol)+("|").rjust(rightcol))
		log.Write(dashedline)


class makeGVformat:
	#Write a file for viewing in Gaussview that contains the low energy conformations in ascending order of energy
	def __init__(self, filein, MolSpec, CSearch, Params, append):
		
		#The following formatting is necessary to make GaussView believe the conformers are the steps of a Gaussian optimization
		gradline = " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad	"
		soline = "                         Standard orientation:                         "
		dashline1 = " ---------------------------------------------------------------------"
		dashline2 = "  ---------------------------------------------------------------------"
		coordline1 = " Center     Atomic     Atomic              Coordinates (Angstroms)"
		coordline2 = " Number     Number      Type              X           Y           Z"
		quoteline1 = " This file enables the use of Gaussview to look at the saved conformers from a FullMonte Search (Select Read Intermediate Geometries when opening...)"
		quoteline2 = " The conformers are ordered by energy, with the most stable first"
		endline = " Normal termination of Gaussian XX."
		
		if os.path.isfile(filein+"_"+append+".out") == 1: os.remove(filein+"_"+append+".out")
		outfile = open(filein+"_"+append+".out", 'w' )
		
		#summarize the FullMonte Search - this info is read to construct QM jobs in the next stage
		formtime = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
		outfile.write("o  FULLMONTE Search Performed on "+filein+" at "+formtime+"\n")
		outfile.write("o  Search Parameters Used ...\n")
		outfile.write("COMP = "+str(Params.COMP)+"\n")
		#outfile.write("o  FIXT = ")
		outfile.write("DEMX = "+str(Params.DEMX)+"\n")
		for equivcoords in Params.EQUI:
			outfile.write("EQUI = "+equivcoords+"\n")
		outfile.write("EWIN = "+str(Params.EWIN)+"\n")
		outfile.write("LEVL = "+str(Params.EWIN)+"\n")
		outfile.write("MCSS = "+str(Params.MCSS)+"\n")
		outfile.write("PROC = "+str(Params.POOL)+"\n")
		outfile.write("RJCT = "+str(Params.RJCT)+"\n")
		outfile.write("STEP = "+str(Params.STEP)+"\n\n")
		
		
		if CSearch.NSAVED > 0:
			outfile.write(gradline+"\n"+gradline+"\n"+soline+"\n"+dashline1+"\n"+coordline1+"\n"+coordline2+"\n"+dashline1+"\n")
			for i in range(0, MolSpec.NATOMS):
				x = "%.6f" % CSearch.CARTESIANS[0][i][0]
				y = "%.6f" % CSearch.CARTESIANS[0][i][1]
				z = "%.6f" % CSearch.CARTESIANS[0][i][2]
				outfile.write(str(i+1).rjust(5)+str(atomicnumber(MolSpec.ATOMTYPES[i])).rjust(11)+"0".rjust(14)+x.rjust(16)+y.rjust(12)+z.rjust(12)+"\n")
			outfile.write(dashline2+"\n\n"+gradline+"\n"+gradline+"\n\n")
			for i in range(0, CSearch.NSAVED-1):
				outfile.write("CONFORMER "+str(i+1)+": "+CSearch.NAME[i]+"  "+str(CSearch.ENERGY[i])+"   "+str(CSearch.TIMESFOUND[i])+"\n")
				outfile.write(soline+"\n"+dashline1+"\n"+coordline1+"\n"+coordline2+"\n"+dashline1+"\n")
				for j in range(0, MolSpec.NATOMS):
					x = "%.6f" % CSearch.CARTESIANS[i][j][0]
					y = "%.6f" % CSearch.CARTESIANS[i][j][1]
					z = "%.6f" % CSearch.CARTESIANS[i][j][2]
					outfile.write(str(j+1).rjust(5)+str(atomicnumber(MolSpec.ATOMTYPES[j])).rjust(11)+"0".rjust(14)+x.rjust(16)+y.rjust(12)+z.rjust(12)+"\n")
				outfile.write(dashline2+"\n\n"+gradline+"\n Step number"+str(i+1).rjust(4)+" out of a maximum of"+str(CSearch.NSAVED).rjust(4)+"\n"+gradline+"\n\n")
			outfile.write("CONFORMER "+str(CSearch.NSAVED)+": "+CSearch.NAME[CSearch.NSAVED-1]+"  "+str(CSearch.ENERGY[CSearch.NSAVED-1])+"   "+str(CSearch.TIMESFOUND[CSearch.NSAVED-1])+"\n")
			outfile.write(gradline+"\n Step number"+str(CSearch.NSAVED).rjust(4)+" out of a maximum of"+str(CSearch.NSAVED).rjust(4)+"\n"+gradline+"\n\n"+soline+"\n"+dashline1+"\n"+coordline1+"\n"+coordline2+"\n"+dashline1+"\n")
			for i in range(0, MolSpec.NATOMS):
				x = "%.6f" % CSearch.CARTESIANS[CSearch.NSAVED-1][i][0]
				y = "%.6f" % CSearch.CARTESIANS[CSearch.NSAVED-1][i][1]
				z = "%.6f" % CSearch.CARTESIANS[CSearch.NSAVED-1][i][2]
				outfile.write(str(i+1).rjust(5)+str(atomicnumber(MolSpec.ATOMTYPES[i])).rjust(11)+"0".rjust(14)+x.rjust(16)+y.rjust(12)+z.rjust(12)+"\n")
			outfile.write(dashline2+"\n"+quoteline1+"\n"+quoteline2+"\n"+endline+"\n")


class makePDBformat:
	#Write a file for viewing in Gaussview that contains the low energy conformations in ascending order of energy
	def __init__(self, filein, MolSpec, CSearch,append):
		pdbfile = open(filein+"_"+append+".pdb", 'w' )
		if CSearch.NSAVED > 0:
			for i in range(0, CSearch.NSAVED):
				#pdbfile = open(CSearch.NAME[i]+"_fm.pdb", 'w' )
				Erel = (CSearch.ENERGY[i]-CSearch.GLOBMIN)*2625.5
				pdbfile.write("COMPND    "+CSearch.NAME[i]+"   E = "+str(Erel)+"\nAUTHOR    FULL MONTE SEARCH - www.patonlab.com")
				for j in range(0, MolSpec.NATOMS):
					x = "%.3f" % CSearch.CARTESIANS[i][j][0]
					y = "%.3f" % CSearch.CARTESIANS[i][j][1]
					z = "%.3f" % CSearch.CARTESIANS[i][j][2]
					pdbfile.write("\nHETATM"+str(j+1).rjust(5)+MolSpec.ATOMTYPES[j].rjust(3)+"   LIG     1".rjust(10)+x.rjust(12)+y.rjust(8)+z.rjust(8))
				pdbfile.write("\nEND    ")
			pdbfile.close()


# Formatting
dashedline = "   ------------------------------------------------------------------------------------------------------------------"
emptyline = "   |                                                                                                                     |"
normaltermination = "\n   -----------------       N   O   R   M   A   L      T   E   R   M   I   N   A   T   I   O   N      ----------------\n"
leftcol=97
rightcol=12


# Ascii Art
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


#Read molecule data from an input file - currently Gaussian *.com and *.pdb supported
class getinData:
	
	def __init__(self, file,log):
		
		def getJOBTYPE(self, inlines):
			if fileformat == ".com":
				for i in range(0,len(inlines)):
					if inlines[i].find("#") > -1: self.JOBTYPE = inlines[i]
		
		
		def getCHARGE(self, inlines):
			if fileformat == ".com":
				for i in range(0,len(inlines)):
					if inlines[i].find("#") > -1:
						if len(inlines[i+1].split()) == 0:
							self.CHARGE = inlines[i+4].split()[0]
							self.MULT = inlines[i+4].split()[1]
						if len(inlines[i+2].split()) == 0:
							self.CHARGE = inlines[i+5].split()[0]
							self.MULT = inlines[i+5].split()[1]
		
		def getMEMREQ(self, inlines):
			if fileformat == ".com":
				for i in range(0,len(inlines)):
					if inlines[i].find("%mem") > -1: self.MEMREQ = inlines[i].split("=")[1].rstrip("\n")
		
		
		def getNPROC(self, inlines):
			if fileformat == ".com":
				for i in range(0,len(inlines)):
					if inlines[i].find("%nproc") > -1: self.NPROC = inlines[i].split("=")[1].rstrip("\n")
		
		
		def getATOMTYPES(self, inlines):
			if fileformat == ".com":
				for i in range(0,len(inlines)):
					if inlines[i].find("#") > -1:
						if len(inlines[i+1].split()) == 0: start = i+5
						if len(inlines[i+2].split()) == 0: start = i+6
						break
				
				self.ATOMTYPES = []
				self.LEVELTYPES = []
				self.MMTYPES = []
				for i in range(start,len(inlines)):
					if len(inlines[i].split()) ==0: break
					else:
						atominfo = inlines[i].split()[0]
						atominfo = atominfo.split("-")[0]
						if len(inlines[i].split()[0].split(atominfo))>1:
							mminfo = inlines[i].split()[0].lstrip(atominfo)
							self.MMTYPES.append(mminfo)
						
						self.ATOMTYPES.append(atominfo)
						level = ""
						for oniomlevel in ["H", "M", "L"]:
							if inlines[i][4:].rfind(oniomlevel)>1:
								level = inlines[i][4:][inlines[i][4:].rfind(oniomlevel):].rstrip("\n")
						self.LEVELTYPES.append(level)
		
		
		def getCONNECTIVITY(self, inlines, natoms):
			if fileformat == ".com":
				for i in range(0,len(inlines)):
					if inlines[i].find("#") > -1:
						if len(inlines[i+1].split()) == 0: start = i+natoms+6
						if len(inlines[i+2].split()) == 0: start = i+natoms+7
						break
				
				if start < len(inlines):
					self.CONNECTIVITY = []
					j = 1
					for i in range(start,len(inlines)):
						if len(inlines[i].split()) != 0:
							try: num = int(inlines[i].split()[0])
							except ValueError: num = 0
						if num == j:
							bond=[]
							neighbors=(len(inlines[i].split())-1)/2
							if neighbors!=0:
								for k in range(0,neighbors): bond.append((inlines[i].split()[1+2*k])+"__"+(inlines[i].split()[2+2*k]))
							self.CONNECTIVITY.append(bond)
							j = j+1
					
					if len(self.CONNECTIVITY) == natoms:
						for i in range(0, natoms):
							for partner in self.CONNECTIVITY[i]:
								info = partner.split("__")
								nextatom = int(info[0])-1
								bondorder = float(info[1])
								nope=0
								for otherpartner in self.CONNECTIVITY[nextatom]:
									otherinfo = otherpartner.split("__")
									othernextatom = int(otherinfo[0])-1
									if othernextatom==i: nope=nope+1
								if nope==0: self.CONNECTIVITY[nextatom].append(str(i+1)+"__"+info[1])
					
					self.OPTIONAL = []
					for i in range(start+j,len(inlines)):
						if len(inlines[i].split()) != 0: self.OPTIONAL.append(inlines[i])
		
		
		def getCARTESIANS(self, inlines, natoms):
			if fileformat == ".com":
				for i in range(0,len(inlines)):
					if inlines[i].find("#") > -1:
						if len(inlines[i+1].split()) == 0: start = i+5
						if len(inlines[i+2].split()) == 0: start = i+6
						break
				
				self.CARTESIANS = []
				for i in range(start,len(inlines)):
					if len(inlines[i].split()) == 0: break
					elif len(inlines[i].split()) == 4: self.CARTESIANS.append([float(inlines[i].split()[1]), float(inlines[i].split()[2]), float(inlines[i].split()[3])])
					elif len(inlines[i].split()) > 4: self.CARTESIANS.append([float(inlines[i].split()[2]), float(inlines[i].split()[3]), float(inlines[i].split()[4])])
		
		def getCONSTRAINED(self, optional):
			if fileformat == ".com":
				self.CONSTRAINED = []
				for line in optional:
					if line.find("X") > -1 and line.find("F") > -1: self.CONSTRAINED.append([int(line.split(" ")[1])-1])
					if line.find("B") > -1 and line.find("F") > -1: self.CONSTRAINED.append([int(line.split(" ")[1])-1,int(line.split(" ")[2])-1])
					if line.find("A") > -1 and line.find("F") > -1: self.CONSTRAINED.append([int(line.split(" ")[1])-1,int(line.split(" ")[2])-1]),int(line.split(" ")[3])-1
					if line.find("D") > -1 and line.find("F") > -1: self.CONSTRAINED.append([int(line.split(" ")[1])-1,int(line.split(" ")[2])-1, int(line.split(" ")[3])-1, int(line.split(" ")[4])-1])
		
		for suffix in [".com", ".pdb"]:
			if os.path.exists(file+suffix): fileformat = suffix

		try: fileformat
		except NameError: log.Fatal("\nFATAL ERROR: Input file [ %s ] does not exist"%file)
		
		infile = open(file+fileformat,"r")
		inlines = infile.readlines()
		self.NAME = file
		getJOBTYPE(self, inlines)
		getCHARGE(self, inlines)
		getMEMREQ(self, inlines)
		getNPROC(self, inlines)
		getATOMTYPES(self, inlines)
		self.NATOMS=len(self.ATOMTYPES)
		getCONNECTIVITY(self, inlines, self.NATOMS)
		getCARTESIANS(self, inlines, self.NATOMS)
		if hasattr(self, "OPTIONAL"): getCONSTRAINED(self, self.OPTIONAL)
		if len(self.ATOMTYPES) == 0 or len(self.CARTESIANS) ==0: log.Fatal("\nFATAL ERROR: Input file [ %s ] cannot be read"%file)


# Reads parameters specified in the params file to be used in the conformational search

class getParams:
	
	def __init__(self, MolSpec, instruct, log):
		#Default Parameters
		self.CSEARCH = "MCMM"
		self.MCSS="Uniform Usage Directed"
		self.ITVL=60
		self.EWIN=20.0
		self.COMP=10.0
		self.DEMX=41.84
		self.RJCT=0.50
		self.POOL=1
		self.FIXT=[]
		self.EQUI=[]
		self.ETOZ=[]
		self.HSWAP=0
		self.STEP=0
		self.LEVL="PM6"
		self.NNBO=0
		
		#The default values are replaced by those in instruct file (if it exists)
		if instruct != "default":
			if not os.path.exists(instruct): log.Fatal("\nFATAL ERROR: Full Monte parameter file [ %s ] does not exist"%instruct)
			
			instructfile = open(instruct,"r")
			instructlines = instructfile.readlines()
			
			
			for line in instructlines:
				if line.find("SUMM") > -1: self.CSEARCH = "SUMM"
				if line.find("ITVL") > -1: self.ITVL = int(line.split("=")[1].rstrip('\n').lstrip())
				if line.find("LEVL") > -1: self.LEVL = line.split("=")[1].rstrip('\n').lstrip()
				if line.find("NNBO") > -1: self.NNBO = int(line.split("=")[1].rstrip('\n').lstrip())
				if line.find("STEP") > -1: self.STEP = int(line.split("=")[1].rstrip('\n').lstrip())
				if line.find("PROC") > -1: self.POOL = int(line.split("=")[1].rstrip('\n').lstrip())
				if line.find("RJCT") > -1: self.RJCT= float(line.split("=")[1].rstrip('\n').lstrip())
				if line.find("COMP") > -1: self.COMP = float(line.split("=")[1].rstrip('\n').lstrip())
				if line.find("EWIN") > -1: self.EWIN = float(line.split("=")[1].rstrip('\n').lstrip())
				if line.find("DEMX") > -1: self.DEMX = float(line.split("=")[1].rstrip('\n').lstrip())
				if 	line.find("FIXT") > -1: self.FIXT.append((line.split("=")[1]).rstrip('\n').lstrip())
				if      line.find("ETOZ") > -1: self.ETOZ.append((line.split("=")[1]).rstrip('\n').lstrip())
				if 	line.find("EQUI") > -1: self.EQUI.append((line.split("=")[1]).rstrip('\n').lstrip())
				if 	line.find("HSWP") > -1: self.HSWAP = int((line.split("=")[1]).rstrip('\n').lstrip())
		
		
		# Check that specified equivalency makes sense
		for equi in self.EQUI:
			atomstring=[]
			equivatoms=equi.split(" ")
			for atom in equivatoms:
				if int(atom) > MolSpec.NATOMS: log.Fatal("\nFATAL ERROR: Full Monte parameter file [ %s ] equivalent atom numbers"%instruct)
				atomstring.append(MolSpec.ATOMTYPES[int(atom)-1])
			
			for i in range(1,len(atomstring)):
				if atomstring[i]!=atomstring[i-1]: log.Fatal("\nFATAL ERROR: Full Monte parameter file [ %s ] equivalent atoms are different elements"%instruct)

		# Check that specified fixed atoms make sense
		self.FIXEDATOMS=[]
		for fix in self.FIXT:
			bondexists=0
			self.FIXEDATOMS.append([int(fix.split(" ")[0]), int(fix.split(" ")[1])])
			for bond in MolSpec.CONNECTIVITY[int(fix.split(" ")[0])-1]:
				partner  = int(bond.split("__")[0])
				if partner == int(fix.split(" ")[1]): bondexists = bondexists+1
			if bondexists == 0: log.Fatal("\nFATAL ERROR: Full Monte parameter file [ %s ] fixed torsion not bonded"%instruct)


#Read molecule data from an output file #######################
class getoutData:
	def __init__(self, MolSpec):
		
		self.NAME = MolSpec.NAME
				
		def getFORMAT(self, outlines):
			for i in range(0,len(outlines)):
				if outlines[i].find("MOPAC") > -1: self.FORMAT = "Mopac"; break
				if outlines[i].find("Gaussian") > -1: self.FORMAT = "Gaussian"; break
		
		def getCHARGE(self, outlines, format):
			if format == "Mopac":
				for i in range(0,len(outlines)):
					if outlines[i].find("CHARGE ON SYSTEM") > -1:
						self.CHARGE = int(outlines[i].split()[5])
						self.MULT = 1
					if outlines[i].find("STATE CALCULATION") > -1:
						if outlines[i].split()[0] == "SINGLET": self.MULT = 1
						if outlines[i].split()[0] == "DOUBLET": self.MULT = 2
						if outlines[i].split()[0] == "TRIPLET": self.MULT = 3
						break
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find("Charge = ") > -1:
						self.CHARGE = int(outlines[i].split()[2])
						self.MULT = int(outlines[i].split()[5].rstrip("\n"))
						break
		
		
		def getATOMTYPES(self, outlines, format):
			self.ATOMTYPES = []
			self.CARTESIANS = []
			if format == "Mopac":
				for i in range(0,len(outlines)):
					if outlines[i].find("CHEMICAL") > -1: standor = i+3
					if outlines[i].find("Empirical Formula") > -1:
						self.NATOMS = int((outlines[i].split("="))[1].split()[0])
				
				if hasattr(self, "NATOMS"):
					for i in range (standor,standor+self.NATOMS):
						outlines[i] = outlines[i].replace("*", " ")
						self.ATOMTYPES.append((outlines[i].split()[1]))
						self.CARTESIANS.append([float(outlines[i].split()[2]), float(outlines[i].split()[3]), float(outlines[i].split()[4])])
			
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find("Standard orientation") > -1:
						standor = i
					if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1:
						self.NATOMS = i-standor-6
				try: standor
				except NameError: pass
				else:
					for i in range (standor+5,standor+5+self.NATOMS):
						self.ATOMTYPES.append(elementID(int(outlines[i].split()[1])))
						self.CARTESIANS.append([float(outlines[i].split()[3]),float(outlines[i].split()[4]),float(outlines[i].split()[5])])
		
		
		def getFREQS(self, outlines, format):
			self.FREQS = []
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find("Frequencies") > -1:
						self.FREQS.append(float(outlines[i].split()[2]))
						if len(outlines[i].split()) > 3: self.FREQS.append(float(outlines[i].split()[3]))
						if len(outlines[i].split()) > 4: self.FREQS.append(float(outlines[i].split()[4]))
				if len(self.FREQS) > 0:
					for i in range(0,len(outlines)):
						if outlines[i].find("Zero-point correction") > -1: self.ZPE = float(outlines[i].split()[2])
						if outlines[i].find("thermal Enthalpies") > -1: self.ENTHALPY = float(outlines[i].split()[6])
						if outlines[i].find("thermal Free Energies") > -1: self.GIBBS = float(outlines[i].split()[7])
		
		def getCPU(self, outlines, format):
			days = 0; hours = 0; mins = 0; secs = 0
			if format == "Mopac":
				for i in range(0,len(outlines)):
					if outlines[i].find("TOTAL CPU TIME") > -1:
						secs = secs + int(float(outlines[i].split()[3]))
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find("Job cpu time") > -1:
						days = days + int(outlines[i].split()[3])
						hours = hours + int(outlines[i].split()[5])
						mins = mins + int(outlines[i].split()[7])
						secs = secs + int(float(outlines[i].split()[9]))
			self.CPU=[days,hours,mins,secs]
		
		
		
		def getENERGY(self, outlines, format):
			if format == "Mopac":
				for i in range(0,len(outlines)):
					if outlines[i].find("TOTAL ENERGY") > -1:
						self.ENERGY = 0.036749309*(float(outlines[i].split()[3]))
			if format == "Gaussian":
				uff = 0
				am1 = 0
				pm3 = 0
				scf = 0
				oniom = 0
				for i in range(0,len(outlines)):
					if outlines[i].find(" UFF") > -1: uff = i
					if outlines[i] .find("AM1") > -1: am1 = i
					if outlines[i].find("PM3") > -1: pm3 = i
					if outlines[i].find("ONIOM") > -1: oniom = i
					if outlines[i].find("SCF Done") > -1: scf = i
				
				calctype = [uff,am1,pm3,oniom,scf]
				for i in range(0,len(outlines)):
					if scf == max(calctype) and outlines[i].find("SCF Done") > -1 and outlines[i].find("Initial convergence to 1.0D-05 achieved")==-1: # Get energy from HF or DFT calculation
						self.ENERGY = (float(outlines[i].split()[4]))
					if oniom == max(calctype) and outlines[i].find("ONIOM: extrapolated energy") > -1: # Get energy from ONIOM calculation
						self.ENERGY = (float(outlines[i].split()[4]))
					if pm3 == max(calctype) or am1 == max(calctype) or uff == max(calctype):
						if outlines[i].find("Energy= ") > -1 and outlines[i].find("Predicted")==-1 and outlines[i].find("Thermal")==-1: # Get energy from Semi-empirical or Molecular Mechanics calculation
							self.ENERGY = (float(outlines[i].split()[1]))
					if outlines[i].find("Total free energy in solution") > -1:
						self.SOLVENERGY = (float(outlines[i+1].split()[7]))

		
		if not os.path.exists(self.NAME+".out") and not os.path.exists(self.NAME+".log"):
			print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)

		if os.path.exists(self.NAME+".out"):
			outfile = open(self.NAME+".out","r")
			outlines = outfile.readlines()
		if os.path.exists(self.NAME+".log"):
			outfile = open(self.NAME+".log","r");
			outlines = outfile.readlines()

		
		getFORMAT(self, outlines)
		if hasattr(self, "FORMAT"):
			getCHARGE(self, outlines, self.FORMAT)
			getENERGY(self, outlines, self.FORMAT)
			getATOMTYPES(self, outlines, self.FORMAT)
			getFREQS(self, outlines, self.FORMAT)
			getCPU(self, outlines, self.FORMAT)
###############################################################

#Create Gaussian com file or Mopac mop file
class writeInput:
	
	def __init__(self, JobSpec, MolSpec):
		
		def writeGaussian(self, JobSpec, MolSpec):
			fileout = open(MolSpec.NAME+".com", "w")
			fileout.write("%chk="+MolSpec.NAME+".chk\n")
			
			if hasattr(JobSpec, "Mem"): fileout.write("%mem = "+JobSpec.Mem+"\n")
			if hasattr(JobSpec, "NPROC"): 
				fileout.write("%nproc = "+str(JobSpec.NPROC)+"\n")
	
			for route in JobSpec.JOBTYPE.split(): fileout.write(route+" ")
			fileout.write("\n\n")
			fileout.write(MolSpec.NAME+"\n\n")
			fileout.write(str(MolSpec.CHARGE)+" "+str(MolSpec.MULT)+"\n")
			
			#print MolSpec.MMTYPES
			for i in range(0,MolSpec.NATOMS):
				fileout.write(" "+MolSpec.ATOMTYPES[i])
				if hasattr(MolSpec, "MMTYPES"): fileout.write(MolSpec.MMTYPES[i])
				for j in range(0,3):
					if math.fabs(MolSpec.CARTESIANS[i][j]) > 0.0001: fileout.write("  "+str(round(MolSpec.CARTESIANS[i][j],6)))
					else: fileout.write("  "+str(0.000000))
				if hasattr(MolSpec,"LEVELTYPES"):
					if len(MolSpec.LEVELTYPES) != 0: fileout.write("  "+MolSpec.LEVELTYPES[i])
				fileout.write("\n")
			fileout.write("\n")
			
			#if len(JobSpec.Optional) > 0:
			#	for option in JobSpec.Optional: fileout.write(option+"\n")
			#	fileout.write("\n")
			
			for i in range(0,MolSpec.NATOMS):
				fileout.write(" "+str(i+1)+" ")
				for connect in MolSpec.CONNECTIVITY[i]:
					fileout.write(connect.split("__")[0]+" "+connect.split("__")[1]+" ")
				fileout.write("\n")
			fileout.write("\n")
			
			if len(JobSpec.CONSTRAINED) > 0:
				for const in JobSpec.CONSTRAINED:
					if len(const) == 1: fileout.write("X "+str(const[0]+1)+" F\n")
					if len(const) == 2: fileout.write("B "+str(const[0]+1)+" "+str(const[1]+1)+" F\n")
				fileout.write("\n")
			
			if hasattr(JobSpec, "Link1"):
				fileout.write("--Link1--\n%chk="+JobSpec.Link0+"\n")
				if hasattr(JobSpec, "Mem"): fileout.write("%mem = "+JobSpec.Mem+"\n")
				if hasattr(JobSpec, "Nproc"): fileout.write("%nproc = "+JobSpec.Nproc+"\n")
				fileout.write("# "+JobSpec.Link1+"\n\n")
				fileout.write(JobSpec.Title+"\n\n")
				fileout.write(str(MolSpec.CHARGE)+" "+str(MolSpec.MULT)+"\n")
				if JobSpec.Link1.find("geom=check")==-1:
					for i in range(0,MolSpec.NATOMS):
						fileout.write(MolSpec.ATOMTYPES[i])
						for j in range(0,3): fileout.write("  "+str(round(MolSpec.CARTESIANS[i][j],6)))
						fileout.write("\n")
					fileout.write("\n")
				else: fileout.write("\n")
			
			if hasattr(JobSpec, "Radii"): fileout.write("radii="+JobSpec.Radii+"\n\n")

		
		def writeMopac(self, JobSpec, MolSpec):
			#print "writing MOPAC"
			fileout = open(MolSpec.NAME+".mop", "w")
			fileout.write(JobSpec.JOBTYPE+" charge="+str(MolSpec.CHARGE))
			
			if MolSpec.MULT == 1: fileout.write(" Singlet ")
			elif MolSpec.MULT == 2: fileout.write(" Doublet ")
			elif MolSpec.MULT == 3: fileout.write(" Triplet ")
			elif MolSpec.MULT == 4: fileout.write(" Quartet ")
			elif MolSpec.MULT == 5: fileout.write(" Quintet ")
			
			fileout.write("\n"+MolSpec.NAME+"\n\n")
			
			freezeindex = []
			for i in range(0,MolSpec.NATOMS):
				if MolSpec.ATOMTYPES[i].find("-"): MolSpec.ATOMTYPES[i]=MolSpec.ATOMTYPES[i].split("-")[0]
				fileout.write(MolSpec.ATOMTYPES[i])
				
				cart = 1; freezecoord = 0; freezedist = 0; freezetorsion = 0; freezeangle = 0
				for const in JobSpec.CONSTRAINED:
					#print const
					if len(const) == 1:
						freezecoord = 1; freezeindex.append(const)
					if len(const) == 2:
						if const[0] == i and const[1] < i: cart = 0; freezedist = 1; internal = const[1]; j = const[0]
						if const[1] == i and const[0] < i: cart = 0; freezedist = 1; internal = const[1]; j = const[0]
					if len(const) == 3:
						if const[0] == i and const[2] < i: cart = 0; freezeangle = 1; internal = const[2]; j = const[1]; k = const[0]
						if const[2] == i and const[0] < i: cart = 0; freezeangle = 1; internal = const[2]; j = const[1]; k = const[0]
					if len(const) == 4:
						if const[0] == i and const[3] < i: cart = 0; freezetorsion = 1; internal = const[3]; j = const[2]; k = const[1]; l = const[0]
						if const[3] == i and const[0] < i: cart = 0; freezetorsion = 1; internal = const[3]; j = const[2]; k = const[1]; l = const[0]
				
				#For all unconstrained coordinates print out Cartesian coordinates
				if cart == 1:
					optmz = "1"
					#print freezeindex
					for frozenatom in freezeindex:
						#print frozenatom
						if frozenatom[0] == i: optmz = "0"
					for n in range(0,3): fileout.write("  "+str(round(MolSpec.CARTESIANS[i][n],6))+"  "+optmz)
				
				#For all constrained coordinates define Internal coordinates
				else:
					#print "constraining", i, j
					if freezedist == 1:
						stop = 0
						for connectj in MolSpec.CONNECTIVITY[j]:
							k = int(connectj.split("__")[0])-1
							if k < i and stop == 0:
								
								for connectk in MolSpec.CONNECTIVITY[k]:
									l = int(connectk.split("__")[0])-1
									if l < i and l != j and stop == 0:
										
										fileout.write(" "+str(round(calcdist(j,i,MolSpec.CARTESIANS),3))+" 0 "+str(round(calcangle(i,j,k,MolSpec.CARTESIANS),3))+" 1 "+str(round(calcdihedral(i,j,k,l,MolSpec.CARTESIANS),3))+" 1  "+str(j+1)+" "+str(k+1)+" "+str(l+1))
										stop = 1
						
						if stop == 0:
							for connecti in MolSpec.CONNECTIVITY[i]:
								k = int(connecti.split("__")[0])-1
								if k < j and stop == 0:
									
									for connectk in MolSpec.CONNECTIVITY[k]:
										l = int(connectk.split("__")[0])-1
										if l != i and l < j and stop == 0:
											
											fileout.write(" "+str(round(calcdist(j,i,MolSpec.CARTESIANS),3))+" 0 "+str(round(calcangle(i,j,k,MolSpec.CARTESIANS),3))+" 1 "+str(round(calcdihedral(i,j,k,l,MolSpec.CARTESIANS),3))+" 1  "+str(j+1)+" "+str(k+1)+" "+str(l+1))
											stop = 1
							
							if stop == 0:
								print ("\nFATAL ERROR: Imposing constraints - try alternative atom order?")
								sys.exit()
					
					if freezetorsion == 1:
						fileout.write(" "+str(round(calcdist(j,i,MolSpec.CARTESIANS),3))+" 1 "+str(round(calcangle(i,j,k,MolSpec.CARTESIANS),3))+" 1 "+str(round(calcdihedral(i,j,k,l,MolSpec.CARTESIANS),3))+" 0  "+str(j+1)+" "+str(k+1)+" "+str(l+1))
					if freezeangle == 1:
						for connectk in MolSpec.CONNECTIVITY[k]: l = int(connectk.split("__")[0])-1
						fileout.write(" "+str(round(calcdist(j,i,MolSpec.CARTESIANS),3))+" 1 "+str(round(calcangle(i,j,k,MolSpec.CARTESIANS),3))+" 0 "+str(round(calcdihedral(i,j,k,l,MolSpec.CARTESIANS),3))+" 1  "+str(j+1)+" "+str(k+1)+" "+str(l+1))
				
				fileout.write("\n")
			fileout.write("\n")
			fileout.close()
		
		
		# Create the input file for optimization - if the file already exists don't overwrite it
		if JobSpec.PROGRAM == "Gaussian":
			if not os.path.exists(MolSpec.NAME+".com"): writeGaussian(self, JobSpec, MolSpec)
		if JobSpec.PROGRAM == "Mopac":
			if not os.path.exists(MolSpec.NAME+".mop"): writeMopac(self, JobSpec, MolSpec)


###############################################################
#                           MultMin                           #
#              Multiple Minimization of Conformers            #
###############################################################
def MultMin(CSEARCH, SEARCHPARAMS,CONFSPEC, MOLSPEC, JOB, start, log):
	
    i = 0
    while i < len(CSEARCH.NAME):
        for k in range(0,SEARCHPARAMS.POOL):
			# Replace geometry and energy of conformer with highly converged values
            if (i+k) < len(CSEARCH.NAME):
                CONFSPEC.NAME = CSEARCH.NAME[i+k]
                CONFSPEC.CARTESIANS = CSEARCH.CARTESIANS[i+k]
                CONFSPEC.CONNECTIVITY = CSEARCH.CONNECTIVITY[i+k]
                CONFSPEC.MULT = MOLSPEC.MULT
                CONFSPEC.MMTYPES = MOLSPEC.MMTYPES
                writeInput(JOB, CONFSPEC)
                submitJob(JOB, CONFSPEC, log)
        
        for k in range(0,SEARCHPARAMS.POOL):
            if (i+k) < len(CSEARCH.NAME):
                CONFSPEC.NAME = CSEARCH.NAME[i+k]
				
                # Make sure optimization is complete
                while(isJobFinished(JOB, CONFSPEC))==0: time.sleep(0.1)
				
                # Successful termination - extract details
                if isJobFinished(JOB, CONFSPEC) == 1:
                    #log.Write("   Extracting optimized structure from "+CONFSPEC.NAME+".out ...")
                    CONFSPEC =  getoutData(CONFSPEC)
                    CONFSPEC.CONNECTIVITY = CSEARCH.CONNECTIVITY[i+k]
                    CSEARCH.ENERGY[i+k] = CONFSPEC.ENERGY
                    CSEARCH.CARTESIANS[i+k] = CONFSPEC.CARTESIANS
                    CSEARCH.CPU[i+k] = CONFSPEC.CPU
                    CSEARCH.TORVAL[i+k] = getTorsion(CONFSPEC)
                    CSEARCH.ALLCPU.append(CONFSPEC.CPU)
                    if CONFSPEC.ENERGY < CSEARCH.GLOBMIN: CSEARCH.GLOBMIN = CONFSPEC.ENERGY
				
                samecheck = 0; toohigh = 0; isomerize = 0
                # Remove all non-essential files associated with this conformer
                CleanAfterJob(JOB, CONFSPEC, samecheck, toohigh, isomerize)
		
		
        i = i + SEARCHPARAMS.POOL
        if round(i*100.0/len(CSEARCH.NAME)) < 100.0 and round(i*100.0/len(CSEARCH.NAME)) > round((i-SEARCHPARAMS.POOL)*100.0/len(CSEARCH.NAME)):
            print "   "+str(round(i*100.0/len(CSEARCH.NAME)))+"% Complete ..."
		
	OrderConfs(CSEARCH, SEARCHPARAMS, start, log)
	
    print "\no  Checking for duplicate conformers ..."
    todel=[]
    for i in range((len(CSEARCH.NAME)-1),-1,-1):
        CONFSPEC.ENERGY = CSEARCH.ENERGY[i]
        CONFSPEC.NAME = CSEARCH.NAME[i]
        CONFSPEC.CARTESIANS = CSEARCH.CARTESIANS[i]
        CONFSPEC.CONNECTIVITY = CSEARCH.CONNECTIVITY[i]
        CONFSPEC.TIMESFOUND = CSEARCH.TIMESFOUND[i]
		
        print "checking",CONFSPEC.NAME
	#Check whether any stereogenic centres have been epimerized
	chircheck = checkchir(CONFSPEC, MOLSPEC, CSEARCH, SEARCHPARAMS)
	if chircheck[0] == 0: CONFSPEC.CONNECTIVITY = MOLSPEC.CONNECTIVITY; isomerize = 0
	else: isomerize = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected: "+chircheck[1]+" has been epimerized").ljust(50))

        #Check whether the molecule has isomerized - usually this is undesirable so filter out structural isomers
        concheck = checkconn(CONFSPEC, MOLSPEC, CSEARCH, SEARCHPARAMS)
        if concheck[0] == 0: CONFSPEC.CONNECTIVITY = MOLSPEC.CONNECTIVITY; isomerize = 0
        else: isomerize = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected: "+concheck[1]+concheck[2]+" has broken from "+concheck[3]+concheck[4]).ljust(50))
        
        if ((CONFSPEC.ENERGY - CSEARCH.GLOBMIN)*2625.5) < SEARCHPARAMS.DEMX: toohigh = 0
        else: toohigh = 1; log.Write("\n   "+(CONFSPEC.NAME+" is rejected due to high energy ... ").ljust(50))
		
        samecheck=0
        # Save or discard the optimized structure - reject if higher than global minimum by DEMX kJ/mol
        if toohigh == 0 and isomerize == 0:
            # Clustering - check if the energy and geometry match a previously saved structure
            for j in range(0, i):
                #print "checking ", CSEARCH.NAME[i], "against", CSEARCH.NAME[j]
                if (CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5 < -0.1: break
                if abs((CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5) < 0.5:
                    print "   "+CONFSPEC.NAME+" cf "+CSEARCH.NAME[j]+" = "+str((CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5)
                    #what if they are suspiciously similar in energy?
                    if checkSame(CONFSPEC, CSEARCH, SEARCHPARAMS, j) > 0 or checkSame(makemirror(CONFSPEC), CSEARCH, SEARCHPARAMS, j) > 0:
                        log.Write("   "+(CONFSPEC.NAME+" is a duplicate of conformer "+CSEARCH.NAME[j]).ljust(50))
                        CSEARCH.TIMESFOUND[j] = CSEARCH.TIMESFOUND[j] + CSEARCH.TIMESFOUND[i]
                        CSEARCH.NREJECT = CSEARCH.NREJECT + CSEARCH.TIMESFOUND[i]
                        todel.append(i)
                        CSEARCH.ENERGY[i] = 999999.9
                        break
		
        # Rejection - discard
        else:
            CSEARCH.NREJECT = CSEARCH.NREJECT + CSEARCH.TIMESFOUND[i]
            CSEARCH.ENERGY[i] = 999999.9
            log.Write("\n   "+(CSEARCH.NAME[i]+" is now discarded ... ").ljust(50))
            todel.append(i)
    #print todel
    if len(todel) !=0:
		OrderConfs(CSEARCH, SEARCHPARAMS, start, log)
		RemoveConformer(CSEARCH, todel)


