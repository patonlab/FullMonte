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
#                           FMC_1.py                          #
#              Monte Carlo Conformational Search              #
#         Dr Robert S Paton, University of Oxford 2010        #
###############################################################
#######  Written by:  Rob Paton ###############################
#######  Last modified:  Mar 20, 2013 #########################
###############################################################

# Python Libraries ############################################
import glob, subprocess, sys, os, random, math, tarfile
from numpy import *
###############################################################

# Full Monte Libaries #########################################
from FMTools import *
###############################################################

if __name__ == "__main__":
	
# An input file must be specified #############################
	instruct = "default"
	interactivemode = 1
	if len(sys.argv)>1: 
		for arg in sys.argv:
			if arg == "-setup": SETUPEXE(MOPAC_EXEC)
			if arg == "-background": interactivemode = 0
		filein = sys.argv[1].split(".")[0]
		if len(sys.argv[1].split("."))>1:
			if sys.argv[1].split(".")[1] == "com": filetype = "com"
			if sys.argv[1].split(".")[1] == "pdb": filetype = "pdb"
		if len(sys.argv)>2 and sys.argv[2] != "-background": instruct = sys.argv[2]
	else: print "\nWrong number of arguments used. Correct format: FullMonte struc.com [params] \n"; sys.exit()
###############################################################
	
# Check if MOPAC.EXEC is defined ##############################
	if not os.path.exists(MOPAC_EXEC): print "\no  Mopac executable cannot be found. Rerun Full Monte with -setup as one of the arguments\n"; sys.exit()
###############################################################

# Initialize the logfile for all text output ##################
	if os.path.exists(filein+"_fm.log") and interactivemode == 1: 
		var = raw_input("\no  Log file already exists! OK to overwrite this file ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": print "   Overwriting ..."
		else: print "\nExiting\n";  sys.exit(1)
	log = FMLog(filein,"log", "fm")
###############################################################	
	
# See if there are any remaining files from a previous run ####
	if os.path.exists(filein+"_fm.tgz") and interactivemode == 1: 
		var = raw_input("\no  Tarfile already exists! OK to overwrite these structures ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": os.remove(filein+"_fm.tgz")
		else: print "\nExiting\n";  sys.exit(1)
        if len(glob.glob(filein+"*step*"))>0:
		if interactivemode == 1:
                	print glob.glob(filein+"*step*")
			var = raw_input("\no  Some optimization ouput files already exist! OK to overwrite these structures ? (Y/N) ")
                	if var.lower() == "y" or var.lower() == "": 
				for file in glob.glob(filein+"*step*"): os.remove(file)
                	else: print "\nExiting\n";  sys.exit(1)
		else:
			for file in glob.glob(filein+"*step*"): os.remove(file)
###############################################################			

# Open the structure file #####################################
	log.Write("\no  Extracting molecule from "+filein+"."+filetype+" ...")
	MOLSPEC = getinData(filein,log)
###############################################################	
	
# Open the specified parameter file for Monte Carlo parameters 
# (default values will be used if not supplied) ###############
	if instruct!="default": log.Write("\no  Extracting conformational search parameters from "+instruct+" ...")
	else: log.Write("\no  No FullMonte parameters specified! Using default values ...")
	SEARCHPARAMS = getParams(MOLSPEC, instruct,log)
###############################################################	
		
# Model Chemistry to be used ##################################
	for level in ["AM1", "PM3", "PM6", "PM7", "PM6-DH2"]:
		if SEARCHPARAMS.LEVL.upper() == level: JOB = JobSpec("Mopac")
	for level in ["UFF"]:
		if SEARCHPARAMS.LEVL.upper() == level: JOB = JobSpec("Gaussian")
	if JOB.PROGRAM != "Mopac" and JOB.PROGRAM != "Gaussian": log.Fatal("\no  "+SEARCHPARAMS.LEVL+" Level of Theory Not Yet Supported ... ")
	JOB.JOBTYPE = SEARCHPARAMS.LEVL
	log.Write("\no  Using "+JOB.JOBTYPE+" level of theory ... ")
###############################################################	
	
# For PM6, a dispersion/H-bond correction is available ########
	if JOB.PROGRAM == "Mopac" and JOB.JOBTYPE.upper() == "PM6" and interactivemode == 1:
		var = raw_input("\no  Use dispersion and H-bonding correction ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": JOB.JOBTYPE = "PM6-DH2"; log.Write("   Dispersion and H-bond correction on ... ") 
###############################################################	
	
# Solvation with CPCM #########################################
	if JOB.PROGRAM == "Gaussian" and interactivemode == 1:
		var = raw_input("\no  Use CPCM solvation ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "":
			EPS = raw_input("\n   Enter solvent name (default=diethylether) ")
			if EPS == "": EPS = "(cpcm,solvent=diethylether)"
			JOB.JOBTYPE = JOB.JOBTYPE+" scrf"+EPS; log.Write("   CPCM solvation correction on ... ")
###############################################################

# Solvation with COSMO ########################################
	if JOB.PROGRAM == "Mopac" and interactivemode == 1:
		var = raw_input("\no  Use COSMO solvation ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": 
			EPS = raw_input("\n   Enter solvent dielectric constant (default = 78.4) ")
			if EPS == "": EPS = "78.4"
			JOB.JOBTYPE = JOB.JOBTYPE+" EPS="+EPS; log.Write("   COSMO solvation correction on ... ") 
###############################################################

# MM correction for N planarity? ##############################
	if JOB.PROGRAM == "Mopac":
		for atom in MOLSPEC.ATOMTYPES:
			if atom == "N": 
				if interactivemode == 1:
					var = raw_input("\no  Use molecular mechanics correction for planar nitrogen atoms ? (Y/N) ")
					if var.lower() == "y" or var.lower() == "": JOB.JOBTYPE = JOB.JOBTYPE+" mmok"; log.Write("   MM correction on ... ") 
					else: log.Write("   No MM correction ... ") 
					break
				else:
					JOB.JOBTYPE = JOB.JOBTYPE+" mmok"; log.Write("   MM correction on ... ")
###############################################################				
	
# Check for any constraints specified #########################
	if hasattr(MOLSPEC, "CONSTRAINED"):
		JOB.CONSTRAINED = MOLSPEC.CONSTRAINED
		#print MOLSPEC.CONSTRAINED
		for const in MOLSPEC.CONSTRAINED: 
			if len(const) == 1: log.Write("\no  The Cartesian position of "+str(const[0]+1)+" will be constrained ...")
			if len(const) == 2: log.Write("\no  The distance "+str(const[0]+1)+"-"+str(const[1]+1)+" will be constrained ...")
			if len(const) == 3: log.Write("\no  The angle "+str(const[0]+1)+"-"+str(const[1]+1)+"-"+str(const[2]+1)+" will be constrained ...")
			if len(const) == 4: log.Write("\no  The dihedral "+str(const[0]+1)+"-"+str(const[1]+1)+"-"+str(const[2]+1)+"-"+str(const[3]+1)+" will be constrained ...")
		if JOB.PROGRAM == "Gaussian":
			if len(MOLSPEC.CONSTRAINED)!=0: JOB.JOBTYPE = "opt(small,modredundant,loose) "+JOB.JOBTYPE
			else: JOB.JOBTYPE = "opt(loose) "+JOB.JOBTYPE
###############################################################

if JOB.PROGRAM == "Gaussian": 
	JOB.JOBTYPE = "# geom=connectivity "+JOB.JOBTYPE
	JOB.NPROC = SEARCHPARAMS.POOL

# Monte Carlo or Systematic (for comparison) ##################
if interactivemode == 1:
	var = raw_input("\no  MCMM (Y) or SUMM (N) ? (Y/N) ")
	if var.lower() == "y" or var.lower() == "": SEARCHPARAMS.CSEARCH = "MCMM"
	if var.lower() == "n": SEARCHPARAMS.CSEARCH = "SUMM" 
###############################################################

# Perform an optimization of the starting geometry ############
MOLSPEC.NAME = MOLSPEC.NAME+"_step_0"
writeInput(JOB, MOLSPEC)
submitJob(JOB, MOLSPEC,log)
while isJobFinished(JOB, MOLSPEC) == 0: time.sleep(0.1)
if isJobFinished(JOB, MOLSPEC) == -1: log.Fatal("\nFATAL ERROR: Optimization of [ %s ] stuck"%file)
if isJobFinished(JOB, MOLSPEC) == 2: log.Fatal("\nFATAL ERROR: Optimization of [ %s ] failed"%file)
###############################################################	
		
# Read the output from the optimization then clean up #########
MOLSPEC.CARTESIANS =  getoutData(MOLSPEC).CARTESIANS
MOLSPEC.ENERGY =  getoutData(MOLSPEC).ENERGY
for suffix in [".com", ".csh", ".mop", ".arc", ".aux", ".joblog", ".errlog", ".chk"]:
	if os.path.exists(MOLSPEC.NAME+suffix): os.remove(MOLSPEC.NAME+suffix)
###############################################################
	
# Assign variable torsions, number of separate molecules and ##
# (eventually when I get round to it) rings ###################
# If number of steps is not assigned use 3^rotatable torsions #
FMVAR = Assign_Variables(MOLSPEC, SEARCHPARAMS, log)
if SEARCHPARAMS.CSEARCH == "MCMM" and SEARCHPARAMS.STEP == 0:
	SEARCHPARAMS.STEP = int(math.pow(3,FMVAR.MCNV))
	SEARCHPARAMS.STEP = SEARCHPARAMS.STEP + int(math.pow(3,len(FMVAR.RING)-1))
	#SEARCHPARAMS.STEP = SEARCHPARAMS.STEP + 20
if SEARCHPARAMS.CSEARCH == "SUMM":
	if interactivemode == 1:
		SEARCHPARAMS.ITVL = raw_input("\no  Required Interval (degrees) for systematic rotations ? ")
		if SEARCHPARAMS.ITVL == "": SEARCHPARAMS.ITVL = "60"
		SEARCHPARAMS.ITVL = int(SEARCHPARAMS.ITVL)
	interval = SEARCHPARAMS.ITVL; ninterval = 360.0/interval
	SEARCHPARAMS.STEP = int(math.pow(ninterval,FMVAR.MCNV))-1
###############################################################


# MONTE CARLO SEARCH ##########################################
start = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
asciiArt(start); Writeintro(MOLSPEC, SEARCHPARAMS, FMVAR, start, log) 

CONFSPEC = MOLSPEC
CSEARCH.NAME.append(MOLSPEC.NAME)
CSEARCH.CARTESIANS.append(MOLSPEC.CARTESIANS)
CSEARCH.TORVAL = [getTorsion(MOLSPEC)]
CSEARCH.CONNECTIVITY.append(MOLSPEC.CONNECTIVITY)
CSEARCH.ENERGY = [MOLSPEC.ENERGY]
CSEARCH.GLOBMIN = MOLSPEC.ENERGY
CSEARCH.CPU = [getoutData(MOLSPEC).CPU]
CSEARCH.ALLCPU = [getoutData(MOLSPEC).CPU]
CSEARCH.LASTFOUND = 0
CSEARCH.CLASH = [0]

# Stop once number of steps exceeded or no new conformers found	
while CSEARCH.STEP*SEARCHPARAMS.POOL <= SEARCHPARAMS.STEP:
	log.Write("o  STEP "+str(CSEARCH.STEP)+": Generating "+str(SEARCHPARAMS.POOL)+" structures ...")


# Setting the geometry that will be altered to generate new conformers - only relevant to MCMM
	if SEARCHPARAMS.CSEARCH == "MCMM":
		for i in range(0, CSEARCH.NSAVED):
			if CSEARCH.ENERGY[i] - CSEARCH.GLOBMIN == 0.0: startgeom = i
		
# Generate new geometries 
	for i in range(((CSEARCH.STEP-1)*SEARCHPARAMS.POOL+1),((CSEARCH.STEP)*SEARCHPARAMS.POOL)+1):
		CONFSPEC.NAME = filein+"_step_"+str(i)
		if SEARCHPARAMS.CSEARCH == "MCMM":
			if SEARCHPARAMS.MCSS == "Uniform Usage Directed":
				for j in range(0, CSEARCH.NSAVED):
					if (CSEARCH.ENERGY[j] - CSEARCH.GLOBMIN) * 2625.5 < SEARCHPARAMS.EWIN:
						if CSEARCH.USED[j] < CSEARCH.USED[startgeom]: startgeom = j
						if CSEARCH.USED[j] == CSEARCH.USED[startgeom] and CSEARCH.ENERGY[j] < CSEARCH.ENERGY[startgeom]: startgeom = j
				CSEARCH.USED[startgeom] = CSEARCH.USED[startgeom] + 1

		NBcontacts = 1

		if SEARCHPARAMS.CSEARCH == "SUMM":
			torsiontwist = [0]* len(FMVAR.TORSION)
			for j in range(0,len(FMVAR.TORSION)-1):
				#print i, math.pow(ninterval, (len(torsiontwist)-j-1))
				#print int(i)/int(math.pow(ninterval, (len(torsiontwist)-j-1))), int(i)%int(math.pow(ninterval, (len(torsiontwist)-j-1)))
				torsiontwist[j] = int(i)/int(math.pow(ninterval, (len(torsiontwist)-j-1))) * interval
				while torsiontwist[j] >= 360.0: torsiontwist[j] = torsiontwist[j] - 360.0
			torsiontwist[len(FMVAR.TORSION)-1] = int(i)%int(math.pow(ninterval, (len(torsiontwist)-j-1))) * interval
			print "   Dihedral twists in degrees:", torsiontwist

			if FMVAR.MCNV != 0: FMVAR.ADJUST = []
			for j in range(0,len(FMVAR.TORSION)): FMVAR.ADJUST.append([int(FMVAR.TORSION[j][0])+1, int(FMVAR.TORSION[j][1])+1, int(torsiontwist[j])])
			if hasattr(FMVAR, "ADJUST"):
				#print FMVAR.ADJUST
				for torsion in FMVAR.ADJUST:
					if torsion[2] != 0: CONFSPEC.CARTESIANS = AtomRot(MOLSPEC, torsion, MOLSPEC.CARTESIANS)
			NBcontacts = checkDists(CONFSPEC, SEARCHPARAMS)
			CSEARCH.CLASH.append(NBcontacts)
		
		if SEARCHPARAMS.CSEARCH == "MCMM":
			
			attempts = 0
			while NBcontacts > 0 and attempts < 100:
				CONFSPEC.CARTESIANS = []
				# The coordinates of the lowest energy, least used structure will be altered
				print "   STARTING FROM GEOMETRY OF", CSEARCH.NAME[startgeom]
				#print CSEARCH.CARTESIANS[startgeom]
				#print "SUM1", sum(CSEARCH.CARTESIANS[startgeom])
				for i in range (0,len(CSEARCH.CARTESIANS[startgeom])):
					CONFSPEC.CARTESIANS.append([])
					#print (CSEARCH.CARTESIANS[startgeom][i])
					for cart in (CSEARCH.CARTESIANS[startgeom][i]):
						#print i, cart
						CONFSPEC.CARTESIANS[i].append(cart)
				CONFSPEC.CONNECTIVITY = CSEARCH.CONNECTIVITY[startgeom]
				CONFSPEC.ATOMTYPES = MOLSPEC.ATOMTYPES
				CONFSPEC.CHARGE = MOLSPEC.CHARGE
				CONFSPEC.MULT = MOLSPEC.MULT
				CONFSPEC.MMTYPES = MOLSPEC.MMTYPES
				nrandom = random.randint(FMVAR.MCNVmin, FMVAR.MCNVmax)
				#print CONFSPEC.CARTESIANS
				#print getTorsion(CONFSPEC)
	
				if FMVAR.MCNV != 0:
					FMVAR.ADJUST = []
					for dihedral in random.sample(FMVAR.TORSION, nrandom): 
						FMVAR.ADJUST.append([int(dihedral[0])+1, int(dihedral[1])+1, random.randint(0,360)])
					
					if len(FMVAR.ETOZ) > 0:
						ezisomerize = random.choice([0,1])
						for dihedral in random.sample(FMVAR.ETOZ,ezisomerize): 
							#print "ETOZ",dihedral,ezisomerize
							FMVAR.ADJUST.append([int(dihedral[0]), int(dihedral[1]), 180])
				
				
# Take input geometry and apply specified torsional changes
				if hasattr(FMVAR, "ADJUST"):
					#print FMVAR.ADJUST
					for torsion in FMVAR.ADJUST: CONFSPEC.CARTESIANS = AtomRot(MOLSPEC, torsion, CONFSPEC.CARTESIANS)
				#print "AFTER ADJUSTMENT:"
				#print getTorsion(CONFSPEC) 
				
# For separate molecules, alter the distances and orientations between a random number of them
				if FMVAR.NMOLS > 1:
					CONFSPEC.CARTESIANS = translateMol(FMVAR, CONFSPEC)
					CONFSPEC.CARTESIANS = rotateMol(FMVAR, CONFSPEC)
				
				if FMVAR.MCRI > 0:
					print "   Detected a ring substructure: atoms", FMVAR.RING, "are connected"
					rcom = find_centroid(FMVAR.RING,CONFSPEC)[3:]
					
					coeffplane, xav, yav, zav, rotated = find_coeffplane(FMVAR.RING,CONFSPEC)
					xcoeff= coeffplane.tolist()[0][0]; ycoeff= coeffplane.tolist()[1][0]; cval= coeffplane.tolist()[2][0]
					
					#print "Equation of best-fit plane:","z="+str(xcoeff)+"x+"+str(ycoeff)+"y+"+str(cval)	#This gives the equation for the plane of best-fit
					####################Make unit vector
					rawvector=array([xcoeff,ycoeff,-1]) #Need to make into unit vector
					x=float(rawvector[0]); y=float(rawvector[1]); z=float(rawvector[2])
					normfactor=1/(x**2+y**2+z**2)**0.5
					x=x*normfactor; y=y*normfactor; z=z*normfactor
					if z<0: z=-z;y=-y;x=-x #Sign flip if z is negative
					#print "   Unit vector:", x, y, z #The length of this vector is 1
					
					if rotated == 1:
						print "************ coordinated system was rotated! ***********"
						old_x = z; old_y = x; old_z = y
						if old_z<0: old_z=-old_z;old_y=-old_y;old_x=-old_x
						print "Unit vector:", old_x, old_y, old_z
						x = old_x; y = old_y; z = old_z
					if rotated == 2:
						print "************ coordinated system was rotated! ***********"
						old_x = y; old_y = z; old_z = x
						if old_z<0: old_z=-old_z;old_y=-old_y;old_x=-old_x
						print "Unit vector:", old_x, old_y, old_z
						x = old_x; y = old_y; z = old_z
					if rotated == 3:
						print "didn't I tell you this was a bad idea?"				
				
					nrandom = random.randint(1, len(FMVAR.RING)/2)
					
					for atomid in random.sample(FMVAR.RING, nrandom):
					
						count = 0; stop=0; currentatom=[]; nextlot=[]
						currentatom.append([atomid])
						#print currentatom
						while count<100 and stop==0:
							nextlot=[]; ringneighbour = []
							for onecurrentatom in currentatom[count]:
								#print "onecurrentatom", (onecurrentatom+1)
								for partners in CONFSPEC.CONNECTIVITY[onecurrentatom]:
									#print "partners", partners
									inf = partners.split("__")
									for n in range(0,len(inf)/2):
										noback=0
										for onepreviousatom in currentatom[count-1]:
											if (int(inf[2*n])-1)==onepreviousatom: noback=noback+1
										for onepreviousatom in currentatom[count]:
											if (int(inf[2*n])-1)==onepreviousatom: noback=noback+1
										for ringatom in FMVAR.RING:
											if (int(inf[2*n])-1)==ringatom:
												#ringneighbour.append(int(inf[2*n])-1)
												noback=noback+1
										if noback==0: nextlot.append(int(inf[2*n])-1)
							count=count+1
							#print count
							if len(nextlot) == 0:stop=stop+1
							currentatom.append(nextlot)
						
						for partners in CONFSPEC.CONNECTIVITY[atomid]:
							inf = partners.split("__")
							for ringatom in FMVAR.RING:
								if (int(inf[2*n])-1)==ringatom:
									ringneighbour.append(int(inf[2*n])-1)
						#print "SUM2", sum(CSEARCH.CARTESIANS[startgeom])
						print "   Neighbours in the ring", ringneighbour
						oldvecA = [CONFSPEC.CARTESIANS[atomid][0] - CONFSPEC.CARTESIANS[ringneighbour[0]][0], CONFSPEC.CARTESIANS[atomid][1] - CONFSPEC.CARTESIANS[ringneighbour[0]][1], CONFSPEC.CARTESIANS[atomid][2] - CONFSPEC.CARTESIANS[ringneighbour[0]][2]]
						oldvecB = [CONFSPEC.CARTESIANS[atomid][0] - CONFSPEC.CARTESIANS[ringneighbour[1]][0], CONFSPEC.CARTESIANS[atomid][1] - CONFSPEC.CARTESIANS[ringneighbour[1]][1], CONFSPEC.CARTESIANS[atomid][2] - CONFSPEC.CARTESIANS[ringneighbour[1]][2]]
						oldnorm = [oldvecA[1]*oldvecB[2]-oldvecA[2]*oldvecB[1], oldvecA[2]*oldvecB[0]-oldvecA[0]*oldvecB[2], oldvecA[0]*oldvecB[1]-oldvecA[1]*oldvecB[0]]
						oldmag = (oldnorm[0]**2 + oldnorm[1]**2 + oldnorm[2]**2) ** 0.5
						oldnorm = [oldnorm[0]/oldmag, oldnorm[1]/oldmag, oldnorm[2]/oldmag]
						
						#print "SUM3", sum(CSEARCH.CARTESIANS[startgeom])
						mag = 1.0; magnitude = random.uniform(-1*mag,mag)
						print CONFSPEC.CARTESIANS[atomid]
						CONFSPEC.CARTESIANS[atomid][0] = CONFSPEC.CARTESIANS[atomid][0] + x * magnitude
						CONFSPEC.CARTESIANS[atomid][1] = CONFSPEC.CARTESIANS[atomid][1] + y * magnitude
						CONFSPEC.CARTESIANS[atomid][2] = CONFSPEC.CARTESIANS[atomid][2] + z * magnitude
						print "      TRANSLATING ATOM", (atomid+1), "BY ", [magnitude*x,magnitude*y,magnitude*z]
						print CONFSPEC.CARTESIANS[atomid]
						
						newvecA = [CONFSPEC.CARTESIANS[atomid][0] - CONFSPEC.CARTESIANS[ringneighbour[0]][0], CONFSPEC.CARTESIANS[atomid][1] - CONFSPEC.CARTESIANS[ringneighbour[0]][1], CONFSPEC.CARTESIANS[atomid][2] - CONFSPEC.CARTESIANS[ringneighbour[0]][2]]
						newvecB = [CONFSPEC.CARTESIANS[atomid][0] - CONFSPEC.CARTESIANS[ringneighbour[1]][0], CONFSPEC.CARTESIANS[atomid][1] - CONFSPEC.CARTESIANS[ringneighbour[1]][1], CONFSPEC.CARTESIANS[atomid][2] - CONFSPEC.CARTESIANS[ringneighbour[1]][2]]
						newnorm = [newvecA[1]*newvecB[2]-newvecA[2]*newvecB[1], newvecA[2]*newvecB[0]-newvecA[0]*newvecB[2], newvecA[0]*newvecB[1]-newvecA[1]*newvecB[0]]
						newmag = (newnorm[0]**2 + newnorm[1]**2 + newnorm[2]**2) ** 0.5
						newnorm = [newnorm[0]/newmag, newnorm[1]/newmag, newnorm[2]/newmag]

						#print "SUM4", sum(CSEARCH.CARTESIANS[startgeom])
						#print oldnorm, newnorm
						rotang = 180.0/math.pi*math.acos(oldnorm[0]*newnorm[0]+oldnorm[1]*newnorm[1]+oldnorm[2]*newnorm[2])
						print "     A ROTATION THROUGH", rotang
						#newvec = [CONFSPEC.CARTESIANS[atomid][0] - rcom[0], CONFSPEC.CARTESIANS[atomid][1] - rcom[1], CONFSPEC.CARTESIANS[atomid][2] - rcom[2]]
						#coords = [oldvec, [0.0,0.0,0.0],newvec]
						#rotang = calcangle(0,1,2,coords)

						#vectorAB =  [oldvec[1]*z-oldvec[2]*y, oldvec[2]*x -oldvec[0]*y, oldvec[0]*y-oldvec[1]*x]
						#print "      OLD VECTOR FROM COM", oldvec
						#print "      NEW VECTOR FROM COM", newvec
						#print "      A ROTATION OF",rotang
						#print "      ABOUT", vectorAB
						#distAB = math.sqrt(vectorAB[0]*vectorAB[0]+vectorAB[1]*vectorAB[1]+vectorAB[2]*vectorAB[2])
						#unitAB = [vectorAB[0]/distAB,vectorAB[1]/distAB,vectorAB[2]/distAB]


						for subst in currentatom[1:]:
							for nextatomid in subst:
								print "         ALSO TRANSLATING ATOM", (nextatomid+1), "BY ", [magnitude*x,magnitude*y,magnitude*z]
								CONFSPEC.CARTESIANS[nextatomid][0] = CONFSPEC.CARTESIANS[nextatomid][0] + x * magnitude
								CONFSPEC.CARTESIANS[nextatomid][1] = CONFSPEC.CARTESIANS[nextatomid][1] + y * magnitude
								CONFSPEC.CARTESIANS[nextatomid][2] = CONFSPEC.CARTESIANS[nextatomid][2] + z * magnitude
								
								print "         ALSO ROTATING ATOM", (nextatomid+1), "THROUGH ANGLE", rotang, "ABOUT THE AXIS", ringneighbour
								print CONFSPEC.CARTESIANS[nextatomid]
								
								atomA = ringneighbour[0]; atomB = ringneighbour[1]
								vectorAB = [CONFSPEC.CARTESIANS[atomB][0]-CONFSPEC.CARTESIANS[atomA][0],CONFSPEC.CARTESIANS[atomB][1]-CONFSPEC.CARTESIANS[atomA][1],CONFSPEC.CARTESIANS[atomB][2]-CONFSPEC.CARTESIANS[atomA][2]]
								distAB = (vectorAB[0]**2+vectorAB[1]**2+vectorAB[2]**2) ** 0.5
								unitAB = [vectorAB[0]/distAB,vectorAB[1]/distAB,vectorAB[2]/distAB]

								dotproduct = unitAB[0]*(CONFSPEC.CARTESIANS[nextatomid][0] - CONFSPEC.CARTESIANS[atomid][0]) + unitAB[1]*(CONFSPEC.CARTESIANS[nextatomid][1] - CONFSPEC.CARTESIANS[atomid][1]) + unitAB[2]*(CONFSPEC.CARTESIANS[nextatomid][2] - CONFSPEC.CARTESIANS[atomid][2])
								centre = [CONFSPEC.CARTESIANS[atomid][0] + dotproduct*unitAB[0], CONFSPEC.CARTESIANS[atomid][1] + dotproduct*unitAB[1], CONFSPEC.CARTESIANS[atomid][2] + dotproduct*unitAB[2]]
								v = [CONFSPEC.CARTESIANS[nextatomid][0] - centre[0], CONFSPEC.CARTESIANS[nextatomid][1] - centre[1], CONFSPEC.CARTESIANS[nextatomid][2] - centre[2]]
								theta = -float(rotang)/180.0*math.pi
								d = (v[0]**2+v[1]**2+v[2]**2) ** 0.5
								px = v[0]*math.cos(theta) + v[1]*math.sin(theta)*unitAB[2] - v[2]*math.sin(theta)*unitAB[1]
								py = v[1]*math.cos(theta) + v[2]*math.sin(theta)*unitAB[0] - v[0]*math.sin(theta)*unitAB[2]
								pz = v[2]*math.cos(theta) + v[0]*math.sin(theta)*unitAB[1] - v[1]*math.sin(theta)*unitAB[0]
								CONFSPEC.CARTESIANS[nextatomid] = [px + centre[0], py + centre[1], pz + centre[2]]
								print CONFSPEC.CARTESIANS[nextatomid]
				
# Check for any VDW contacts smaller than specified limits
				NBcontacts = checkDists(CONFSPEC, SEARCHPARAMS)
				if NBcontacts != 0:
					attempts = attempts + 1
					print "   EXTREME CONTACTS: TRYING AGAIN..."
			CSEARCH.CLASH.append(NBcontacts)

# Write input and optimize geometry
		#print CONFSPEC.CARTESIANS
		if NBcontacts == 0 and attempts < 100: print "   SUCCESSFULLY ALTERED GEOM..."
		if attempts > 99 and NBcontacts != 0: print "   EXCEEDED MAXIMUM ATTEMPTS TO ALTER GEOM..."
		writeInput(JOB, CONFSPEC); submitJob(JOB, CONFSPEC, log)

# Filter after optimization
	for i in range(((CSEARCH.STEP-1)*SEARCHPARAMS.POOL+1),((CSEARCH.STEP)*SEARCHPARAMS.POOL)+1):
		CONFSPEC.NAME = filein+"_step_"+str(i)
		
# Make sure optimization is complete
		#print i, CSEARCH.CLASH
		if CSEARCH.CLASH[i] == 0:
			while isJobFinished(JOB, CONFSPEC) == 0: time.sleep(0.1)
			
# Successful termination - extract details
			if isJobFinished(JOB, CONFSPEC) == 1:
				CONFSPEC =  getoutData(CONFSPEC)
				CONFSPEC.ATOMTYPES = MOLSPEC.ATOMTYPES
				CSEARCH.ALLCPU.append(CONFSPEC.CPU)
				
#Check whether the molecule has isomerized - usually this is undesirable so filter out structural isomers
				concheck = checkconn(CONFSPEC, MOLSPEC, CSEARCH, SEARCHPARAMS)
				if concheck[0] == 0: CONFSPEC.CONNECTIVITY = MOLSPEC.CONNECTIVITY; isomerize = 0
				else: isomerize = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected: "+concheck[1]+concheck[2]+" has broken from "+concheck[3]+concheck[4]).ljust(50))
			
#Check whether any stereogenic centres have been epimerized
				chircheck = checkchir(CONFSPEC, MOLSPEC, CSEARCH, SEARCHPARAMS)
				if chircheck[0] == 0: CONFSPEC.CONNECTIVITY = MOLSPEC.CONNECTIVITY; isomerize = 0
                                else: isomerize = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected: atom "+str(chircheck[1])+" has been epimerized").ljust(50))

	
#Check whether the molecule has high energy				
				if ((CONFSPEC.ENERGY-CSEARCH.GLOBMIN)*2625.5) < SEARCHPARAMS.DEMX: toohigh = 0
				else: toohigh = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected due to high energy ... ").ljust(50))
				samecheck=0

# Save or discard the optimized structure - reject if higher than global minimum by DEMX kJ/mol
				if toohigh == 0 and isomerize == 0: 
					for j in range(0, CSEARCH.NSAVED): 
						#if (CONFSPEC.ENERGY-CSEARCH.GLOBMIN)*2625.5 < -0.1: break
						if abs((CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5) < ECOMP:
							print "   COMPARING   "+CONFSPEC.NAME+"   "+str(CONFSPEC.ENERGY)+" cf "+CSEARCH.NAME[j]+"   "+str(CSEARCH.ENERGY[j])+": ediff = "+str((CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5)
							if checkSame(CONFSPEC, CSEARCH, SEARCHPARAMS, j) > 0 or checkSame(makemirror(CONFSPEC), CSEARCH, SEARCHPARAMS, j) > 0:
								log.Write("   "+(CONFSPEC.NAME+" is a duplicate of conformer "+CSEARCH.NAME[j]+" ... ").ljust(50))
								CSEARCH.TIMESFOUND[j] = CSEARCH.TIMESFOUND[j] + 1
								CSEARCH.NREJECT = CSEARCH.NREJECT + 1
								samecheck = samecheck + 1
								break
					
# Unique conformation with low energy! ########################			
					if samecheck == 0:
						if CONFSPEC.ENERGY < CSEARCH.GLOBMIN:
							CSEARCH.GLOBMIN = CONFSPEC.ENERGY
							log.Write("   "+(CONFSPEC.NAME+" is a new Global Minimum!").ljust(80)+("E = "+str(CSEARCH.GLOBMIN)).rjust(rightcol))
						else : log.Write("   "+(CONFSPEC.NAME+" is saved").ljust(80)+("E = "+str(CONFSPEC.ENERGY)).rjust(rightcol))
						AddConformer(CSEARCH, CONFSPEC)
						if (CONFSPEC.ENERGY-CSEARCH.GLOBMIN)*2625.5 < SEARCHPARAMS.EWIN: CSEARCH.LASTFOUND = CSEARCH.STEP*SEARCHPARAMS.POOL
###############################################################
			
# Rejection - discard #########################################
				else: CSEARCH.NREJECT = CSEARCH.NREJECT + 1
			else: 
				log.Write("\n   Unsuccessful optimization of "+CONFSPEC.NAME+".out ...")
				CSEARCH.NFAILED = CSEARCH.NFAILED + 1

			CleanAfterJob(JOB, CONFSPEC, samecheck, toohigh, isomerize)
			OrderConfs(CSEARCH, SEARCHPARAMS, start, log)
###############################################################

		else: log.Write("\n  "+CONFSPEC.NAME+" not run due to extreme steric crowding ...")
			
# End of step - update the search statistics ##################
	if (CSEARCH.STEP*SEARCHPARAMS.POOL-CSEARCH.NFAILED) != 0: CSEARCH.AERATE = float(CSEARCH.STEP*SEARCHPARAMS.POOL-CSEARCH.NREJECT-CSEARCH.NFAILED)/float(CSEARCH.STEP*SEARCHPARAMS.POOL-CSEARCH.NFAILED)*100.0
	else: CSEARCH.AERATE = 0.0
	if len(CSEARCH.TIMESFOUND) > 0:
		for dup in CSEARCH.TIMESFOUND:
			if dup < CSEARCH.DMIN: CSEARCH.DMIN = dup
	else: CSEARCH.DMIN = 0
###############################################################

#Tidy up the geometries and output the results ################
	if CSEARCH.STEP % 100 == 0:
		todel=[]
		#print (CSEARCH.NAME)
		for i in range(0,len(CSEARCH.NAME)):
			if ((CSEARCH.ENERGY[i] - CSEARCH.GLOBMIN)*2625.5) > SEARCHPARAMS.DEMX or (i > 199):
				#print CSEARCH.NAME[i], "earmarked for deletion"
				todel.append(i)
		if len(todel) !=0: RemoveConformer(CSEARCH, todel)
		
		WriteSummary(CSEARCH, SEARCHPARAMS, start, log)	
		CleanUp(CSEARCH, SEARCHPARAMS, filein, log)
		makeGVformat(filein, MOLSPEC, CSEARCH, SEARCHPARAMS, "fm"); makePDBformat(filein, MOLSPEC, CSEARCH, "fm")	
		
	#End of step - update step number
	CSEARCH.STEP = CSEARCH.STEP + 1
###############################################################
		
#Summary of completed Full Monte search #######################
CSEARCH.COMPLETE = 1
if SEARCHPARAMS.CSEARCH == "MCMM":
		if (CSEARCH.STEP*SEARCHPARAMS.POOL-CSEARCH.LASTFOUND) >100: log.Write("\no  Full Monte stopped finding new conformers ...")
		CSEARCH.STEP = SEARCHPARAMS.STEP
else:log.Write("\no  Full Monte completed all "+str(SEARCHPARAMS.STEP)+" steps...")
WriteSummary(CSEARCH, SEARCHPARAMS, start, log)	
CleanUp(CSEARCH, SEARCHPARAMS, filein, log)
makeGVformat(filein, MOLSPEC, CSEARCH, SEARCHPARAMS, "fm"); makePDBformat(filein, MOLSPEC, CSEARCH, "fm")
###############################################################

# Multiple Minimization with higher convergence criterion ####
if SEARCHPARAMS.MMIN >0:
	log.Write("\no  Reoptimizing conformers with strict convergence criteria ...")
	if JOB.PROGRAM == "Mopac": JOB.JOBTYPE = JOB.JOBTYPE+" gnorm=0.0 "
	if JOB.PROGRAM == "Gaussian": JOB.JOBTYPE = JOB.JOBTYPE.replace("loose", "")
	MultMin(CSEARCH, SEARCHPARAMS,CONFSPEC, MOLSPEC, JOB, start, log)
	#OrderConfs(CSEARCH, SEARCHPARAMS, start, log)	
	CSEARCH.GLOBMIN = CSEARCH.ENERGY[0]
###############################################################

# Final Summary of Full Monte search ##########################
WriteSummary(CSEARCH, SEARCHPARAMS, start, log)	
if os.path.isfile(filein+"_fm.tgz") == 1: os.remove(filein+"_fm.tgz")
CleanUp(CSEARCH, SEARCHPARAMS, filein, log)
makeGVformat(filein, MOLSPEC, CSEARCH, SEARCHPARAMS, "fm"); makePDBformat(filein, MOLSPEC, CSEARCH, "fm")
end = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
asciiArt(end); log.Write(normaltermination); log.Finalize()	
###############################################################					
