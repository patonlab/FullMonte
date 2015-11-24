#!/usr/bin/python

                  ###                     ###          ###      
                  ###                     ###          ###      
                  ###                     ###          ###      
#####b.   ####b.  ###### .d##b.  #####b.  ###  ####b.  #####b.  
### "##b     "##b ###   d##""##b ### "##b ###     "##b ### "##b 
###  ### .d###### ###   ###  ### ###  ### ### .d###### ###  ### 
### d##P ###  ### Y##b. Y##..##P ###  ### ### ###  ### ### d##P 
#####P"  "Y######  "Y### "Y##P"  ###  ### ### "Y###### #####P"  
###                                                             
###                                                             
###

###############################################################
#                           FMC_2.py                          #
#             QM opts on Monte Carlo Conformers               #
#         Dr Robert S Paton, University of Oxford 2010        #
###############################################################
 

# Python Libraries ############################################
import subprocess, sys, os, random, math, tarfile
###############################################################

# Full Monte Libaries #########################################
from FMTools import *
###############################################################

if __name__ == "__main__":

    # A FullMonte output file is required
    if len(sys.argv)>=4:
        filein = sys.argv[1].split("_fm.")[0]
        comfile = sys.argv[2].split(".")[0]
        # Energy range for submitting QM jobs (kJ/mol)
        CSEARCH.ENERGYRANGE = float(sys.argv[3])
	
    else: print "\nWrong number of arguments used. Correct format: FullMonte struc_fm.out struc.com energyrange\n"; sys.exit()
#    if not os.path.exists(filein+"_fm.out") or not os.path.exists(filein+"_fm.sdf"): print "\nFullMonte output not found! Exiting ...\n";  sys.exit(1)
	
    # Initialize the logfile for all text output
    print "\no  Creating Full Monte QM log file for",filein+"_qm ..."
    log = FMLog(filein, "dat", "qm")
    qmtgz = filein+"_qm.tgz"
    qmsuffix = "dft"
    if len(sys.argv)==5: qmsuffix = sys.argv[4]
	
    # Open the specified Gaussian structure file
    log.Write("\no  Extracting molecule data from "+comfile+".com ...")
    MOLSPEC = getinData(comfile, log)
    CONFSPEC = MOLSPEC
		
    # Perform Gaussian TS optimizations - need to automate the job type here somehow!!!
    JOB = JobSpec("Gaussian")
    JOB.JOBTYPE = MOLSPEC.JOBTYPE
    if hasattr(MOLSPEC, "MEMREQ"): JOB.Mem = str(MOLSPEC.MEMREQ)
    if hasattr(MOLSPEC, "NPROC"):
        JOB.NPROC = str(MOLSPEC.NPROC)
		

	# Read the Monte Carlo parameters (from the FullMonte outfile)
    #log.Write("\no  Extracting conformational search parameters from "+filein+"_fm.out"+" ...")
	#SEARCHPARAMS = getParams(MOLSPEC, filein+"_fm.out",log)
	
	
    # DEFINE NEW VARIABLES
    OLDSEARCH= JobSpec("Gaussian")
    OLDSEARCH.ENERGY = []
    OLDSEARCH.CARTESIANS = []
    OLDSEARCH.NAME = []
    OLDSEARCH.TIMESFOUND = []
    OLDSEARCH.GLOBMIN = 999.9
    CSEARCH.GLOBMIN = 999.9
    CSEARCH.CPU = []
    CSEARCH.ALLCPU = []
    CSEARCH.TORVAL = []
    CSEARCH.NSAVED = 0
    CSEARCH.ENERGY = []
	

    if os.path.isfile(filein+"_fm.out"):
        # Parse all saved constrained optimizations and create new input files for TS optimization
        log.Write("\no  Extracting previously optimized structures from "+filein+"_fm.out ...")
        conffile = open(filein+"_fm.out","r") 
        conflines = conffile.readlines()
        

        for line in conflines:
            if line.find("CONFORMER") > -1: 
                if float(line.split()[3]) < OLDSEARCH.GLOBMIN: OLDSEARCH.GLOBMIN = float(line.split()[3])
                                    
        for line in conflines:
            if line.find("CONFORMER") > -1:
                if (float(line.split()[3])-OLDSEARCH.GLOBMIN)*2625.5 < CSEARCH.ENERGYRANGE:
                    OLDSEARCH.NAME.append(line.split()[2])
                    OLDSEARCH.ENERGY.append(float(line.split()[3]))
                    OLDSEARCH.TIMESFOUND.append(int(line.split()[4]))
                    if float(line.split()[3]) < OLDSEARCH.GLOBMIN: OLDSEARCH.GLOBMIN = float(line.split()[3])
                
                
        for i in range(0,len(conflines)):
            if conflines[i].find("CONFORMER")> -1:
                if (float(conflines[i].split()[3])-OLDSEARCH.GLOBMIN)*2625.5 < CSEARCH.ENERGYRANGE:
                    if conflines[i+1].find("Standard orientation:")> -1: skip = 6
                    else: skip = 10
                    j = 0
                    
                    MOLSPEC.CARTESIANS = []
                    while conflines[i+skip+j].find("----------------") == -1:
                        if len(conflines[i+skip+j].split()) == 6:
                            MOLSPEC.CARTESIANS.append([float(conflines[i+skip+j].split()[3]), float(conflines[i+skip+j].split()[4]), float(conflines[i+skip+j].split()[5])])
                            j = j +1
                    OLDSEARCH.CARTESIANS.append(MOLSPEC.CARTESIANS)


    if os.path.isfile(filein+"_fm.sdf"):
        # Parse all saved constrained optimizations and create new input files for TS optimization
        log.Write("\no  Extracting previously optimized structures from "+filein+"_fm.sdf ...")
        conffile = open(filein+"_fm.sdf","r")
        conflines = conffile.readlines()
            
        for line in conflines:
            if line.find("E=") > -1:
                if float(line.split(("="))[1]) < OLDSEARCH.GLOBMIN: OLDSEARCH.GLOBMIN = float(line.split(("="))[1])
        
        for i in range(0,len(conflines)):
            if conflines[i].find("E=") > -1:
                if (float(conflines[i].split(("="))[1])-OLDSEARCH.GLOBMIN)*4.184 < CSEARCH.ENERGYRANGE:
                    OLDSEARCH.NAME.append(conflines[i-1].rstrip().lstrip().split()[0])
                    OLDSEARCH.ENERGY.append(float(conflines[i].split(("="))[1]))
                    OLDSEARCH.TIMESFOUND.append(1)
                    MOLSPEC.CARTESIANS = []
                    j = 0
                    while conflines[i+3+j].find("0  0  0  0") > -1:
                        MOLSPEC.CARTESIANS.append([float(conflines[i+3+j].split()[0]), float(conflines[i+3+j].split()[1]), float(conflines[i+3+j].split()[2])])
                        j = j +1
                    OLDSEARCH.CARTESIANS.append(MOLSPEC.CARTESIANS)
				
				
	if not os.path.isfile(qmtgz):
		log.Write("\no  Creating QM optimization for structures within "+str(CSEARCH.ENERGYRANGE)+" kJ/mol of global minimum ...")
		comtar = tarfile.open(qmtgz, "w|gz")
		
		print OLDSEARCH.NAME
		for i in range(0, len(OLDSEARCH.NAME)):
			# Read coordinates from constained optimization
			MOLSPEC.NAME = OLDSEARCH.NAME[i]+"_"+qmsuffix
			MOLSPEC.CARTESIANS =  OLDSEARCH.CARTESIANS[i]
			
		
			# Write new input file for QM optimization
			writeInput(JOB, MOLSPEC)
			comtar.add(MOLSPEC.NAME+".com")	
			os.remove(MOLSPEC.NAME+".com")
			
		comtar.close()	
	#except: log.Fatal("\nFATAL ERROR: Conformers from file [ %s ] cannot be read"%(filein+".out"))	 
		

	# Now update tar to point to TS optimization files
	tar = tarfile.open(qmtgz, "r:gz")
	names = tar.getnames()

	for name in names: CSEARCH.TODO.append(name.split(".")[0])
	CSEARCH.NJOBS = len(names)
	
	# Begin conformational search - print out the various parameters to terminal and log file
	start = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
	#asciiArt(start)

	for todo in CSEARCH.TODO:
		if os.path.isfile(todo+".log")==1 or os.path.isfile(todo+".out")==1 or os.path.isfile(todo+".chk")==1:
			
			if os.path.isfile(todo+".log")==1 or os.path.isfile(todo+".out")==1:
				CONFSPEC.NAME = todo
				if isJobFinished(JOB, CONFSPEC) > 0:
					CSEARCH.DONE = CSEARCH.DONE + 1
					CSEARCH.DONELIST.append(todo)
						
				else:
					CSEARCH.RUNNING = CSEARCH.RUNNING + 1
					CSEARCH.RUNNINGLIST.append(todo)
				
			else:
				CSEARCH.RUNNING = CSEARCH.RUNNING + 1
				CSEARCH.RUNNINGLIST.append(todo)
								
	# Summary at the beginning				
	log.Write("\no  "+str(CSEARCH.NJOBS)+" optimizations to be performed ...")
	log.Write("   "+str(CSEARCH.DONE)+" have already terminated ...")
	#for item in CSEARCH.DONELIST: print item,
	
	log.Write("   "+str(CSEARCH.RUNNING)+" are presumably still running ...\n")
	#for item in CSEARCH.RUNNINGLIST: print item,
	

	#If there are still jobs to be run:
	if CSEARCH.DONE < CSEARCH.NJOBS:
		
		#run job not already in the list of running jobs or done jobs
		for job in CSEARCH.TODO:
			already = 0
			for running in CSEARCH.RUNNINGLIST:
				if job == running: already = 1
			for finished in CSEARCH.DONELIST:
				if job == finished: already = 1
				
			if already == 0:
				tar.extract(job+".com", path="")
				MOLSPEC.NAME = job
				log.Write("o  Please submit "+job+".com for QM optimization ")	
				CSEARCH.RUNNING = CSEARCH.RUNNING + 1
				CSEARCH.RUNNINGLIST.append(job)
					
	print "\no  Checking QM conformers for duplicates..."
	#Filter after optimization - check for high energy, duplicate of previously saved, isomerization and presence of imaginary frequency if applicable.
	if CSEARCH.DONE > 0:
		for outfile in CSEARCH.DONELIST:
			CONFSPEC.NAME = outfile
			
			# Make sure optimization is complete 
			exitjob = isJobFinished(JOB, CONFSPEC)
			while exitjob == 0: exitjob = isJobFinished(JOB, CONFSPEC)
			
			
			# Successful termination - extract details
			if exitjob == 1:
				CONFSPEC =  getoutData(CONFSPEC)

				CSEARCH.ALLCPU.append(CONFSPEC.CPU)
				
				
				#Check whether the molecule has isomerized - usually this is undesirable so filter out structural isomers
				concheck = checkconn(CONFSPEC, MOLSPEC, OLDSEARCH, SEARCHPARAMS)
				if concheck[0] == 0: CONFSPEC.CONNECTIVITY = MOLSPEC.CONNECTIVITY; isomerize = 0
				else: isomerize = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected: "+concheck[1]+concheck[2]+" has broken from "+concheck[3]+concheck[4]).ljust(50))
			
				#If it is a TS, check the forming/breaking bond of interest
				if len(MOLSPEC.CONSTRAINED) > 0:
					const0 = calcdist(MOLSPEC.CONSTRAINED[0][0],MOLSPEC.CONSTRAINED[1][0],MOLSPEC.CARTESIANS)
					const1 = calcdist(MOLSPEC.CONSTRAINED[0][0],MOLSPEC.CONSTRAINED[1][0],CONFSPEC.CARTESIANS)				
					#print const0, const1, const1-const0
					if math.fabs(const1 - const0) > 0.3: 
						isomerize = 1	
						log.Write("   "+(CONFSPEC.NAME+" is rejected: "+str(MOLSPEC.CONSTRAINED[0][0])+" has broken/formed with "+str(MOLSPEC.CONSTRAINED[1][0])).ljust(50))
					
				# Check for high energy
				if ((CONFSPEC.ENERGY - CSEARCH.GLOBMIN)*2625.5) < SEARCHPARAMS.DEMX * 2: toohigh = 0
				else: toohigh = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected due to high energy ... ").ljust(50))
				
				
				#check for an imaginary frequency
				negfreq = 1
				if len(CONFSPEC.FREQS)>1: 	
					if CONFSPEC.FREQS[0] > 0.0: negfreq = 0; log.Write("\n   "+(CONFSPEC.NAME+" is rejected due to no imaginary freqs ... ").ljust(50))
					if CONFSPEC.FREQS[1] < -30.0: negfreq = 2; log.Write("\n   "+(CONFSPEC.NAME+" is rejected due to more than one imaginary freqs ... ").ljust(50))
				else: conffreq=[0.00]
				#log.Write("o  "+(CONFSPEC.NAME+": unable to read freqs ... ").ljust(50))
				
				# Save or discard the optimized structure - reject if higher than global minimum by DEMX kJ/mol
				if toohigh ==0 and isomerize == 0 and negfreq == 1: 
					samecheck=0
					
					# Clustering - check if the energy and geometry match a previously saved structure
					esame = 0.5
					for j in range(0, CSEARCH.NSAVED):
						if (CONFSPEC.ENERGY-CSEARCH.GLOBMIN)*2625.5 < esame*-1.0: break
						if abs((CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5) < 0.5:
                                                        #print "   "+CONFSPEC.NAME+" cf "+CSEARCH.NAME[j]+" = "+str((CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5)
							if checkSame(CONFSPEC, CSEARCH, SEARCHPARAMS, j) > 0 or checkSame(makemirror(CONFSPEC), CSEARCH, SEARCHPARAMS, j) > 0:
								log.Write("   "+(CONFSPEC.NAME+" is a duplicate of conformer "+CSEARCH.NAME[j]+" ... ").ljust(50))
								CSEARCH.TIMESFOUND[j] = CSEARCH.TIMESFOUND[j] + 1
								CSEARCH.NREJECT = CSEARCH.NREJECT + 1
								samecheck = samecheck + 1
								break
								
					
					# Unique conformation with low energy!				
					if samecheck == 0:
						# Check to see if conformer is now the lowest in energy
						if CONFSPEC.ENERGY < CSEARCH.GLOBMIN:
							CSEARCH.GLOBMIN = CONFSPEC.ENERGY
							log.Write("   "+(CONFSPEC.NAME+" is a new Global Minimum!").ljust(80)+("E = "+str(CSEARCH.GLOBMIN)).rjust(rightcol))
						else : log.Write("   "+(CONFSPEC.NAME+" is saved").ljust(80)+("E = "+str(CONFSPEC.ENERGY)).rjust(rightcol))
						AddConformer(CSEARCH, CONFSPEC)
						
						
						# Check when last low-energy conformer was found
						if (CONFSPEC.ENERGY - CSEARCH.GLOBMIN) * 2625.5 < SEARCHPARAMS.EWIN: CSEARCH.LASTFOUND = CSEARCH.STEP*SEARCHPARAMS.POOL
							
							
				# Rejection - discard
				else: CSEARCH.NREJECT = CSEARCH.NREJECT + 1
		
			# Unsuccessful termination - discard
			else: 
				log.Write("\n   Unsuccessful optimization of "+CONFSPEC.NAME+".out ...")
				
				
				# Update statistics
				CSEARCH.NFAILED = CSEARCH.NFAILED + 1


		#Order the low energy conformers by energy
		OrderConfs(CSEARCH, SEARCHPARAMS, start, log)
		CSEARCH.GLOBMIN = CSEARCH.ENERGY[0]
	
		# Final Summary of Full Monte search
		WriteQMSummary(CSEARCH, OLDSEARCH, SEARCHPARAMS, start, log, qmsuffix)	
	
		makeGVformat(filein, MOLSPEC, CSEARCH, SEARCHPARAMS, "qm")
		makePDBformat(filein, MOLSPEC, CSEARCH, "qm")
	
	else: log.Write("\no  No complete jobs to report on ...")
	#Write the saved conformers to the log file for viewing in GaussView

	end = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
	asciiArt(end)
	log.Write(normaltermination)
	log.Finalize()	

