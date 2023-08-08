#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <math.h>
#include <string.h>
#include "fitsio.h"
//#include "export.h"

// define constants
const int KEYLEN = 40;      // Maximum keyword length 
const int FIDLEN = 6;       // Maximum keyword length

// Define keyword structure
struct lbti_key {
	 char  objname[KEYLEN], timeobs[KEYLEN], instrume[KEYLEN], nomicfw1[KEYLEN], nomicfw2[KEYLEN], lmirfw1[KEYLEN], lmirfw2[KEYLEN], lmirfw3[KEYLEN], lmirfw4[KEYLEN];
	 float exptime, lbt_lxos, lbt_lyos, lbt_rxos, lbt_ryos, lambda;
	 int   fileid, datatype, obstype, naxis1, naxis2, naxis3, pid;
};

void lbti_masterlog(int argc, void *argv[])
{

	//lbti_masterlog.c
	//Last updated 23 Oct 2014
	//Created by Denis Defrère
	//University of Arizona
	//Steward Observatory -- ddefrere@email.arizona.edu

	//To compile this code via IDL:
    //c_path = '/Users/Denis/IDLWorkspace82/nodrs/c/'
    //file_c = 'lbti_masterlog'
    //MAKE_DLL, file_c, file_c, file_c, INPUT_DIRECTORY=c_path, OUTPUT_DIRECTORY=c_path, EXTRA_LFLAGS='-lcfitsio'

	//-------------------- DECLARE -------------------- 
	
	//Definitions that shouldn't be touched
    fitsfile *fptr;         
    char *comment=NULL, text[FIDLEN];
	char **inputfiles, *ptr;
    int i, j, k, nodid, chpid, deltaj, deltajmax, status;
    float lbt_rxos0, lbt_ryos0, lbt_lxos0, lbt_lyos0; 
    struct lbti_key keywords;

    time_t mytime;
    mytime = time(NULL);
    
    //------------------------------------------------- 
    
    
	//-------------------- INPUT --------------------

	//Here are the variables passed by the IDL wrapper.
	//Note that IDL is not neccessary--you can create your own wrapper if you want.
	//Note: all scalars are interpreted below as passed-by-value, so lbti_masterlog.c won't change them
	char *masterlogfile = (char *) argv[0]; //Masterlog file
	char *inputfilelist = (char *) argv[1]; //Concatenated strings of input files from IDL
	int numfiles = *(int *) argv[2];        //number of input file strings to expect  
	int verbose_flag = *(int *) argv[3];   //Print a lot of info
	
	//-----------------------------------------------
	
	
	
    //-------------------- REWORK SOME INPUT --------------------

	//Here we reinterpret the list of input files given by IDL.
	//IDL just gives us a series of concatenated strings, separated by a null terminating
	//character.  It's more elegant to deal with an array of strings, so we manipulate it a bit.

	//Find length of the longest filename, deltajmax
	j=0;
	deltajmax=0;
	for(i=0;i<numfiles;i++) {
		deltaj = (1 + strcspn (&inputfilelist[j],"\0"));
		j += deltaj;
		if(deltaj > deltajmax) deltajmax=deltaj;
	}
	deltajmax += 10;

	//Now create an array of inputfiles, each with length deltajmax
	inputfiles = (char **) malloc(numfiles * sizeof(char *));   // allocate storage for an array of pointers   
	if (inputfiles==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating inputfile array.\n"); exit (1);}
	for (i = 0; i < numfiles; i++) {
		inputfiles[i] = (char *) malloc(deltajmax * sizeof(char));   // for each pointer, allocate storage for an array of chars
		if (inputfiles[i] == NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating inputfile array.\n"); exit (1);} 
	}

	//Now print the concatenated strings into an array of strings
	j=0;
	for(i=0;i<numfiles;i++) {
		sprintf(inputfiles[i],"%s",&inputfilelist[j]);    
		j += (1 + strcspn (&inputfilelist[j],"\0"));
	}
	
	//-------------------------------------------------------------    
	
	    
    
    //-------------------- COMPUTATION -------------------- 
    // Open output file
	FILE *f = fopen(masterlogfile, "w");
	if (f == NULL)
	{
		printf("LBTI_MASTERLOG.C: error opening file!\n");
		exit(1);
	}
	//printf("Writing in masterlog file: %s", masterlogfile);
	fprintf(f,"Masterlog created on %s ", ctime(&mytime));
    fprintf(f,"ID;name;time_obs;lam_eff[m];int_time[s];NOD_ID;CHP_ID;OBSTYPE;DATATYPE;N_XPIX;N_YPIX;N_FRAME;PID\n");
    
    // Init variables
    lbt_lxos0 = -1;
    lbt_lyos0 = -1;
    lbt_rxos0 = -1;
    lbt_lyos0 = -1;
    nodid = -1;
    
    // Loop over the files  
    for(i=0;i<numfiles;i++) {
    	//Open file (make sure comment is NULL and status is 0)
    	status = 0;
    	comment=NULL;
		fits_open_file(&fptr, inputfiles[i], READONLY, &status);
		
		// Derive file ID (6 digits for NOMIC, 5 for LMIRCam)
		memset(text, '\0', sizeof(text));
		deltaj = 1 + strcspn(inputfiles[i],".");
		strncpy(text, inputfiles[i]+deltaj-7, 6);
		if (strcspn(text,"_") == 0)
		{
			strncpy(text, inputfiles[i]+deltaj-6, 5);
			text[5] = 0;  // Add null terminator when only 5 digits
		}
		
		// Convert file ID to long
		keywords.fileid = strtol(text, &ptr, 10); 
		
		// Read main file keywords 
		fits_read_key(fptr, TSTRING, "OBJNAME" , &keywords.objname,  NULL, &status);
		fits_read_key(fptr, TSTRING, "TIME-OBS", &keywords.timeobs,  NULL, &status);
		fits_read_key(fptr, TFLOAT,  "EXPTIME" , &keywords.exptime,  NULL, &status);
		fits_read_key(fptr, TINT,    "NAXIS1"  , &keywords.naxis1,   NULL, &status);   
		fits_read_key(fptr, TINT,    "NAXIS2"  , &keywords.naxis2,   NULL, &status);  
		fits_read_key(fptr, TINT,    "NAXIS3"  , &keywords.naxis3,   NULL, &status); 
		fits_read_key(fptr, TINT,    "OBSTYPE" , &keywords.obstype,  NULL, &status);	
		if (status != NULL) 
		{
		 	keywords.obstype = 0;
		}
		fits_read_key(fptr, TINT,    "DATATYPE", &keywords.datatype, NULL, &status);	
		if (status != NULL)  
		{
		 	keywords.datatype = 0;
		}
		fits_read_key(fptr, TINT,    "PID"     , &keywords.pid,      NULL, &status);
		if (status != NULL) 
		{
			keywords.pid = -1;
		}
				
		// Compute NOD_ID
		fits_read_key(fptr, TFLOAT,  "LBT_LXOS", &keywords.lbt_lxos, NULL, &status);   
		fits_read_key(fptr, TFLOAT,  "LBT_LYOS", &keywords.lbt_lyos, NULL, &status);
		fits_read_key(fptr, TFLOAT,  "LBT_RXOS", &keywords.lbt_rxos, NULL, &status);
		fits_read_key(fptr, TFLOAT,  "LBT_RYOS", &keywords.lbt_ryos, NULL, &status);
		if ((keywords.lbt_lxos != lbt_lxos0 || keywords.lbt_lyos != lbt_lyos0) 
			|| (keywords.lbt_rxos != lbt_rxos0 || keywords.lbt_ryos != lbt_ryos0)) 
			nodid += 1;
   	    lbt_lxos0 = keywords.lbt_lxos; 
   	    lbt_lyos0 = keywords.lbt_lyos;
    	lbt_rxos0 = keywords.lbt_rxos;
    	lbt_ryos0 = keywords.lbt_ryos;
    	
    	// Compute CHOP_ID (not yet defined)
    	chpid = 0;
    	
    	// Compute wavelength and bandwidth (not working yet)
		// fits_read_key(fptr, TSTRING, "INSTRUME", &keywords.instrume, NULL, &status);
		// if (strcmp(keywords.instrume, "NOMIC") == 0) 
		// {
  	    //   	fits_read_key(fptr, TSTRING, "NOMICFW1", &keywords.nomicfw1, NULL, &status);
		//	fits_read_key(fptr, TSTRING, "NOMICFW2", &keywords.nomicfw2, NULL, &status);
		//	//printf(" NOMIC!");
		//} 
		//else if (strcmp(keywords.instrume, "LMIRCAM") == 0)
		//{
 		//	fits_read_key(fptr, TSTRING, "LMIR_FW1", &keywords.lmirfw1, NULL, &status); 
		//	fits_read_key(fptr, TSTRING, "LMIR_FW2", &keywords.lmirfw2, NULL, &status);  
		//  fits_read_key(fptr, TSTRING, "LMIR_FW3", &keywords.lmirfw3, NULL, &status);
		//   fits_read_key(fptr, TSTRING, "LMIR_FW4", &keywords.lmirfw4, NULL, &status);
		//}
		/* more else if clauses */
		//else /* default: */
		//{
		//	 printf("Unknown instrument");
		//}
		keywords.lambda=0.0;  

		//if (verbose_flag >= 2) {
	    //  printf("    - %s\n",inputfiles[i]);
	    //	printf(" File %d", nodid);
  		// }
		
		//CLose file and print error
		fits_close_file(fptr, &status);
		if (status)          /* print any error messages */
			fits_report_error(stderr, status);
				
		// Print information to masterlog file
		
		/* print some text */
		fprintf(f,"%06d;%s;%s;%2.2E;%1.4f;%d;%d;%d;%d;%d;%d;%d;%d\n", keywords.fileid, keywords.objname,
				keywords.timeobs, keywords.lambda, keywords.exptime, nodid, chpid,
				keywords.datatype, keywords.obstype, keywords.naxis1, keywords.naxis2,
				keywords.naxis3, keywords.pid);
    }
    
    // Close text file
    fclose(f);
    
    //---------------------------------------------------
}
