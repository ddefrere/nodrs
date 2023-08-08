/* uvklib.c read uko5v.dat, uko9v.dat ... ukm2i.dat
 * write flux to flam.uvklib, and rms to rms.uvklib
 * lamA (1150:25000:5) in row 0, rows 1-131 stellar spectra
 *
 * compile with gcc -lm -lc -o uvklib uvklib.c (type make in this directory)
 * and run to write the files flam.uvklib and rms.uvklib
 */
#include </usr/5include/stdio.h>
#include </usr/5include/stdlib.h>

#define WAVE 4771
#define STAR 131      /* 0row wavelength, 1-131rows spectra */

#define DEBUG 1

float wave[WAVE], flam[WAVE], rms[WAVE];

extern int strcat();

void main()
{
  FILE    *fi, *fo, *fr;
  int     w=WAVE, s=STAR+1, i, k;
  char    string[111], fname[14];
  char    name[STAR][8] = {"o5v"    ,"o9v"    ,"b0v"    ,"b1v"    ,"b3v", 
			   "b57v"   ,"b8v"    ,"b9v"    ,"a0v"    ,"a2v",
			   "a3v"    ,"a5v"    ,"a7v"    ,"f0v"    ,"f2v",
			   "f5v"    ,"wf5v"   ,"f6v"    ,"rf6v"   ,"f8v",
                           "wf8v"   ,"rf8v"   ,"g0v"    ,"wg0v"   ,"rg0v",
                           "g2v"    ,"g5v"    ,"wg5v"   ,"rg5v"   ,"g8v",
                           "k0v"    ,"rk0v"   ,"k2v"    ,"k3v"    ,"k4v",
                           "k5v"    ,"k7v"    ,"m0v"    ,"m1v"    ,"m2v",
                           "m2.5v"  ,"m3v"    ,"m4v"    ,"m5v"    ,"m6v",
			   "b2iv"   ,"b6iv"   ,"a0iv"   ,"a47iv"  ,"f02iv",
                           "f5iv"   ,"f8iv"   ,"g0iv"   ,"g2iv"   ,"g5iv",
                           "g8iv"   ,"k0iv"   ,"k1iv"   ,"k3iv"   ,"o8iii",
                           "b12iii" ,"b3iii"  ,"b5iii"  ,"b9iii"  ,"a0iii",
                           "a3iii"  ,"a5iii"  ,"a7iii"  ,"f0iii"  ,"f2iii",
                           "f5iii"  ,"g0iii"  ,"g5iii"  ,"wg5iii" ,"rg5iii",
                           "g8iii"  ,"wg8iii" ,"k0iii"  ,"wk0iii" ,"rk0iii",
                           "k1iii"  ,"wk1iii" ,"rk1iii" ,"k2iii"  ,"wk2iii",
                           "rk2iii" ,"k3iii"  ,"wk3iii" ,"rk3iii" ,"k4iii",
                           "wk4iii" ,"rk4iii" ,"k5iii"  ,"rk5iii" ,"m0iii",
                           "m1iii"  ,"m2iii"  ,"m3iii"  ,"m4iii"  ,"m5iii",
                           "m6iii"  ,"m7iii"  ,"m8iii"  ,"m9iii"  ,"m10iii",
                           "b2ii"   ,"b5ii"   ,"f0ii"   ,"f2ii"   ,"g5ii",
                           "k01ii"  ,"k34ii"  ,"m3ii"   ,"b0i"    ,"b1i",
                           "b3i"    ,"b5i"    ,"b8i"    ,"a0i"    ,"a2i",
                           "f0i"    ,"f5i"    ,"f8i"    ,"g0i"    ,"g2i",
                           "g5i"    ,"g8i"    ,"k2i"    ,"k3i"  ,"k4i", "m2i"};

  if( (fo = fopen("flam.uvklib", "w+b")) == NULL) {
      printf("Cannot open flam.uvklib\n");
      exit(1);
  }
  if( (fr = fopen("rms.uvklib", "w+b")) == NULL) {
      printf("Cannot open rms.uvklib\n");
      exit(1);
  }
  fwrite(&w, sizeof(int), 1, fo);    fwrite(&w, sizeof(int), 1, fr);
  fwrite(&s, sizeof(int), 1, fo);    fwrite(&s, sizeof(int), 1, fr);

  for(k=0; k<STAR; k++) { strcpy(fname, "uk");
    strcat(fname, name[k]);strcat(fname, ".dat");printf("%i %s\n", k+1,fname);
    if( (fi = fopen(fname, "r")) == NULL) {
    printf("Cannot open %s\n", fname);     exit(1); }
    fgets(string, 110, fi); if(DEBUG) printf("%s", string);
    fgets(string, 110, fi); if(DEBUG) printf("%s", string);
    fgets(string, 110, fi); if(DEBUG) printf("%s", string);
    for(i=0; i<WAVE; i++) {
      fgets(string, 110, fi); if(DEBUG>2) printf("%s", string);
      sscanf(string, "%f %f %f", &wave[i], &flam[i], &rms[i]);
if(DEBUG>1)printf("%i %7.1f %10.7f %10.7f\n",i,wave[i],flam[i],rms[i]); }
    fclose(fi);  
    if(k==0) {fwrite(wave, sizeof(float), WAVE, fo); 
	      fwrite(wave, sizeof(float), WAVE, fr);}
    fwrite(flam,sizeof(float),WAVE,fo);
    fwrite(rms, sizeof(float),WAVE,fr);
  }
		       
  fclose(fo); fclose(fr);
}

