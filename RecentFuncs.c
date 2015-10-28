/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright  1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: RecentFuncs.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

/* The use of this software (except that by Harald T. Alvestrand, which is described in 'udping.c')
   and that by James Knight, which is described in 'seqio.c' is bound by the notice that appears
   in the file 'tacg.h' which should accompany this file.  In the event that 'tacg.h' is not bundled
   with this file, please contact the author.
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pcre.h>
#include <stddef.h>
#include <sys/time.h>
#include <time.h>   /* for the timing for regex study routine */
#include <locale.h> /* for the timing for regex study routine */
#define LOOPREPEAT 30 // 3 times the # of subsequences at each match, if interested
#define OVECCOUNT 300    /* for pcre should be a multiple of 3 */

#include "seqio.h"
#include "tacg.h"

//Function to spit out debug error messages at the given level of verbosity (1-3) Use like this:
// tacg_debug(1, __FILE__, __LINE__, "This is the error message I want printed:  Please do it well!");

void tacg_debug(int dbg_lvl, char *fle_nm, int lin_nbr, char *err_msg) {
    if (F.Verbose >= dbg_lvl) {
        fprintf(stderr, "DEBUG:[%s]:%d: %s\n", fle_nm, lin_nbr, err_msg);
    }
}



/* tacg_Time_Ops is a developer's function to test the speed at which a
	system can execute different types of operations:
	floating point ops - we'll start here.
	integer ops
	memory ops
 */

double tacg_TimeOps (void) {

    double result=0, result2=0, result3=0; //, djunk=0;
    long max=31623, maxsq=0, counter=0, counter2=0, iresult, iresult2;
    double t1, t2, t3=2.3453, sum_all, fjunk, elapsed=0, cnt_time=0, dcounter, dcounter2, logtarget=38736483.345234;
    //int flarr_sz = 1000000;

    struct timeval tp;
    unsigned int seed;
    int i, foo, bar, step_sz, rtn, *ptr; //
    int step[] = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096}; /*  sizes for malloc */
    maxsq = max * max;
    printf("\nstarting TimeOps now..\n");
    //float flarr[flarr_sz];

    // fill a 1,000,000 el array with random #s
//     printf("start filling 1,000,000 els with random #s\n");
//     for (i=0; i<flarr_sz; i++) {
//         flarr[i] = (float)rand();
//     }
//     printf("finished\n");

    /************************************************************************
     *                             random()                                 *
     ************************************************************************/
//get the time for the loop itself"
    bar = 3000000;
    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (counter=1; counter<bar; counter++) {
        //do nothing
    }
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    t3=(t2-t1);
    printf("\nrandom()\nTime for 30M count loop: %f sec\n",t3);
    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    seed=(unsigned int)t1;
    srandom(seed);
    for (counter=1; counter<bar; counter++) {
        //seed=(unsigned int)tp.tv_sec+(1.e-6)*tp.tv_usec;
        //srandom(seed);
        iresult = random();
        //printf("%ld ", iresult);
    }
    //djunk = result * 1.12432;
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    elapsed=(t2-t1-t3);
    printf("Time for %2d M each log & sqrt ops (- loop time): %6.3f \n for %6.3f M Tr_Op/s\n",(bar/1000000), elapsed, (((2*bar)/elapsed)/1000000));



    /************************************************************************
     *                             Count to 2B                              *
     ************************************************************************/
// count to 2B
    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (counter=0; counter<max; counter++) {
        for (counter2=0; counter2<max; counter2++) {
            // do nothing
        }
    }
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    cnt_time=t2-t1;
    printf("\nInteger Ops\nTime to int count to %ld: %f sec\n", maxsq, cnt_time);

// integer ops
    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (counter=0; counter<max; counter++) {
        for (counter2=0; counter2<max; counter2++) {
            iresult = (counter2 - (counter%6)) << 4;
            iresult2 +=  (counter2 * (iresult/counter2));
            //iresult += 7;  1sub, 1 add, 2 assigns, 1 mul, 1 div, 1 mod
        }
    }
    fjunk = (double)iresult2 * 343.984; // make sure that result gets used.
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    elapsed=(t2-t1)-cnt_time;
    printf("Time for %f integer ops (1sub, 1 add, 2 assigns, 1 mul, 1 div, 1 mod)\n(- loop time): %f sec: %6.2f Mop/s\n",(7*(double)maxsq), elapsed, (((7*(double)maxsq)/elapsed)/1000000));


    /************************************************************************
     *                        floating point ops                            *
     ************************************************************************/
    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (dcounter=1; dcounter<max; dcounter++) {
        for (dcounter2=1; dcounter2<max; dcounter2++) {
            //do nothing
        }
    }
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    t3 = t2-t1;
    printf("\nFloating Point Ops\nTime to float count to  %ld: %f s\n",maxsq, t3);

    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (dcounter=1; dcounter<max; dcounter++) {
        for (dcounter2=1; dcounter2<max; dcounter2++) {
            result = dcounter / (dcounter2 - 67) ;
            result3 = dcounter * ((result + 876) / dcounter);
            //result3 = result / result2; 2 assigns, 2 divs, 1 mul, 1 add 1 sub
        }
    }
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    elapsed=(t2-t1)-t3;
    fjunk = result3 * 343.984; // make sure that result gets used.
    printf("Time for 7B floating point ops (2 assigns, 2 divs, 1 mul, 1 add 1 sub)\n(- loop time): %7.3f sec: %6.3f Mflop/s\n", elapsed, ((7*(double)(max*max))/elapsed)/1000000 );


    /************************************************************************
     *                             transendental ops                        *
     ************************************************************************/
//get the time for the loop itself"
    bar = 30000000;
    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (counter=1; counter<bar; counter++) {
        //do nothing
    }
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    t3=(t2-t1);
    printf("\nTranscendentals\nTime for 30M count loop: %f sec\n",t3);

    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (counter=1; counter<bar; counter++) {
        result = log(logtarget);
        result = sqrt(logtarget);
    }
    //djunk = result * 1.12432;
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    elapsed=(t2-t1-t3);
    printf("Time for %2d M each log & sqrt ops (- loop time): %6.3f \n for %6.3f M Tr_Op/s\n",(bar/1000000), elapsed, (((2*bar)/elapsed)/1000000));


    /************************************************************************
     *                Memory Ops: callocs and frees                         *
     ************************************************************************/
    printf("\nMemory ops - calloc() and free()\n");
    for (step_sz=0; step_sz<11; step_sz++) {
        rtn=gettimeofday(&tp, NULL);
        t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        for (counter=1; counter<1000000; counter++) {
            ptr = (int *)calloc(step[step_sz], sizeof(int));
            if (ptr == NULL) BadMem("Can't init mem for ptr", 1);
            free(ptr);
        }
        rtn=gettimeofday(&tp, NULL);
        t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        elapsed=(t2-t1);
        printf("Time for 1M callocs and frees of size %5d x int(%d bytes): %6.3f sec\n",step[step_sz], (int)sizeof(int), elapsed);
    }

    /************************************************************************
     *                Memory Ops: mallocs and frees                         *
     ************************************************************************/
    printf("\nMemory ops - malloc() and free()\n");
    for (step_sz=0; step_sz<11; step_sz++) {
        rtn=gettimeofday(&tp, NULL);
        t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        for (counter=1; counter<1000000; counter++) {
            ptr = (int *)malloc(step[step_sz]);
            if (ptr == NULL) BadMem("Can't init mem for ptr", 1);
            free(ptr);
        }
        rtn=gettimeofday(&tp, NULL);
        t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        elapsed=(t2-t1);
        printf("Time for 1M mallocs and frees of size %5d x int(%d bytes): %6.3f sec\n",step[step_sz], (int) sizeof(int), elapsed);
    }



    /************************************************************************
     *                Memory Ops: assignment & accesses                     *
     ************************************************************************/
    printf("\nMemory ops - assignment and access\n");
    // assign a 100M array
    bar = 10000000;
    ptr = (int *)calloc(bar, sizeof(int));
    if (ptr == NULL) BadMem("Can't init mem for ptr", 1);
    // for (bar=0; bar<1000000; bar++){	ptr[bar] = bar; }; //and fill it

    //count the empty loop
    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (counter=1; counter<bar; counter++, counter2=counter-1) {
        counter2 = counter-1;
    }
    foo = counter2;
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    elapsed=(t2-t1);
    printf("Time to int count to %2d M with 1 extra add: %6.3f\n", bar/1000000, elapsed);


    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    for (counter=1; counter<bar; counter++, counter2=counter-1) {
        ptr[counter] = counter;
        ptr[counter2] = ptr[counter];
    }
    rtn=gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    foo = ptr[counter2] - 17;
    elapsed=(t2-t1);
    printf("Time %2d M int accesses & %2d M int assigns: %f sec\n", ((3*bar)/1000000), ((2*bar)/1000000),elapsed);

    sum_all = (double)fjunk - result + (double)ptr[1876871] - (double)foo - result2 + result3 - (double)iresult;
    printf("End of TimeOps\n");

    return(sum_all);
}


/* NmerAnalysis() takes the sequence (and the antipar sequence?) and returns the number of
    degens enountered.  Or it could return a pointer to a struct or array that actually holds
    the data.
*/
long NmerAnalysis(char *sequence, long primes[]) {
    /* similar to Cutting(), but have to make code handle variable length words instead of the default
        hexamer
      NB: in this version note that 'sequence' is the UNPADDED version (only the original sequence)  */
    //DECLARATIONS
    int nmer_len, DCnt=0, basize;
    double nmer_base;
    long CSP=0, clean, i, lp, lag=0, lead, NM_end;
    char nmer[21], *Dgen_ba;  // Dgen_ba is the bit aray for counting degens
    unsigned long int *NM, Key; // tkey;
    short skiplead = 0; // skiplead is a switch that skips the shift/add bits if we hit a degen
    unsigned long a_lag, c_lag,  g_lag, t_lag;

//  ASSIGNMENTS
    clean = 0;
    nmer_len = F.nmer;
    nmer_base = (long)pow(4,nmer_len-1);
    a_lag = 0;
    c_lag = nmer_base;
    g_lag = 2 * nmer_base;
    t_lag = 3 * nmer_base;
    NM_end = F.Seq_Len - nmer_len;  // where in sequence[] to end   for nmer stuff
    //CSP = 0;

    // F.Seq_Len is the ACTUAL UNBUFFERED sequence length

    // INITIAL ALLOCATIONS
    if ((NM = calloc(F.Seq_Len, sizeof(unsigned long int))) == NULL) BadMem("NmerAnalysis:NM", 1);
    // to convert this over to shift/add, have to call nmer_hash only 1x at start or if hit a degen
    // that knocks out the base hash and it has to be re-calc'ed.

    basize = (F.Seq_Len / 8) + 1;
    if ((Dgen_ba = calloc(basize, sizeof(char))) == NULL) BadMem("NmerAnalysis:Dgen_ba", 1);
    // Code starts
    // Initialize the 1st nmer and has
    strncpy(nmer,sequence, nmer_len);
    nmer[nmer_len] = '\0'; //grab the 1st nmer & terminate
    while (nmer_is_not_clean(nmer,nmer_len) && clean < NM_end) {  // check for degens in the defined nmer
        clean++; //has to incr by 1 bc degen may not be in end position; ie have to clean slowly.
        fprintf(stdout,"%.*s @ clean: %ld\n", nmer_len, nmer, clean);
        strncpy(nmer,sequence+clean, nmer_len);
        nmer[nmer_len] = '\0';  //null terminate!
    }
    if (clean >= NM_end) { //only if we've skipped all the way to the bottom
        return(CSP);
    }
    Key = nmer_hash(nmer); // calc hash
    //fprintf(stdout,"%.*s @ clean: %ld and Key = %ld\n", nmer_len, nmer, clean, Key);
    NM[clean] = Key;  // and log it

    // now do the rest as add/shift as much as possible.
    //fprintf(stdout, "nmer_len = %d, a=%ld c=%ld g=%ld t=%ld\n", nmer_len, a_lag, c_lag, g_lag, t_lag);
    i = 0;
    for (CSP=clean+1; CSP<NM_end; CSP++) { // NM_begin+1 cuz we just did the 1st nmer above
//        fprintf(stdout,"nmer = %.*s @ %ld\n", (int)nmer_len, sequence+CSP, CSP);
        lag = 0;
        lp = CSP+nmer_len-1;  /* lp points to the 'lagging' base */
        /* calculate 'lag' value 1st to check if there's incoming degenerates */
        switch (sequence[lp]) {    /* a=0, c=1, g=2, t=3,  */
        case 'a':
            lag = a_lag;
            break;   /* lag values are numeric value * 1024 */
        case 'c':
            lag = c_lag;
            break;
        case 'g':
            lag = g_lag;
            break;
        case 't':
            lag = t_lag;
            break;
            // default detects incoming degens.  should emit a stderr warning and bump the CSP around it.
        default:
            fprintf(stdout, "%c @ %ld, CSP @ %ld \n", sequence[lp], lp, CSP);
            DCnt++;
            BitArray(Dgen_ba, (int)CSP, 1);  // ... set the bit in the Degen log bitarray
            skiplead = 1;  // if degen, don't have to do 2nd part of the add/shift
            // if nmer is NOT free of degens then keep skipping down seq by nmer nucs
            // now skip ahead, looking for next clean nmer, but 1st have to define the bad nmer
            while (nmer_is_not_clean(sequence+lp, nmer_len) && CSP < NM_end) {
                //fprintf(stdout, "loop: %.*s @ CSP= %ld\n", nmer_len, sequence+CSP, CSP);
                CSP += nmer_len;
                lp = CSP;
                strncpy(nmer,sequence+CSP, nmer_len);
                nmer[nmer_len] = '\0';  //null terminate!
            }
            lp = CSP+nmer_len-1; //now bump it to look for the lead
            if (CSP >= NM_end) {
                return(CSP);
            }
            Key = nmer_hash(nmer);
            break;
        }
        if (!skiplead) {
            /* calculate 'lead' value */
            lead = 0;
            switch (sequence[CSP-1]) {     /*  a=0, c=1, g=2, t=3,  */
            case 'a':
                lead = 0;
                break;   /* lead values are numeric value * 1 */
            case 'c':
                lead = 1;
                break;
            case 'g':
                lead = 2;
                break;
            case 't':
                lead = 3;
                break;
            default:
                break;  // bad character detection - already handled above
            }
            /* instead of calling 'hash' on each new hexamer, can calc it incrementally as below */
            Key = ((Key-lead) >> 2) + lag;  //  bit shift version
        }
        fprintf(stdout, "%ld\n", Key);
//        fprintf(stdout, "%ld \t %ld \t %.*s\n", Key, CSP, nmer_len, sequence+CSP);
        //fprintf(stdout, "from shift/add, Key = %ld @ CSP = %ld (%.*s)\n", Key, CSP, nmer_len, sequence+CSP);
        //strncpy(nmer,sequence+CSP,nmer_len); nmer[nmer_len] = '\0';
        //tkey = nmer_hash(nmer);
        //fprintf(stdout, "from nmer_hash, Key = %ld @ CSP = %ld (%s)\n", tkey, CSP, nmer);
        NM[CSP] = Key; // mark each nmer in the sequence
        skiplead = 0;
    }
    // note that in the output below, degens can be distinguished from runs of a's because the
    // degens will go to zero suddenly and run for at least nmer_len and then return to a higher key
    // whereas an incoming run of a's will cause the key to decrease in value gradually, bottoming out at
    // zero.
    for (i=0; i<=F.Seq_Len; i++) {
        //fprintf(stdout, "%ld \t@\t %ld\n", NM[i], i);
    }
    return CSP;
}

// nmer_is_not_clean(nmer, nmer_len) does a scan of the nmer for any degen
// returns 1 if nmer contains a degen, 0 if it's clean, (so true if its dirty)
short nmer_is_not_clean(char *nmer, int nmer_len) {
    short dirty = 0;
    int i=0;
    while (dirty == 0 &&  i < nmer_len) {
        switch (nmer[i]) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
            break;
        default:
            dirty = 1 ;
            break;
        }
        i++;
    }
    //fprintf(stdout, "ninc: %.*s, %d (0 is clean)\n", nmer_len, nmer, dirty);
    return dirty;
}


/*************************  Function nmer_hash  **************************************************
*   nmer_hash() takes a pointer to the n-mer string and generates the                          *
*   integer equivalent.  If there are degens in the string, it will exit with an                 *
*   error code of 0; otherwise it returns the long equiv of the string.  If the verbose option   *
*   is set it will tell you what it bonked on.                                                   *
*************************************************************************************************/
unsigned long int nmer_hash(char *nmer) {
    /* nmer  pointer to begin of array that holds the nmer to be 'hashed'
       num   length of the nmer to be hashed;
    */
// to go to 32-mers, have to use the -std=c89 gcc flag (which allows the use of 'unsigned long long int'
// which would allow hashes to 32-mers, but the c89 flag causes all kinds of errors in the
// compilation of other functions - probably a good thing to corect them but no time to correct them now.

    unsigned long int key[] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304,
                                16777216, 67108864, 268435456, 1073741824
                              };
    //4294967296, 17179869184, 68719476736,
    //            274877906944, 1099511627776 };
//                   9223372036854775807
    unsigned long int sum = 0;
    int ptr, el, num;
    ptr = 0;
    //degen = N_degen = 1;
    //memset(sum,0,sizeof(int)*256);  /* set all of 'sum' to 0 - should work but any faster?*/
    num = strlen(nmer);
    if (num < 2 || num > 20) {
        fprintf (stderr, "nmer_hash(): input string %s is too long (>20) or too short (<2), ", nmer);
        return 0;
    }

    /* Big 'for' loop that calculates the hexamer, gets executed at every RE site, */
    for (el = 0;  el < num;  el++) {       /* and at every overlapping n-mer */
        switch (nmer[el])  { /* pointer passed to function already offset to starting position */
            /* a=0, c=1, g=2, t=3, degenerates are handled below */
        case 'a':
            break;   /* not really needed - (a) = 0 so no change in sum */
        case 'c':
            sum +=   key[el];
            break;
        case 'g':
            sum += 2*key[el];
            break;
        case 't':
            sum += 3*key[el];
            break;

        default:  /* bad character detection */
            if (F.Verbose > 2) {
                fprintf (stderr, "nmer_hash(): don't like %c (at el = %d of %d) of %s! \n",
                         nmer[el], el, num, nmer);
            }
            return 0;
        }  /* end of switch/case statement */
        //degen=N_degen;
    }  /* end of big 'for' loop */
    return sum;
}



/* RegexMatch() searches the sequences for a regex expression read in by ReadRegex() and
   pseudo munged into a 'real' regex by DNA_IUPACtoRegex() for use by the pcre library by
   Phil Hazel.  It differs from the posix regex lib in that the posix lib requires
   additional levels of escaping for the regex to be digested:
   gy(tt|gc)nc{2,3}m -> g[ct]\(tt\|gc\).c\{2,3\}[ca]
   whereas the pcre lib is a bit more sophiticated:
   gy(tt|gc)nc{2,3}m -> g[ct](tt|gc).c{2,3}[ca]
   It also seems to be considerably faster and is MUCH better documented, altho as
   is the case for all regex stuff this is quite a relative term..

   This routine searches BOTH strands for it (actually searches one strand for all the regexes
   then Anti_Par()s the sequence and searches it again with all the regexes (stored in the
   RE struct, same as for the other patterns  */

int RegexMatch(long SEQlen, char *sequence, int NumREs) {

    long sptr, D; /*  search over the boundary, note it if it happens to occur here */
    long ovlapped_seq_len;
    int RxCode = 0, hits = 0, ori, RC, strand, base;
    char *APseq, *SP;
//    vars for pcre regex stuff
    pcre *compiled_regex; /* this is the pcre struct that describes the compiled regex */
//   pcre_extra *extra;
    const char *error;
    int erroffset, rc=0, study_options = 0, options=0;
    int ovector[OVECCOUNT];

    /* note that to buffer correctly, have to add the 2 extra chars */
    ovlapped_seq_len = SEQlen + 2 + (2 * BASE_OVERLAP);

    if ((APseq = calloc(ovlapped_seq_len, sizeof(char))) == NULL) {
        BadMem("Can't init mem for RegexMatch:APseq", 1);
    }
    SP = sequence;  /* SP now can point to the same array as sequence */
    D = DD->Dat[0];

    /* this loop needs to be repeated for ALL the regex's derived from the commandline or
    	from a file of regex's.  Those regex's will all be stored in RE in the std form and this
    	loop runs thru them all in order.  UNLIKE the fast hashing that gets done in Cutting(),
    	this will be a simple iterative process that searches for 1 regex thru the entire
    	sequence, then repeats it for the next one, until the end of the RE entries.  It SHOULD
    	be possible since they're stored as std RE entries to do things like do fragments or
    	multiple digests.  Since the regex's will overwhelmingly NOT be palindromes, I'll just do
    	the 'reverse' of the regex by default (and all of which will indicate in DD->Dat that
    	they are in the reverse strand..unless the regex DOES turn out to be a pal (tt[at]aa, as
    	a simple example), but if users were looking for those kind of simple patterns, they
    	should look for something else.  */

    /* for the 1st pass, I'm going to fake the inverse search by going thru the RE entries in
    	the forward direction and then Anti_Par()'ing the sequence, then searching again thru it.
    	This approach, while perhaps not as elegant as going thru the same sequence with
    	'reversed' regex's permits me to skip writing a regex parser and inverter.  I can spend
    	the time writing a DNA IUPAC->regex converter that will convert gcnytvmtk ->
    	gc.[ct]t[^t][ca]t[gt] */

//..see pcretest.c for the code to time the study routines - no need to do this yet

    base = 0; /* lagging counter to make sure that hits get counted correctly */
    strand = 1;


    for (ori=0; ori<2; ori++) { /* outer loop that Anti_Par()'s the sequence itself */
        /* note that on pass 2 the coordinates have to be inverted to match the forward numbers */
        // for each regex that might be in the queue
//      fprintf(stderr,"\n\n\t\t\tseq orientation = %d\n", ori);
        for (RC=F.RE_ST; RC < NumREs; RC++) {   /* start at RE_ST, not 1, not 0 */
			// fprintf(stderr,"RC= %d, Regex Name = %s, regex = %s\n",RC, RE[RC].E_nam, RE[RC].E_wsit[0]);
            /* get data required from RE so it can feed these loops */
            /* and need to generate the compiled_regex for each pattern */

            compiled_regex = pcre_compile(
                                 RE[RC].E_wsit[0],     /* the pattern */
                                 study_options,        /* default options */
                                 &error,               /* for error message */
                                 &erroffset,           /* for error offset */
                                 NULL);                /* use default character tables */

            if (compiled_regex == NULL) {
                fprintf(stderr, "Error # %s.\n", error);
                exit(1);
            }
            sptr = DD->BOS; /* corrected for BASE_OVERLAP, topology */
//         fprintf(stderr, "\nDD->BOS: %ld", DD->BOS);
//		 fprintf(stderr, "\nDD->EOS: %ld", DD->EOS);
//		 fprintf(stderr, "\nSEQlen: %ld\n", SEQlen );
            RxCode = 0;  /* zero those that need it (RxCode becomes nonzero at the end of the seq) */

			while (sptr < (DD->EOS) && rc != PCRE_ERROR_NOMATCH) {
//				fprintf(stderr, "sptr: %ld \n", sptr);
                rc = pcre_exec(
                         compiled_regex,   /* the compiled pattern */
                         NULL,             /* no extra data - we didn't study the pattern */
                         SP,               /* the subject string */
						 DD->EOS,           /* the adjusted length of the subject */
/*                         SEQlen,            the length of the subject */
                         sptr,             /* starting offset in the subject */
                         options,          /* options */
                         ovector,          /* output vector for substring information */
                         OVECCOUNT);       /* number of elements in the output vector */
//            if (rc >= 0) { printf("\nMatch succeeded at offset %d\n", ovector[0]); }
                if (rc > 0) { /* if there's a hit of at least 1 substring...*/
                    /* need to start stuffing the hit data in DD->Dat in the std way */
                    DD_Mem_Check();   /* checks DD->Dat and realloc's more mem if it needs it */
                    if (ori == 0) {
                        /* top strand get the index */
                        DD->Dat[(DD->Dat[0])++] = ovector[0] - BASE_OVERLAP; //- BASE_OVERLAP;
                    } else {/* bottom strand has to be inverted */
                        DD->Dat[(DD->Dat[0])++] = (SEQlen + BASE_OVERLAP + 1 - ovector[0]);
                    }
                    /* enter the RE index, negating if the bottom strand */
                    DD->Dat[(DD->Dat[0])++] = RC * strand;
                    RE[RC].E_Ncuts++;
                    hits++;
                }
//            fprintf(stderr, "\nat end of loop: ovector[0] = %d, sptr = %d ", ovector[0], sptr);
                sptr = ovector[0] + 1;
            }



//         fprintf(stderr, "Finished with pattern %d\n", RC);
            rc = 0;
        }
        if (ori == 0) {
            /* reverse the sequence so it can be searched again */
            Anti_Par(sequence, APseq, ovlapped_seq_len);
            sptr = DD->BOS;
            strand = -1;
            rc = 0;  /* zero those that need it (RxCode becomes nonzero at the end of the seq) */
            SP = APseq;  /* point SP at the Anti_Par()'ed sequence, and repeat... */
//            fprintf(stderr, "sequence: %s \n\n", sequence);
//            fprintf(stderr, "APseq: %s \n", APseq);
        }
    }
//    free(APseq);
    DD->Dat[DD->Dat[0]] = -22222; /* end it politely */
    return hits;
}


/* DNA_IUPACtoRegex() takes a pseudo regex expression in DNA IUPAC codes and
   translates it into REAL regex form for exapmple:
   gcnytvmtk -> gc.[ct]t[^t][ca]t[gt]
	the main problem seems to be getting the expression to this function without being creamed
	by the shell or other intermediate routines that are trying rabidly to interpret the
	characters as meaning something.  This is what they're paid to do, but it doesn't make
	things easy for the user or programmer.  This fn() assumes that the entering DNA string has
	been **downcased and pre-filtered** to remove any NON-IUAC characters from the string -
	maybe a faulty assumption.  */

void DNA_IUPACtoRegex(char *Regex, char *IUPAC) {
    /* this expects to have the incoming IUPAC string NULL terminated */
    int IP = 0, RP = 0;
    while (IUPAC[IP] != '\0') {
        switch (IUPAC[IP]) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
            Regex[RP++] = IUPAC[IP++];
            break;
        case 'r':
            sprintf(Regex+RP, "[ag]");
            RP += 4;
            IP++;
            break;
        case 'y':
            sprintf(Regex+RP, "[ct]");
            RP += 4;
            IP++;
            break;
        case 'w':
            sprintf(Regex+RP, "[at]");
            RP += 4;
            IP++;
            break;
        case 's':
            sprintf(Regex+RP, "[gc]");
            RP += 4;
            IP++;
            break;
        case 'm':
            sprintf(Regex+RP, "[ca]");
            RP += 4;
            IP++;
            break;
        case 'k':
            sprintf(Regex+RP, "[gt]");
            RP += 4;
            IP++;
            break;
        case 'b':
            sprintf(Regex+RP, "[^a]");
            RP += 4;
            IP++;
            break;
        case 'd':
            sprintf(Regex+RP, "[^c]");
            RP += 4;
            IP++;
            break;
        case 'h':
            sprintf(Regex+RP, "[^g]");
            RP += 4;
            IP++;
            break;
        case 'v':
            sprintf(Regex+RP, "[^t]");
            RP += 4;
            IP++;
            break;
        case 'n':
            sprintf(Regex+RP, ".");
            RP++;
            IP++;
            break;
            /* these char's need to be escaped to work with the POSIX regex, but not with
               the pcre routines, Now that tacg is using pcre, it won't be used anymore,
               but I'll leave them in for anyone who wants to re-use it
               NB:CANNOT escape '[' and ']' or they will try to be found as  [ ] literals   */
//         case '{': case '}': case '(': case ')': case '|':
//         case '?': case '*': case '+':
//             Regex[RP++] = '\\';  Regex[RP++] = IUPAC[IP++];
            break;
            /* or should the default be used to filter non-IUPACs like 'efj' ? */
        default:  /* everything else ><[]{}()+*?^\/,. etc just copies over */
            Regex[RP++] = IUPAC[IP++];
            break;
        }
    }
    Regex[RP] = '\0';
    if (F.Verbose > 0) {
        fprintf(stderr, "Original String is: %s; Converted String is: %s.\n", IUPAC, Regex);
    }
}



/* This function takes some basic argumnets and does essentially what PrintGelLadderMap() does -
   or at least the ladder part of it, but in postscript, so it can generate more info per inch
   and also more accurately.  It should take only about 1/8 - 3/16 to get the entire line for
   an RE onto the page.

   -takes hints as to width, etc from Flags struct.
 */
/* long *Degen_Log, int NumProtos(?) , float FrACGT[4] - needed to be passed in if going to mark degens as in plasmid code */
void tacg_psladder(int *Protos, char *Description,  char datestamp[80]) {
    int need_prolog=0, iBasesPerDiv, Cnt, iNumOfDivs, i, Pi;

    long flabels[15], seq_len = F.Seq_Len, itmp, li;

    double IN_FOR_LABEL, IN_FOR_SEQ, IN_FOR_NUMCUTS, pt_for_label, pt_for_seq, pt_for_numcuts,
    x1, y1, x2, y2, rBasesPerPt, rtmp, ftics[15],  rXPt, rYPt, cur_XPt,
    cur_YPt, bot_YPt, rPtPerMajDiv, rNumOfDivs, rPtsPerDiv, rBasesPerDiv ;

    char  *MapFileName, *tmpstr;

    FILE *PS;     /* output postscript file */
    struct stat *Fbuf=0x0;

    /* this is for US letter page initially */
    IN_FOR_LABEL   = 1.0;
    IN_FOR_SEQ     = 6.0;
    IN_FOR_NUMCUTS = 0.5;

    pt_for_label   = IN_FOR_LABEL * 72;
    pt_for_seq     = IN_FOR_SEQ * 72;
    pt_for_numcuts = IN_FOR_NUMCUTS * 72;


    if ((MapFileName = calloc(128, sizeof(char))) == NULL) {
        BadMem("RecentFuncs:MapFileName", 1);
    }
    /* mangle the file name to be compatible with a cgi call if it was so called */
    if (F.HTML != -1) { /* then we need the /tmp-named files */
        snprintf(MapFileName, 128, "%s/tacg_Map.ps", F.tmppath);
    } else {
        strcpy(MapFileName, "tacg_Map.ps\0"); /* name by itself to write into the $PWD */
        /* MapFileName = "tacg_Map.ps";   */
    }
    /*
       X Pts printing origin =  25.0;  ~1/3" up from the margin to start
       X Pts right border    = 540.0;  ~8" right margin
       y Pts printing top    = 720.0;  ~10.5" top margin
     */

    /* if file doesn't exist, need to write in prolog */
    if (stat(MapFileName, Fbuf) == -1) {
        need_prolog = 1;
    }

    if ((PS=fopen(MapFileName,"a")) == NULL) {  /* or the standard/optional one from SetFlags() */
        fprintf(stderr,"Cannot open the output postscript file \"%s\" for writing!!\n", MapFileName);
        exit(EXIT_FAILURE); /* print an error and die gracefully */
    }

    if (need_prolog == 1) {
        /* use the Illustrator PS prolog as it has a number of useful defs, so copy PROLOG into PS */
        /* rstr = PS_Prolog(); */
        tmpstr = tacg_PS_Prolog();
        fprintf(PS, "%s", tmpstr);
    }

    /* calc the scale and tic placements for the sequence */
    rBasesPerDiv = (float) seq_len / PREFERED_NUM_DIVS;
    iBasesPerDiv = (int) rBasesPerDiv;
    rBasesPerPt = rBasesPerDiv / pt_for_seq;
    rPtPerMajDiv = pt_for_seq / PREFERED_NUM_DIVS;

    rtmp = rBasesPerDiv;
    for (Cnt=1; rtmp >= 10.0; Cnt++) rtmp = rtmp /10; /* 1st part of finding a suitable div size */
    iBasesPerDiv = (int) (rtmp);
    for (i=1; i<Cnt;i++) iBasesPerDiv *= 10; /* and once found, multiply it back up to the right size */

    rNumOfDivs = (float)seq_len / (float)iBasesPerDiv; /* real number of divs to cover whole sequence */
    iNumOfDivs = (int) rNumOfDivs;  /* truncate to get the whole seq covered in the available space */

    /* want #pts per rNumOfDivs */
    rPtsPerDiv = pt_for_seq / rNumOfDivs;
    rBasesPerPt = (float) seq_len / pt_for_seq;
    iBasesPerDiv =  (int) seq_len / iNumOfDivs;

    /* set ftics to 0, then fill out ftics[] */
    ftics[0] = pt_for_label; /* 1st el is offset for the label space */
    flabels[0] = 0;
    for (i=1; i<15; i++) {
        ftics[i] = 0;    /* rest set to 0 */
        flabels[i] = 0;
    }
    rtmp = 0;
    i = 1, itmp = 0;
    /*       while (rtmp - pt_for_seq < DBL_EPSILON && i < 15) {  */
    /* fill out both ftics[] with incr'ing XPt values up to 10 or so and
       flabels[] up to the corresponding BASE values  */
    while (rtmp < pt_for_seq  && i < 15) {
        rtmp += rPtsPerDiv;
        itmp += iBasesPerDiv ;
        fprintf(stderr, "at i = %d, rtmp =  %f, itmp = %ld\n", i, rtmp, itmp);
        ftics[i] = rtmp;
        flabels[i++] = itmp;
    }

    /* Header for the page */
    tacg_psladder_header(PS, flabels, ftics, pt_for_seq, Description, datestamp);

    /* (FILE *PS, long flabels[15], double ftics[15], double pt_for_seq, char *Description, char datestamp[80])  */
    /* calc the starting point for the labeling &  set some PS basics like font info and stuff  */
    rXPt = 20; /* 1" from left side */
    rYPt = 520; /* 10" from bottom */
    fprintf(PS, "\n\n/Helvetica-Bold findfont 12 scalefont setfont\n");
    fprintf(PS, "newpath %f %f moveto\n 100 (Pattern Name) rightshow\n", rXPt, rYPt);
    fprintf(PS, "%f %f moveto\n 100 (Pattern Sequence) rightshow\n", rXPt, rYPt-15);
    /* still missing a separator line, but now write the seq line across */
    rYPt -= 7;
    rXPt += pt_for_seq;
    fprintf(PS, "%f %f moveto\n2 setlinewidth %f %f lineto stroke show\n",
            rXPt, rYPt, rXPt, rYPt);
    rtmp =  rYPt;
    rXPt += 36;
    fprintf(PS, "%f %f moveto\n(Number of Hits) show\n", rXPt, rYPt);

    /* now drop in the labels, tics, vertical hairlines */
    for (i=0; ftics[i] !=0; i++) {
        /* the numeric labels */
        fprintf(PS, "newpath %f  %f moveto\n (%ld) show\n",
                ftics[i], rYPt+15, flabels[i]) ;
        /* the tics */
        fprintf(PS, "newpath %f %f moveto\n1 setlinewidth %f %f lineto stroke show\n",
                ftics[i], rYPt+10,  ftics[i], rYPt-10);
        /* the vertical hairlines */
        fprintf(PS, "newpath %f  %f moveto\n0.1 setlinewidth %f  0 lineto stroke show\n",
                ftics[i], rYPt,  ftics[i]);
    }

    /* set up the limits */
    cur_XPt = 72.0;
    cur_YPt = 648.0;
    bot_YPt = 36.0;

    /* go to starting point & translate the ladder format to PS */
    fprintf(PS, "newpath %f  %f  moveto \%\% move to the ladder start point\n", cur_XPt, cur_YPt);

    /* and now do the repeating rungs, checking that the whole next rung
       will fit in the remaining space (est ~1/4" per rung) */

    Pi = 0; /* Pi is Protos index */

    while (Protos[Pi] != 0) { /* 1st entry past the good ones should be a 0  */
        if  (RE_Filters_OK(Protos[Pi]) == 1) { /* if the RE passes all the filters */
            /* write it out (is Protos re-sorted to honor the -c flag?) */
            /* if there's not enuf space, start a new page, write out a new header & continue*/
            if (cur_YPt < (bot_YPt + 18)) {
                fprintf(PS, "showpage %% new page\n");
                tacg_psladder_header(PS, flabels, ftics, pt_for_seq, Description, datestamp);
            }
            /* OK, write out the rung - 1st the labels  */
            cur_XPt -= 7; /* move the labels a bit to the left. */
            fprintf(PS, "newpath %f %f moveto \n(%s) rightshow \n", cur_XPt, cur_YPt, RE[Protos[Pi]].E_nam);
            cur_YPt -= 10; /* move down a bit more to write the site under the name */
            fprintf(PS, "newpath %f %f moveto \n(%s) rightshow \n", cur_XPt, cur_YPt, RE[Protos[Pi]].E_raw_sit);
            fprintf(PS, "newpath %f %f moveto \n(%d) rightshow \n",
                    pt_for_label + pt_for_seq + 7, cur_YPt + 10, RE[Protos[Pi]].E_Ncuts);

            /* draw the sequence line from the start of the label to the end and endpoints */
            cur_XPt = 20; /* this is a static point */
            cur_YPt += 10;
            fprintf(PS, "newpath %f %f moveto \n0.3 setlinewidth %f %f lineto stroke \n",
                    cur_XPt, cur_YPt, cur_XPt+pt_for_seq, cur_YPt);
            fprintf(PS, "newpath %f %f moveto \n0.1 setlinewidth %f %f lineto stroke \n",
                    pt_for_label, cur_YPt+10, pt_for_label, cur_YPt-10);
            cur_XPt = pt_for_label+pt_for_seq;
            fprintf(PS, "newpath %f %f moveto \n0.1 setlinewidth %f %f lineto stroke show\n",
                    cur_XPt, cur_YPt+10, cur_XPt, cur_YPt-10);


            /* now the tics; need to make it into a fn().  cur_YPt is still set to the right place */

            for (li=1;li<RE[Protos[Pi]].Sites[0];li++) { /* the 0th el of Sites is the end */
                /* calc the Xoffset in Pts - all of the following could be fn() but to what end?*/
                cur_XPt = pt_for_label + (pt_for_seq * (float)RE[Protos[Pi]].Sites[li]/(float)seq_len);
                y1 = cur_YPt+3;
                y2 = cur_YPt-3;
                if (RE[Protos[Pi]].E_olap > 0)      {
                    x1 = cur_XPt-2;
                    x2 = cur_XPt+2;
                } else if (RE[Protos[Pi]].E_olap < 0) {
                    x1 = cur_XPt+2;
                    x2 = cur_XPt-2;
                } else {
                    x1 = x2 = cur_XPt;
                }
                fprintf(PS, "newpath %f %f moveto 0.1 setlinewidth %f %f lineto stroke show\n",
                        x1, y1,                        x2, y2 );
            }
        } else {
            fprintf(stderr, "\nRE # %d filtered out", Pi);
        }
        Pi++;
    }
    free(MapFileName);
    /* free(tmpstr);  */
}

/* tacg_psladder_header() writes the header for the psladder page */
void tacg_psladder_header(FILE *PS, long flabels[15], double ftics[15], double pt_for_seq, char *Description, char datestamp[80]) {

//   int i;
    double rXPt, rYPt, x_left, y_bot, x_right, y_top, yoffset, xoffset;  //rtmp;

    /* calc the starting point for the labeling &  set some PS basics like font info and stuff  */
    rXPt = 72; /* 1" from left side */
    rYPt = 720; /* 10" from bottom */

    x_left = y_bot = 25; /* ~1/3" up from the margin to start */
    x_right = 500; /* ~8" right margin */
    y_top = 720; /* ~10.5" top margin */

    /* ------------------------  write the title, date, etc --------------------------------------- */
    /* write a header at top left */
    yoffset = y_top;
    fprintf(PS, "\n\n/Helvetica-Bold findfont 10 scalefont setfont\n");
    fprintf(PS, "newpath %f %f moveto\n (Description: %s) show      \n",
            x_left, yoffset, Description);
    yoffset -= 12;
    fprintf(PS, "%f %f moveto\n (Date: %s) show\n", x_left, yoffset, datestamp );

    fprintf(PS, "\n\n/Helvetica findfont 8 scalefont setfont\n");
    yoffset -= 10;
    fprintf(PS, "newpath %f %f moveto\n (Map by tacg v. %s (http://tacg.sf.net)) show\n",
            x_left, yoffset, TACG_VERSION);
    fprintf(PS, "\n\n/Helvetica findfont 8 scalefont setfont\n");
    yoffset -= 10;
    fprintf(PS, "newpath %f %f moveto\n (Suggestions to Harry Mangalam (hjm@tacgi.com)) show\n",
            x_left, yoffset);

    /* print a little legend at top right */
    yoffset = 0;
    xoffset = x_right-150;
    fprintf(PS, "\n/Helvetica findfont 7 scalefont setfont\n");
    fprintf(PS, "newpath %f %f moveto\n (Label = Name[base#]overhang (*), where: ) show\n", xoffset, (y_top-yoffset));
    yoffset += 8;
    fprintf(PS, " %f %f moveto\n (Name = Name of the RE or pattern,  ) show\n", xoffset, (y_top-yoffset));
    yoffset += 8;
    fprintf(PS, " %f %f moveto\n ([base#] = cut site of the RE or ~center of the pattern) show\n",  xoffset, (y_top-yoffset));
    yoffset += 8;
    fprintf(PS, " %f %f moveto\n (overhang: \\\\ = 5', / = 3', | = blunt. '*' = 1 hit only.) show\n",  xoffset, (y_top-yoffset));

}
