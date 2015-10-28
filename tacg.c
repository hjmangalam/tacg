/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright  1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847)  */

/* $Id: tacg.c,v 1.4 2005/01/26 22:26:13 mangalam Exp $  */

/* The use of this software (except that by James Knight, which is described in
   'seqio.c', and by Gray Watson, which is described in strsep.c) is bound by
   the notice that appears in the file 'tacg.h' which should accompany this
   file.  In the event that 'tacg.h' is not bundled with this file, please
   contact the author.
*/

#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <errno.h>
#ifdef IRIX
extern char *sys_errlist[];
#else
extern const char *const sys_errlist[];
#endif

#include "seqio.h"
#include "tacg.h"   /* contains all the defines, includes, function prototypes for both main() and functions;
simple enuf not to have a version #; I'll just keep it up to date as I go...*/

/* Declare the key GLOBALS */
struct RE_struct *RE;
struct flag_struct F;
struct ORF_struct *ORFs[6];
struct Digest_Data *DD;
struct Prox_struct *PP;
struct PW_struct *PosWin;
struct Clone_struct *Clone_S;  /* need to zero it still */
struct Rule_struct *Rules; /* like RE, this has to be reallocatable */


int main(int argc, char *argv[])
{

/* Declarations */

int i, ii, j, d, k, b, e, mm=0, m,  thisframe=0, moo=0, *hashtable[4096], NumREs, codon_Table=0,
NumProtos=0, Cur_RE, gel=0, HTML=0, max_okline, rule=0,  sys_err = 0, fr,
ProtoCutters, *Protos, topo,  basesPerLine, newsize,
*metstop[6], //array for logging the mets & stops for each frame
ORFstep=100,  // size of each increment in alloc size
//         ORF_End=100, nitial 'end-of-alloced-memory' for metstop
msp, aap, reps, AI=0, SEQlen=0, SEQNum=0, NumPosSlidWin=0, RE_dam=0, RE_dcm=0;
//VERBOSE,


short MoreSequence;

long  /* DEL *BackCheck,  BackCheck[] checks proto entries in the Linear map output to remove duplicates */
l, Base_cutsite,  seq_len = 0, tot_seq_Cnt = 1, Xn=0, Buf_Seq_Len,
Pad_Seq_Len, Bare_Seq_Len, /* ProtoNameHash[RE_NAME_TAB_SIZ],  */
*Degen_Log; /* Degen_Log holds the degen locations.  both singles and multiples are noted by start/end positions
[0] = currently calloc'ed end
[1] = current end of data
      [2] = start of degen and all starts
            [3] = end of degens  and all ends, etc   */
                  /* primes[26] is for setting sizes of realhash tables: from
                  http://planetmath.org/encyclopedia/GoodHashTablePrimes.html */
// long primes[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157,
// 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
// 50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};

                  float FrACGT[4];

char  *sequence=0x0,  *SEQseq=0x0, *TmpSeq, *Prot4ORF, *s, *EnzFileName, *inputstream,
       *RTDNA, *rulestring, datestamp[80], *gs_bin, gs_cmd[200];

/* the full set of environment variables that might need to be grokked from main() */
/*   static char *Env_Vars[] = {"TACG_GS_BIN", "$TACG_GS (ghostscript binary)",
                              "PWD",         "$PWD (Current Directory)",
                              "HOME",        "$HOME (Home Directory)",
                              "TACGLIB",     "$TACGLIB (TACG library)"};
*/

char Codons[MAX_ORGS][N_CODONS][7];
char *Codon_Labels[MAX_ORGS];
struct SE_struct *SelEnz;
//    struct tm *date;
struct stat *Fbuf;
time_t now;

SEQFILE *SEQfp=0x0;
/* SEQfp is a pointer to a regular file, in this case, it will be
   repointed to stdin to handle input the way God intended.
   SEQDBfp will be used to handle pointing to entire DATABASES ..later */
/* SEQFILE *SEQDBfp;    SEQDBfp is a pointer to a DATABASE,
   not a regular file or even a file of sequences - they
   can be handled by the normal run of affairs.  It requires
   that the user sets up the BIOSEQ environment as described in
   the SEQIO docs.  This may use a different struct than SEQFILE, tho */
SEQINFO *SEQinfo;
SEQinfo = 0x0;
FILE *tmpfp;

tmpfp = tmpfile();
if (tmpfp == NULL) {
    fprintf(stderr,"!!tempfile creation failed - que??\n");
    exit (0);
}
fprintf(tmpfp, "Temporary file created by tacg.\n..\n");

/* Initialize whatever vars are needed */
//SEQformat = " ";
tacg_SetF();

F.PID = getpid(); /* get Process ID of run for various self-ID reasons */

/* structs needed for datestamping */
time(&now);                           /* Get the current calendar time.   */
strftime(datestamp, 80, "%c", localtime(&now));  /* Convert that to a string.  */

/* Probably don't need any of this for the nmer stuff */
/* using RE as a realloc-able pointer */

if ((RE = calloc(INIT_RES, sizeof(*RE))) == NULL) {
    BadMem("Can't init mem for RE", 1);
}

if ((Fbuf = calloc(1, sizeof(*Fbuf))) == NULL) {
    BadMem("tacg:Fbuf",1);
}

/* And doing the same thing for SelEnz (if use it as a realloc-able) */
if ((SelEnz = calloc(MAX_SEL_ENZ, sizeof(*SelEnz)+1)) == NULL) {
    BadMem("Can't init mem for SelEnz", 1);
}

/* and give some mem to rulestring (MAX_LEN_RULESTRING=2000 at this time) */
if ((rulestring = calloc(MAX_LEN_RULESTRING, sizeof(char))) == NULL) {
    BadMem("Can't init mem for rulestring", 1);
}

// init *Degen_Log to 1000 - enough for 500 degens (0=end of alloc'ed mem,
// 1=current end, 2=start, 3=end, etc
if ((Degen_Log = calloc(1000, sizeof(long))) == NULL) BadMem("Can't init mem for Degen_log", 1);
Degen_Log[0] = 1000;
Degen_Log[1] = 2;  /* set up the initial counters */

/* since DD is pointer to struct, have to explicitly get mem for it */
if ((DD = calloc(1, sizeof(*DD))) == NULL) BadMem("Can't init mem for DD", 1);
DD->Dat = (long *) calloc(10000, sizeof(long)); /* init the pointer to 10,000 */
if (DD->Dat== NULL) BadMem("calloc failed -> DD->Dat mem!", 1);
DD->Dat[0] = 2;      /* this is Cnt and it has to start at 2 to skip these two unholy counters */
DD->Dat[1] = 10000;  /* this is Siz */

s = (char *) calloc(256, sizeof(char));
if (s == NULL) BadMem("failed to calloc mem for 's'", 0);

/* Figure out what the program should do, by parsing the commandline options */
/* if the program is invoked by name alone, tell the user how to use it */
if (argc < 2 ) {
    fprintf(stderr, "type 'tacg -h' or 'man tacg' for more help on program use or type: \n"
            "tacg -n6 -slLc -S -F2 < [your.input.file]\nfor an example of what tacg can do.\n");
    exit(0);
}
/* Now get and process the flags and return the alt REBASE file, if any */
EnzFileName = SetFlags(argc, argv, SelEnz, tmpfp, rulestring, datestamp);

if (F.Pat == 1 && F.Output == 0) {
    F.Output =1;  // in case they used -p without -S
}
if (F.Output == 0) { /* if it's still 0, then no output flags have been set, ergo no output will be generated */
    if (F.Output == 0)
        fprintf(stderr, "You haven't requested any output.  You need to explicitly ask for it.\n"
                "Use any of: -F -g -G -l -L -O -P -s -S -X --rule --rulefile. \nType 'tacg -h' for a brief usage.\n");
    exit(1);
}

if (F.Timeops) {
    tacg_TimeOps();
    exit(0);
}

/* NumREs is set to 1 as the default, use negation to mark doppels, changes dep'g on where dam, dcm */
NumREs = 1;
F.RE_ST = 1; /* F.RE_ST is the global flag that tells the various loops where to start */

/* if F.dam is set, need to set it up in RE */
if (F.dam == 1) {
    RE_dam = NumREs;  /* set the RE offset to use */
    F.RE_ST++;  /* incr NumREs, F.RE_ST  */
    Add_dam(NumREs++); /* only add it if we need it */
}

/* if F.dcm is set, need to set it up in RE */
if (F.dcm == 1) {   /* only add it if we need it */
    RE_dcm = NumREs;  /* set the RE offset to use */
    F.RE_ST++;  /* incr NumREs, F.RE_ST  */
    Add_dcm(NumREs++); /* only add it if we need it */
}

/* Get and process the Enzymes from the rebase or alternate file or regex pattern*/
if (F.Matrix == 0 && F.Regex == 0) { /* want the standard rebase type file */
    NumProtos = ReadEnz(EnzFileName, tmpfp, &NumREs, SelEnz, hashtable, seq_len);
} else if (F.Matrix > 0) { /* then we must want a matrix match */
    NumProtos = ReadMatrix(EnzFileName, &NumREs, SelEnz);
} else if (PCRE == 1) { /* there must be a regex file waiting to be read */
    NumProtos = ReadRegex(EnzFileName, tmpfp, &NumREs, SelEnz);
} else {
    fprintf(stderr, "tacg compiled without regex support - install pcre and re-compile\n");
    exit(1);
}

if (F.Verbose > 0) fprintf(stderr, "\n# of Prototype patterns = %d; Total # RE Entries = %d\n", NumProtos, NumREs);

/* Now calloc the Protos[] array using the NumProtos; used to crosscheck Sites output and to sort/re-sort RE entries */
Protos = (int *)calloc((size_t)NumProtos+2, sizeof(int));
if (Protos == NULL) BadMem("Can't init mem for Protos", 1);

/* Get the Codon Usage Table from the (only slightly) modified NCBI Codon Table file*/
Read_NCBI_Codon_Data(Codons, Codon_Labels);


/*  now calc the max length of the names of all the patterns stored in RE */
F.MaxPatNameLen = 0;
for (i=0; i<=NumREs; i++) {
    j = RE[i].E_nam_l;
    if (F.MaxPatNameLen < j) {
        F.MaxPatNameLen = j;
//        index = i;
    }
}
if (F.MaxPatNameLen < O_LMARG) {
    F.MaxPatNameLen = O_LMARG;
}

/* ***************************************************************************************** */
/* HERE IT IS - THE ENTRY POINT FOR SEQIO - IF WANTED, SEQIO GETS THE SEQUENCE ETC, AND
HANDS IT OFF TO TACG; IF NOT (NO -d FLAG), WE DO IT THE REGULAR WAY. The patterns don't
change; none of the options change; the only thing that changes is that new sequences are
read in via SEQIO, they're padded end for end via GetSequence2(), and then they're submitted
for searching in the same way as the old version.  Have to check: if F.DB == 1, there's
a DATABASE to read in via SEQIO; if it's == 0, read the sequence the normal way from stdin.
This can be either a single file containing a single sequence (as before) or a file of
sequences, like a FASTA-formatted file:

 >comments for sequence 1
 gcatcgatcgcatagcac...
 >comments for sequence 2
 gcatcgagcatacgcgca...
 >comments for sequence 3
 atgggtcagggctcttcac...

So.... the whole thing from immediately below to almost the VERY end, will be in a 'while' loop that
keeps going as long as there's sequence coming in... and CAN document the output stanzas with the 'info' bits
that get sucked up by SEQIO as it goes along...

For 1st pass, handle only the FILE-oriented input; then add the DATABASE-handling
routines afterwards; may be able to use the seqfopen2() later to do either file
OR database reading with the same fn.
*/
/* ******************************************************************************************* */
//if (F.infile == NULL || F.infile == "")  //then we get seq on stdin, else, we use the filename
if (F.infile == NULL) { //then we get seq on stdin, else, we use the filename
    //      if ((inputstream = (char *) calloc((size_t)5, sizeof(char))) == NULL) { BadMem("inputstream", 1); } // this shouldn't require a calloc to get mem - just assign it as below
    inputstream = "-";
} else { //already verified that it exists in SetFlags, so we just try to read it
    inputstream = F.infile; //just point inputstream at F.infile, no calloc needed
}
MoreSequence = 0;
if (F.RawSeq ==1 ) {
    sequence = GetSequence(&tot_seq_Cnt, &seq_len, &Degen_Log);
    SEQNum = 0;
    MoreSequence = 1;
} else if ((SEQfp = seqfopen(inputstream, "r", NULL)) == NULL) { /* "-", "r" = read from stdin  */
    fprintf(stderr, "\nWhere's the input of this thang?!?  Either the file or input stream was empty,\n"
            "it was in a format that seqio did't grok (was there an error message\n"
            "indicating that just above?), or the input may have been damaged by editing \n"
            "or transfer so that it doesn't appear like a known sequence format anymore.\n"
            "Any of these sound applicable?  You might also try the '--raw' option to \n"
            "read raw, unformatted, or badly formatted sequence (assumes all IUPAC \n"
            "characters are valid sequence but will ignore numbers, newlines, tabs, etc.\n");
    exit(1);
    // nmer - can use following code to loop thru sequence to load char * array for nmer analysis
} else if ((F.RawSeq != 1) && (SEQseq = seqfgetseq(SEQfp, &SEQlen, 0)) != NULL) {
    /* get the sequence via SEQIO and munge it into the standard format with PadEndForEnd(),
       pass it thru the rest of tacg, and then come back for more until there's no more sequences */
    /* entry can be either thru seqio or via raw sequence input  */
    SEQNum = 0;
    SEQinfo = seqfinfo(SEQfp, 0);         /* this one doesn't cause a sequence increment */
    /* to protect on some OSs that explode when printing a NULL, catch them here
       SHOULD actually be fixed in the seqio code! */
//         printf("size of SEQinfo: %d\n", sizeof(SEQinfo));
    if (SEQinfo != NULL) {
        if (SEQinfo->format == NULL)      SEQinfo->format      = "Unidentified or Raw";
        if (SEQinfo->description == NULL) SEQinfo->description = "None";
        if (SEQinfo->idlist == NULL)      SEQinfo->idlist      = "None";
        if (SEQinfo->dbname == NULL)      SEQinfo->dbname      = "N/A";
    }
    seq_len = (long)SEQlen;
    //hjm519 - sequence is the padded sequence but seq_len has NOT been adjusted. tot_seq_Cnt is the padded length
    sequence = GetSequence2(SEQseq, &tot_seq_Cnt, &seq_len, &Degen_Log);
    //free(SEQseq); /* as soon as it's moved into sequence, don't need the original anymore */
    //SEQseq=0x0;
    MoreSequence = 1;
//	printf("\nSizeof tot_seq_Cnt: %ld\n", tot_seq_Cnt);

    /* Check the Degen_Log to see if it has picked up the degens */
    if (F.Verbose > 0) {
        fprintf(stdout, "INFO:tacg:Degen_Log has gone to %ld\n", Degen_Log[1]);
        for (i=0; i<=Degen_Log[1]; i+=2) {
            fprintf(stdout, "INFO:tacg:Start[%d]=%ld, End[%d]= %ld, # of consecs = %ld\n", i, Degen_Log[i], i+1, Degen_Log[i+1],
                    (Degen_Log[i+1]-Degen_Log[i]+1));
        }
    }
}

/* now sequence should be OK, regardless of which way it came in, but it needs to be
   checked again at the END of tacg to either read in more sequence or fail to  */

Buf_Seq_Len = strlen(sequence);   /* seq_len + (2*BASE_OVERLAP); */

F.Seq_Len = seq_len; /* the actual sequence length (NOT padded) note it in F... to pass it around, if nec */

if (F.Rev == 1) {
    Reverse(sequence); /* works in-place, no length measure required */
} else if (F.Comp == 1) {
    Complement_In_Place(sequence, Buf_Seq_Len); /* works in-place; note that Rev_Comp() (erroneously) does the same thing, but not in-place */
} else if (F.RevComp == 1) {
    TmpSeq = calloc((size_t)Buf_Seq_Len+1, sizeof(char)); /* dmalloc: need the '+1' to avoid smashing pickets */
    if (TmpSeq == NULL) BadMem("tacg.c:~308 TmpSeq calloc err", 1);
    Anti_Par(sequence, TmpSeq, Buf_Seq_Len);
    free(sequence);
    sequence = TmpSeq;
}

while (MoreSequence) {
    if (F.RawSeq == 1) MoreSequence = 0;
    SEQNum++;
    if (F.Seq_Len < (F.XtEnd - F.XtBeg) && F.Verbose > 0) {
        //issue warning that size of extract is larger than entire sequence length
        fprintf(stderr, "\nWARN: Sequence: %s,\nextracted seq length (specified with -X option) is greater than entire\nsequence length.  Sequence  will be truncated\n", SEQinfo->description);
    }

    for (i=0; i<4; i++) {
        FrACGT[i] = (float)F.NBases[i]/(float)seq_len;
    }
    /* If want to get Estimated sites for each RE for each sequence, should do it here, after sequence is got -
       move that code OUT of ReadEnz() - could add this behind a -E (estimated sites) flag so that it would
       only be calc'ed on demand - for large pattern databases searching large DNA databases, this could take a
       a significant amount of time, no?  Well, for the 1st pass, just add it in and try it for a bad case scenario */
    if (F.Matrix == 0 && F.Regex == 0) { /* only for std RE mapping - otherwise it's unpredictable */
        for (i=F.RE_ST; i<NumREs; i++) {      /* calculate est # of hits based on ACGT distrib */
            /* below incr the PROTO RE entry EstSites with the EstSites of each of its derivatives */
            RE[RE[i].proto].EstSites += (float) (seq_len * CalcEstSites(RE[i].E_wsit[1], RE[i].E_len, FrACGT));
        }
    }
    /* start SEQIO rearrangement */
    /* set up more vars, based on flags - could use only array to pass vars, but this is more readable */
    topo = (int)F.Topo;
    basesPerLine = (int)F.Width;
    if (F.GelLO != 0) gel = 1;
    HTML = (int)F.HTML;
//    VERBOSE = (int)F.Verbose;

    /* Decide from flags how to handle translation options ... */
    if (F.XLinFr != 0) {
        Xn = F.XLinFr;
//        n_letters = (int)F.Lab1or3;
    }

    max_okline =  O_SEQ_LINE + 4 + (int)labs(Xn);   /* how many more lines of output are required for the translation */

    if (F.Strands == 1)  max_okline -= 1; /* if only 1 strand, max_okline is decr by 1 */
    /*   if (F.Tics == 0) max_okline -= 1; if don't want tics, max_okline is decr by 1 */

    codon_Table = (int)F.CodTab; /* which Codon Table to use in Translate() */

    if (topo == 0) {           /* if the sequence is circular */
        DD->BOS = 1;                /* Beginning Of Sequence starts at the very beginning of 'sequence' */
        DD->EOS = tot_seq_Cnt;      /* and End Of Sequence ends at very end  */
    } else {                  /* but if it's linear */
        DD->BOS = BASE_OVERLAP + 1; /* BOS starts after the beginning buffer */
        DD->EOS = tot_seq_Cnt - BASE_OVERLAP;   /* and EOS is before the ending buffer */
    }
    /* the '--silent' sites search goes here, as it has to Translate the DNA in,
       then reverse translate it to Degenerate DNA out, THEN and ONLY THEN, subject
       it to all the analyses down the pike */
    if (F.Silent == 1) {
        if ((Prot4ORF = (char *)calloc((seq_len/3 + 2), sizeof(char))) == NULL) BadMem("tacg:Prot4ORF\n",1);
        if ((SEQseq = (char *)calloc((seq_len+(2*BASE_OVERLAP)+2), sizeof(char))) == NULL) BadMem("tacg:SEQseq\n",1);

        Translate(sequence+DD->BOS, Prot4ORF, seq_len, 2, Codons, codon_Table);
        SEQseq = ReverseTranslate(Prot4ORF, Codons, codon_Table);
        RTDNA = PadEndForEnd(SEQseq, &Pad_Seq_Len, &Bare_Seq_Len);
        free(SEQseq);
        free(sequence);
        sequence = RTDNA;  /* but can't free(RTDNA) cus now sequence is pointing at it!! */
        free(Prot4ORF);
    }

    /************   Here's where the rubber hits road.. Go and cut/match the sequence   *****************/
    /* as tacg expands it repertoire, there are now a few functions that do not require the
       restriction analysis, so there should be a test to skip it if it's not necessary - but not yet ...*/

    /* Note that each of these will will require the insertion of the '--extract' code if they are each
       going to be able to extract sequence OR it will have to be done after the fact, using the
       data from DD->Dat (which would make a lot more sense frankly) I will have to be able to check
       what the orientation is tho, but this should not be too much of a problem */

    if (F.nmer > 0) { /* if we want the nmer analysis */
        /* what do we need to pass?? What do we return?  The number of degens?*/
        // SEQseq is the unpadded sequence direct from seqio - probably easier to handle
        fprintf(stdout, "nmer analysis done; exiting now\n");
        free(SEQseq); //if we're going to use it here - otherwise free it in other instances
        exit (1);
    }
    if (F.Matrix == 0 && F.Regex == 0) { /* if there isn't a MatrixMatch cutoff value stored, do regular Cutting() */
        Cutting(sequence, hashtable); /* this is all there is left in main() */
    } else if (F.Matrix > 0) {
        MatrixMatch(sequence, NumREs);  /* is this all it really needs??!? */
    } else if (F.Regex != 0) {
        if (PCRE == 1) {
            b = RegexMatch(seq_len, sequence, NumREs);  /* that was simple wasn't it? */
        } else {
            fprintf(stderr, "tacg compiled without regex support - install pcre and re-compile\n");
            exit(1);
        }
    }

    // translate the data from DD->Dat to Dig_Sits for calc'ing fragment sizes
    // calloc the space for *RE[].Sites and sister array *Frag_Siz[] and init
    //  the 0th el to point to 1st free spot
    // can also alloc the space for the bit array that will hold the orientation data
    for (i=1; i<NumREs; i++) { /* start at 1, NOT F.RE_ST, NOT 0; have to get space for RE[1,2].Sites */
        if (RE[i].proto == i) {
            /* calloc Sites */
            RE[i].Sites = (long *) calloc((size_t)RE[i].E_Ncuts+2, sizeof(long));
            if (RE[i].Sites == NULL) BadMem("!! calloc failed: RE[i].Sites\n", 1);
            RE[i].Sites[0] = 1;      /* the [i][0] el points to the next free element, so init here */
            /* calloc Frags  - not needed for dam/dcm, but it's not a huge mem waste...*/
            RE[i].Frags = (long *) calloc ((size_t)RE[i].E_Ncuts+2, sizeof(long)); /* need '+2' for extra frag */
            if (RE[i].Frags == NULL) BadMem("!!calloc failed: RE[i].Frags\n", 1);
        }
    }

    /* The following loads RE[].Sites, the array that tracks the digestion sites  */
    i=2; /* REMEMBER - DD->Dat real data starts at [2]; [0] is Cnt, [1] is Siz */
    while (DD->Dat[i] != -22222) { /* -22222 is the end point indicator */
        if (F.Verbose > 2) {
            fprintf(stderr,"\n");
            for (ii=0; ii<10; ii++) fprintf(stderr, "%6d", i+ii);
            fprintf(stderr, "  [i]\n");
            for (ii=0; ii<10; ii++) fprintf(stderr, "%6ld", DD->Dat[i+ii]);
            fprintf(stderr,"  DD->Dat[i]\n");
        }

        b = 1; /* this will transfer the sign of the orientation from the index of the pattern to the site */
        Base_cutsite = DD->Dat[i++]; /* i now points to RE# */
        if (F.Verbose > 1) fprintf(stderr, "Base_cutsite=%ld\n", Base_cutsite);
        e = i; /* saves a '-' */

        if (DD->Dat[i] < 0 ) b = -1;
        /* mm to shorten next line, incr's i to point to next position - remember, it's the
           PATTERN index that's negated in Cutting(), not the offset */
        mm = (int)labs(DD->Dat[i++]);
        if (F.Verbose > 2) {
            for (ii=0; ii<10; ii++) fprintf(stderr, "%6d", ii);
            fprintf(stderr, "  [i]\n");
            for (ii=0; ii<10; ii++) fprintf(stderr, "%6ld", RE[1].Sites[ii]);
            fprintf(stderr,"  RE[1].Sites[i]\n");
        }

        if (F.Verbose > 1) {
            fprintf(stderr, "1: RE[%d].Sites[0] = %ld, RE[%d].E_Ncuts=%d\n", mm, RE[mm].Sites[0], mm,RE[mm].E_Ncuts);
            fprintf(stderr, "1: %ld = %ld ?\n", RE[mm].Sites[RE[mm].Sites[0]-1], Base_cutsite);
        }
        if (Base_cutsite == 1 ||  ((RE[mm].Sites[RE[mm].Sites[0]-1] != Base_cutsite) && (Base_cutsite <= seq_len))) { /* if it's not the same as the one before (-P, degens) */
            RE[mm].Sites[RE[mm].Sites[0]++] = Base_cutsite * b;   /* weird, self-referential definition, but it works.. */
            /* 'b' above is the switch that makes the site neg if the pattern is in the rc */
        }
        /* adjustment for the above adjustment - this doesn't decr the count
        it just makes it (pointer to the next el - 1) */
        RE[mm].E_Ncuts = (int) RE[mm].Sites[0] -1; // -1
        if (F.Verbose > 1) fprintf(stderr, "2: RE[%d].Sites[0] = %ld, RE[%d].E_Ncuts=%d\n", mm, RE[mm].Sites[0], mm,RE[mm].E_Ncuts);
    }    /*   [RE#] [Pointer to next free el++ ]     */

    /* RE[].Sites now loaded with REAL WORLD Coords so calc'n of Frag_Siz should be accurate
       (except that some are negative) and this will have to be taken account of */

    /* Need to do this here so that we can sort and de-jitter the combined one as well as the
       individual ones */
    /* Now for NEW '-x' option, where we want ALL the specified enzymes to cut in the same digest,
       we need to have an RE struct appended to the rest of them, dedicated to this digest.
       This struct is automatically added in ReadEnz() fn - there will be a minimum of 12 extra
       structs (a bit wasteful, but they're left for this reason  */
    /*   The index to use is NumREs+1 to be sure of clearing the last one. RE[] now includes .Sites, .Orient
       and .Frags and so all the required data is now in RE[] */
    /* Now, if needed by -x flag, copy all the adjusted sites to the last RE struct, one after another
       so that we can calculate the fragments in the same stanza below...*/
    /* FUTURE: F.Xplicit reflects the number of iterations of '-x' that need to be calculated */
    /* On 1st pass, I'll do this for a SINGLE group of sites and then include the loopy stuff later - it has to
       loop over each group of REs, so that each one needs to have another RE[] struct alloc'ed to it */

    /* if need to handle pattern names passed in via -x flag and PERHAPS
       calculate the data for the ALL CUTS IN 1 DIGEST option ... */
    /* there has to be more than one pattern requested via -x flag, AND you have to WANT
       the multiple digest */

    /* this should ALSO be functionized to get it out of tacg.c */
    if (F.Xplicit == 1 && F.AllInc == 1) {
        /* NumREs has to be incr to point to the end of the RE struct  */
        AI = NumREs++;
        k = 1;
        strcpy(RE[AI].E_nam, "Combined");
        RE[AI].E_Ncuts = 0;
        /* for all the REs in the -x list (aka MulDigest[]), sum the # of all their sites to alloc mem */
        for (i=1; i<NumREs; i++) {
            if (RE[i].proto == i) RE[AI].E_Ncuts += RE[i].E_Ncuts;
        }
        /* and now ask for the AI .Sites, .Frags arrays to be alloc'ed to that size */
        RE[AI].Sites = (long *) calloc ((size_t)RE[AI].E_Ncuts + 2, sizeof(long)); /* '+ 2' for extra site */
        if (RE[AI].Sites == NULL) BadMem("!! calloc failed: RE[AI].Sites\n", 1);
        RE[AI].Frags = (long *) calloc ((size_t)RE[AI].E_Ncuts + 2, sizeof(long)); /* '+ 2' for extra frag */
        if (RE[AI].Frags == NULL) BadMem("!! calloc failed: RE[AI].Frags\n", 1);

        for (i=F.RE_ST; i<AI; i++) { /* for all the REs in the -x list, copy their sites to RE[AI].Sites */
            if (RE[i].proto == i) { /* '==' was '=' until 10-23-98 */
                for (j=1; j < (int) RE[i].Sites[0]; j++) {
                    RE[AI].Sites[k++] = RE[i].Sites[j];
                }
            }
        }

        RE[AI].Sites[0] = k;   /* just to make it like all the others [0] points to the next free space */
        RE[AI].E_Ncuts = k-1;  /* and adjust the number */
        RE[AI].proto = AI;     /* and make sure that it thinks it a proto */
        RE[AI].E_mag = 8;      /* and make sure the mag is enough to let it pass the default checks */
        RE[AI].E_pal = 1;      /* and make sure that it thinks it's a pal */
        RE[AI].E_raw_sit = "All Sites";   /* and give it a label so that it can output something */
        // Now sort the Sites so that they're in order (Check whether the # of cuts is correct or not)
        qsort(RE[AI].Sites+1, (size_t)RE[AI].E_Ncuts, sizeof(long), abscompare);
    } /* this should end with all the REs separate sites being appended in one long array */

    /* OK - try to put the dam/dcm stuff in here by zeroing those sites are are contextually sensitive to dam or dcm */
    /* ALSO!  if a contextually sensitive RE is also palindromic, AND the dam/dcm site extends outside of the recog patten,
       methylation can block it on EITHER of the external/extended sites: ie:
       for BsaBI: gatnnnnatc  can be blocked at EITHER a 'gatCnnnatc' or a 'gatnnnGatc', so have to check them both. */

    /* dam is a palindrome, so may interact oddly with nonpal's on the other strand, so have to do this extra weirdness */
    if (F.dam == 1) { /* only if we want to do dam sensitivity */
        for (j = F.RE_ST; j < NumREs; j++) {
            if (RE[j].dam < 100) { /* only if the RE is contextually sensitive */
                i = d = 1;
                while (i < RE[j].Sites[0] && d < RE[RE_dam].Sites[0]) {
                    if (RE[j].Sites[i] < 0) e = (RE[j].Sites[i] * -1) + (RE[j].dam * -1) + RE[j].E_olap - 4 ;
                    else e = RE[j].Sites[i] + RE[j].dam;
                    while (e > RE[RE_dam].Sites[d]) {
                        d++;
                    } /* now either equal or .Sites[d] > .Sites[i] */
                    if (e == RE[RE_dam].Sites[d]) { /* if equal, zero the entry */
                        RE[j].Sites[i] = 0;
                    }
                    i++; /* bump up i to the next el and do it again */
                }
            }
        }
    }   /* collapsing the zero entries is done in the next stanza */

    /* treat dcm the same way, even tho it's not quite the same kind of palindrome (degenerate)... - it works...  */
    if (F.dcm == 1) { /* only if we want to do dcm sensitivity */
        for (j = F.RE_ST; j < NumREs; j++) {
            if (RE[j].dcm < 100) { /* only if the RE is contextually sensitive */
                i = d = 1;
                while (i < RE[j].Sites[0] && d < RE[RE_dcm].Sites[0]) {
                    /*                              RWC        +  inv dcm offset         + the overlap  - len(dcm)  */
                    if (RE[j].Sites[i] < 0) e = (RE[j].Sites[i] * -1) + (RE[j].dcm * -1) + RE[j].E_olap - 5 ;
                    else e = RE[j].Sites[i] + RE[j].dcm;
                    while (e > RE[RE_dcm].Sites[d]) d++; /* now either equal or .Sites[d] > .Sites[i] */
                    if (e == RE[RE_dcm].Sites[d]) RE[j].Sites[i] = 0; /* if equal, zero the entry */
                    i++; /* bump up i to the next el and do it again */
                }
            }
        }
    }   /* collapsing the zero entries is done in the next stanza */


    /* This next stanza takes care of removing the last jitter from the RE[].Sites.  If there is substantial
       overlapping cutting, such as '-pname,wwwwwwwwwwwwwwww,1' (thanks Wes), there will be a very large amount of
       jitter that remains even after the check in Cutting() that checks to see if the previous cut point is the
       same as the next (ie .. 213 214 213 215 213 216 214 213 ...-> .. 213 213 213 213 214 214 215 216 ... ->
       .. 213 214 215 216 ..) so this next stanza definitively trues the array, leaving a max of one name per base. */

    for (i=F.RE_ST; i<NumREs; i++) { /* start at F.RE_ST, NOT 0, fool! (and NumREs has been incr'ed above to account for AI) */
        /* if it's a proto && there are non-protos immediately following || a gonzo degeneracy (-D > 2)  */
        if ((RE[i].E_Ncuts > 0) && (F.dam == 1 || F.dcm == 1 || (F.Degen > 2 || F.Regex != 0 || (i == RE[i].proto && RE[i+1].proto == i)))) {
            /* sort the cutsites to remove the jitter, ALLOWING NEGATIVES (rev orientation)
               to be sorted as abs(values) in abscompare so I can map them correctly on the output */
            qsort(RE[i].Sites+1, (size_t)RE[i].E_Ncuts, sizeof(long), abscompare);  /* now they're sorted */

            if ((F.dam == 1 || F.dcm == 1) && RE[i].Sites[1] == 0) {
                j = 1;
                while (RE[i].Sites[j] == 0) j++;
                k = 0;
            } else {
                j = 2;
                k = 1;
            }
            while (j <= RE[i].E_Ncuts) {
                if (RE[i].Sites[j] != RE[i].Sites[k]) {
                    RE[i].Sites[++k] = RE[i].Sites[j++];
                } else {
                    j++;
                }
            }

            RE[i].Sites[0] = k + 1;    /* correct RE[].Sites[0] to the new values */
            RE[i].Sites[k+1] = 0;      /* just put a 0 there JIC */
            RE[i].E_Ncuts = k;         /* reset the Ncuts counter to correct for the de-jittering */
        }
        /* calloc Orient w/ size Sites/8 + 1 cuz it's a bit array  */
        RE[i].Orient = (char *) calloc(((int)(RE[i].E_Ncuts/8))+1, sizeof(char));
        if (RE[i].Orient == NULL) BadMem("!! calloc failed: RE[i].Orient\n", 1);
        /* and may as well transfer the orientation data to RE[].Orient */
        for (e=1; e<=RE[i].E_Ncuts; e++) {
            if (RE[i].Sites[e] < 0) { /* set the corresponding bit in RE[].Orient */
                m = BitArray(RE[i].Orient, e-1, 1); /* the bit get SET if if's in the oppo orientation */
                /* and then swap the negative value with the abs() value so that further
                   calcs don't get screwed up */
                RE[i].Sites[e] = labs(RE[i].Sites[e]);  /* now everything's back to the way it was ... */
                /* ..except that we've recorded the orientation data in the bitarray RE[].Orient */
            }
        }
    }

    qsort(DD->Dat+2, (size_t)(((DD->Dat[0])/2) - 1), (2*sizeof(long)), abscompare);

    /* Now calc fragment sizes, sort in increasing size, now using REAL WORLD Coords so frags should be
       exactly right size - also have to consider topology for 1st and last cut*/
    /* NumREs-1 because the frags for RE[H] were already calc'ed above  */
    for (i=F.RE_ST; i<NumREs; i++) { /* the lower index point starts @ 'F.RE_ST' */
        if (RE[i].proto == i) { /* '==' was '=' until 10-24-98 */
            if (RE[i].E_Ncuts != 0) { /* using 'l' as index -> Ncuts+1 frags in linear, counting the 0th el*/
                l = RE[i].E_Ncuts;  /* need a modifiable variable;  */
                if (topo == 1) {  /* handling topo diffs - if it's linear... */
                    RE[i].Frags[0] = RE[i].Sites[1]; /* 1st frag will be magnitude of 1st cut */
                    RE[i].Frags[l] = seq_len - RE[i].Sites[l];/* last frag will be length of sequence - last cut */
                } else { /* otherwise it's circular and the 1st frag is the size of the above 2 combined ...*/
                    RE[i].Frags[0] = RE[i].Sites[1] + seq_len - RE[i].Sites[l];
                    RE[i].Frags[l] = 0; /* and the 'last' doesn't exist (size 0) */
                }
                if (RE[i].E_Ncuts > 1) {
                    for (j=RE[i].E_Ncuts-1; j>0; j--) {
                        RE[i].Frags[j] = labs(RE[i].Sites[l] - RE[i].Sites[l-1]); /* frag = site nearer end - next site  */
                        l--;  /* need the 'labs' because using REAL WORLD coords, it's possible to get negative values by */
                    }        /* same offset cutting RE site on diff strands in close proximity */
                }
            } else RE[i].Sites[0] = seq_len; /* if it doesn't cut, the fragment is the size of sequence */
        }
        if (topo == 0) RE[i].E_Nfrags = RE[i].E_Ncuts + 1;
        else           RE[i].E_Nfrags = RE[i].E_Ncuts ;
    }

    /******************************************************************************************
                                    Now, starting the OUTPUT phase
      for ALL output, need to wrap it in HTML, both for entire pages (HTML=0) or just the
      internal bits, with (=1) or without (=2) TOC.  There should also be some eror checking
      to detect when things like --idonly are called, so it doesn't go off and generate an HTML
      page that contains only a list of names.. altho that wouldn;t be too bad...
     *******************************************************************************************/

    /* If HTML == 0, want to generate a complete Web page so need the regular page headers */
    if (HTML == 0) {
        fprintf(stdout, " <!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\" "
                " \"http://www.w3.org/TR/REC-html40/loose.dtd\"> "
                " <HTML><HEAD><TITLE>tacg Ver 4 Output</TITLE></HEAD><body>\n");
    }
    if (HTML == 0 || HTML == 1) {  /* figure out and print the TOC (or figure it out elsewhere )  */
        /* go thru the flags - which were changed or should be represented in a TOC.  */
        /* have to check only those which result in a major output stanza, not all */
        /* Linear Map(-L), Gel(-g), Ladder(-l), Sites(-S), Frags(-F), ORFs(-O), orfmap (--orfmap), Graphics(-G),
           Proximity(-P), rule (--rule), Rules(--rules), Silent Sites(--silent), Extraction(-X), AFLP(-x) */
        GenerateTOC(datestamp);
    }

    /* Load Protos with the names of those REs that have > 0 cuts  */
    // FIXME Protos could be org'ed so that el [0] is ProtoCutters
    j=0;
    for (i=F.RE_ST; i<NumREs; i++) { /* start at F.RE_ST, stupid */
        if (RE[i].proto == i && RE[i].E_Ncuts > 0) {
            Protos[j++] = i; /* if it's a proto && has more than 0 cuts */
        }
    }
    ProtoCutters = j;

    /* ALL Output is dependent on whether the --idonly flag is set.  If set, ONLY the name, ID line is
       emitted; all the rest of the output is ignored. */

    if ((F.IDonly > 0 && ProtoCutters > 0) || F.IDonly == 0) { /* filter out non-hits if -i says so */
        if (HTML > -1) {
            fprintf(stdout, "<font color=\"#ee0000\"><big>");
        }
        if (SEQinfo != NULL) {  /* ID info */
            if (F.IDonly < 2)  fprintf(stdout, "\n\n");
            if (SEQinfo->idlist == 0x0) SEQinfo->idlist = "None"; /* protection on solaris */
            if (strcmp(SEQinfo->format,"Plain") == 0 || F.RawSeq) {
                fprintf(stdout, "## Sequence: #%d;  Raw sequence; no description\n",
                        SEQNum);
            } else if (F.infile != NULL) { // if there's a file name available then punt it out.
                fprintf(stdout, "## Sequence: #%d; from file: [%s] \n   Format: %s; ID: %s; Description: %s. \n",
                        SEQNum, F.infile, SEQinfo->format, SEQinfo->idlist, SEQinfo->description);
            } else {
                // else just punt to 'UNAVAILABLE'
                fprintf(stdout, "## Seq: #%d; %s %s",
                        SEQNum, SEQinfo->idlist, SEQinfo->description);
            }
            if (F.Verbose >= 1 ) {
                for (i=0; i<argc; i++) fprintf(stdout, "argv[%d]= %s\n",i,argv[i]);
                fprintf(stderr, "## Sequence: #%d; File name: [%s] DB: %s; Format: %s; Description: %s. \n",
                        SEQNum, "UNAVAILABLE", SEQinfo->dbname, SEQinfo->format, SEQinfo->description);
            }
            if (F.IDonly < 2)  fprintf(stdout, "\n");
        }   /* end of ID info */
        if (HTML > -1) {
            fprintf(stdout, "</big></font>");
        }
    }

    /* This is the BEGIN of the scope of the --idonly flag  */
    if ((F.IDonly == 1 && ProtoCutters > 0) || F.IDonly == 0) {
        /* Start with the sequence info line regardless */
        if (HTML > -1) fprintf(stdout, "<H3>");
        //fprintf (stdout, "== Sequence info:\n");
        if (seq_len > (l = F.NBases[0] + F.NBases[1] +F.NBases[2] + F.NBases[3])) {
            fprintf (stdout, "\n    NB: sequence length > A+C+G+T due to -> %ld <- IUPAC degeneracies.\n", seq_len-l);
            fprintf (stdout,
                     "    # of:  N:%ld  Y:%ld  R:%ld  W:%ld  S:%ld  K:%ld  M:%ld  B:%ld  D:%ld  H:%ld  V:%ld \n\n",
                     F.NBases[14], F.NBases[5], F.NBases[6], F.NBases[7], F.NBases[8], F.NBases[9], F.NBases[10],
                     F.NBases[10], F.NBases[11], F.NBases[12], F.NBases[13]);
        }
        if (HTML > -1) fprintf(stdout, "</A></H3>\n");
        //fprintf (stdout, "   #s below are for top strand; 'sites exp' values calculated on the basis of both strands.\n");
        fprintf (stdout, "   %ld bases; %ld A(%.2f %%)  %ld C(%.2f %%)  %ld G(%.2f %%)  %ld T(%.2f %%)",
                 seq_len, F.NBases[0], FrACGT[0]*100, F.NBases[1], FrACGT[1]*100, F.NBases[2], FrACGT[2]*100, F.NBases[3], FrACGT[3]*100);

        /* EXAMPLE: This is just the test of whether the example fn() flag is set; if it is, it just emits a short note to
           that effect.  If not, it skips it. No big deal either way...  The only reason for putting the test right
           here as oppo to anywhere else is that at this point, all the data has been generated and it's only whether
           you want to do anything special with the data at this point (via the code that you might insert here) */
        if (F.Example != 0) {
            fprintf(stderr, "\nEXAMPLE output: The example flag value is: %d \n\tand the number of bases divided by it is: %f .\n\n", F.Example, (float)(seq_len/F.Example));
        }

        /* Output a Summary table of REs that do not cut, if -s flag is set */
        if (F.Summary == 1) {  /* only output the summary if really wanted */
            fprintf(stdout, "\n\n");
            if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"nocuts\">");
            fprintf (stdout,"== Enzymes that DO NOT MAP to this sequence:\n");
            if (HTML > -1) fprintf(stdout, "</A></H3>\n");

            reps = (int)((O_LMARG+O_RMARG+basesPerLine)/10 - 1);
            i=0;
            j=0;
            Cur_RE = F.RE_ST;  /* start off the index at F.RE_ST, not 0, not 1 */
            while (Cur_RE < NumREs ) {
                if (RE[Cur_RE].E_Ncuts == 0 && RE[Cur_RE].proto == Cur_RE) {
                    fprintf(stdout, "%*s  ", F.MaxPatNameLen, RE[Cur_RE].E_nam);
                    j++;
                    i++;
                }
                if (j == reps) {
                    fprintf (stdout, "\n");
                    j=0;
                }
                Cur_RE++;
            }
            if (i==0) fprintf(stdout, "\n\tThere were NO NON-matches - ALL patterns matched at least ONCE.\n");

            /* Now, just output a table of the number of cuts by enzyme if they have 1 or more cuts */
            fprintf(stdout, "\n\n");
            if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"allcuts\">");
            fprintf (stdout, "== Total Number of Hits per Enzyme:\n");
            if (HTML > -1) fprintf(stdout, "</A></H3>\n");

            reps = (int) (((O_LMARG+O_RMARG+basesPerLine)/15)-1);
            m = (int) (ProtoCutters/reps)+1;
            for (i=0; i<m; i++) {
                for (j=0; j<reps && i+(j*m)<ProtoCutters; j++) {
                    fprintf(stdout, "  %*s %7d", F.MaxPatNameLen, RE[Protos[i+(j*m)]].E_nam, RE[Protos[i+(j*m)]].E_Ncuts);
                }
                fprintf(stdout,"\n");
            }
        }

        /* If hookey() has been requested, spit it out. */
        if (F.Hookey != -1) { /* if we've detected a '=' in the '-x' stream */
            /* want to check that the 'tagged' RE has cut at least once, else no detectable frags */
            if (RE[F.Hookey].E_Ncuts > 0) Hookey(NumREs++, AI, seq_len);
        }
        reps = (int) ((O_LMARG+O_RMARG+basesPerLine)/7); /* needed for both PrintSites() and PrintFrags() */
        /* Use Protos instead of GREs - should be close to a 1:1 correspondence */
        /* Use this gruntsort routine to sort GREs in place according to # cuts for each enzyme-
           don't have to change a thing then for the site/fragment/ladder/gel output!!  Otherwise,
           it's just dumped by order of their listing in REBASE */
        if (F.StriderSort==1) { /* if want the frag, sites, etc listed by number of cuts (-c flag)*/
            for (j=0; j<ProtoCutters; j++) {
                for (k=j+1; k<ProtoCutters; k++) {
                    if (RE[Protos[j]].E_Ncuts > RE[Protos[k]].E_Ncuts) {
                        mm = Protos[j];
                        Protos[j] = Protos[k];
                        Protos[k] = mm;
                    }
                }
            }

            /* and now re-sort on basis of RE names so that it will be alphabetic within cut bins */
            j=k=b=e=0;
            while (e < ProtoCutters) {
                while (RE[Protos[e]].E_Ncuts == RE[Protos[j]].E_Ncuts) e++;
                for (j=b; j<e; j++) {
                    for (k=j+1; k<e; k++) {
                        if (strcmp(RE[Protos[j]].E_nam , RE[Protos[k]].E_nam) > 0) {
                            mm = Protos[j];
                            Protos[j] = Protos[k];
                            Protos[k] = mm;
                        }
                    }
                }
                b=e;
            }
        }

        /* If wanted, output the extracted sequences (-X/--extract b,e,revcomp) */
        if (F.Xtract != 0) ExtractSeqAtHit(SEQinfo, seq_len, sequence, ProtoCutters, Protos);

        /* And if wanted (-S), print out all the cut sites, filtered as below in PrintFrags() */
        /* What does ProtoCutters equal??  Does it include the combo cuts?  the Hookey data? */
        if (F.Sites > 0) PrintSitesFrags(ProtoCutters, reps, Protos, 1, seq_len);

        /* if asked for -P matching, do the name site matching calculations */
        if (F.Prox > 0) {
            Proximity();
            PrintGelLadderMap((int)F.Prox, seq_len, Protos, 4);
        }

        /* Now if wanted, print out the UNsorted fragments generated by the cuts, listed across the page.
           This routine also uses a very large number of printf calls to generate the output format and
           might be better written to 'compose and dump' entire lines rather than doing it piecewise
           as it does here - advice is welcome. */

        mm = (int)F.Frags;
        if (mm==1 || mm==3) { /* if want *unsorted* or *both sorted and unsorted * frags */
            /* can use the 'sites' (var 4) to indicate not only that whether its to print Fragments or Sites, but
                also sorted or unsorted fragments sites = 1, unsorted frags = 0, sorted frags = -1  */
            PrintSitesFrags(ProtoCutters, reps, Protos, 0, seq_len);
        }

        if (F.ORFs != -1) { /* if want a table of ORF data */
            fprintf(stdout, "\n\n");
            if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"orfs\">");
            fprintf (stdout, "== Open Reading Frame Analysis:\n");
            if (HTML > -1) fprintf(stdout, "</A></H3>\n");

            for (i=0; i<6 && F.ORFs2Do[i] != 0; i++) { /* jjj - did have a ',' instead of the '&&' */
                if (F.ORFs2Do[i] >3 && F.ORFs2Do[i] <7) {
                    /* have to get mem for and antipar sequence into a temp holder to submit to Translate, etc */
                    /* l = strlen(sequence+1); */ /* commented out 10-23-98 */
                    TmpSeq = calloc((size_t)tot_seq_Cnt+1, sizeof(char)); /* dmalloc: need the '+1' to avoid smashing pickets */
                    if (TmpSeq == NULL) BadMem("tacg.c:~822 TmpSeq calloc err", 1);
                    Anti_Par(sequence, TmpSeq, tot_seq_Cnt); /* now TmpSeq is pointing at the AP'ed sequence */
                    /* also have to do some muddling with frame because have to convert frames 456 to 123 and
                       pass them correctly */
                    if (F.ORFs2Do[i] == 4) moo = 3;
                    else if (F.ORFs2Do[i] == 5) moo = 1;
                    else moo = 2;
                } else {
                    TmpSeq = sequence;  /* else point TmpSeq at sequence as a temp pointer */
                    moo = F.ORFs2Do[i];
                }
                thisframe = (int)F.ORFs2Do[i] ;
                Prot4ORF = (char *)calloc((seq_len/3 + 10), sizeof(char));
                if (Prot4ORF == NULL) BadMem("tacg.c:~653 Prot4ORF calloc err", 1);
                if (moo==3) moo = 0; /* horible hack */
                Translate(TmpSeq+BASE_OVERLAP+moo, Prot4ORF, seq_len-moo, 2, Codons, codon_Table);

                if (F.Orfmap > 0) { // if want met/stop map, process these ORFs into a met/stop structure
                    fr = thisframe -1;
                    // need to store position, whether it's a met (tho some orgs will start with a non-met) or a stop.
                    //simplest structure is a 3D array frame x position x start/stop - unfortunately, will have to be
                    //an int array to capture # of aas int *metstop[6] - allocate it 1st here (only if we need it)
                    //and only for those frames that were called for. The array only notes positions -
                    // '+' if start, '-' if stop
                    if ((metstop[fr] = calloc(ORFstep, sizeof(int))) == NULL) {
                        BadMem("Can't init mem for metstop", 1);
                    }
                    metstop[fr][0] = ORFstep; //the end of allocated mem
                    msp = 1; // the next free el.
                    aap = 0;
                    while (Prot4ORF[aap] != '\0') {
                        // this offsets the real index by 1 so that there will never be a '0' start
                        if (Prot4ORF[aap] == 'M') metstop[fr][msp++] = aap+1;
                        if (Prot4ORF[aap] == '*') metstop[fr][msp++] = (aap+1) * -1;
                        if ((metstop[fr][0] - msp) < 5) {
                            newsize = metstop[fr][0] + ORFstep;
                            if ((metstop[fr] = (int *) realloc(metstop[fr], sizeof(int) * (metstop[fr][0] + ORFstep))) == NULL) {
                                BadMem("Can't realloc mem for metstop", 1);
                            }
                            metstop[fr][0] = newsize;
                        }
                        aap++;
                    }
                    metstop[fr][msp] = 0; //terminate the fucker so it can't go on and on and on
                }

                ORF_Analysis(Prot4ORF, (long)strlen(Prot4ORF), (short)thisframe);
                free(Prot4ORF);
                if (i >2 && i <6) free(TmpSeq);
            }
        }

        if (F.Orfmap > 0) {  /* if want a text-ladder ORF map */
            if (F.ORFs != -1) {
                PrintORFMap(seq_len, metstop);  /* separate fn() for printing ORFs ...?  or insert code into PrintGelLadderMap()? */
            } else {
                fprintf(stderr, "\n\nSyntax error: No ORFs (-O) were specified for analysis, so '--orfmap' cannot be provided\n\n");
                exit(1);
            }
        }


        if (mm==2 || mm==3 || (gel==1)) { /* if want _sorted_ or if need sorted for gel */
            for (i=0; i<ProtoCutters; i++) { /* the lower index point starts @ 'F.RE_ST', not '10' because we're using GREs[] */
                if (RE_Filters_OK(Protos[i])) {
                    b = 0;
                    if (topo == 0) b = 1;
                    qsort(RE[Protos[i]].Frags, (size_t)RE[Protos[i]].E_Nfrags-b, sizeof(long), compare);
                }
            }

            if (gel) {  /* use sorted frags to generate the gel */
                PrintGelLadderMap(ProtoCutters, seq_len, Protos, gel);
            }

            if (mm==2 || mm==3) { /* Now the frags are sorted! Just use PrintFrags again (w/ sites = -1) to dump them out */
                PrintSitesFrags(ProtoCutters, reps, Protos, -1, seq_len);
            }
        }

        if (F.Ladder == 1) {  /* if want a text-ladder map of the sequence */
            PrintGelLadderMap(ProtoCutters, seq_len, Protos, 0);  /* subst 'ProtoCutters' for 'NumGREs', 'Protos' for 'GREs' */
        }

        /* if want the RE data streamed out as data appropriate for graphing */
        if (F.GrfBBin > 0)  DumpDataForPlot(seq_len, ProtoCutters, Protos);

//
        /* if user wants to do the SW routines and the seq > window size, then do it, but am
           also assuming this indirectly from the setting of F.Rules, as well as F.SlidWin.
           users might want to use the sliding window feature as part of the regular tacg
           features and NOT with a set of rules from '--rule/--rulefile'.  Will this support Sliding
           Windows WITHOUT any rules???  I think not.  Either the fn() has to be broken into 2
           explicitly or logically to enable the functionality to be split.
           They might also want to use  the '--rule(s)' feature with entire sequences, not
           particularly with the sliding window feature (this is easily addressed by setting the
           window to the sequence length */
        /* this needs to be neatened up a bit; there's still cruft left over form the massive surgery, but
           it'll get straightened out as debuggin goes along.. */

        if ((F.Rules == 1 || F.SlidWin != 0) && seq_len >= F.SlidWin) {
            if (F.SlidWin == 0) F.SlidWin = seq_len;
            /* should now write this completely as a fn() from the ground up, so... */
            NumPosSlidWin = tacg_SlidingWindow(0, NumREs, 1, SelEnz, Protos, NumProtos);
            printf("\n### NumPosSlidWin = %d \n", NumPosSlidWin);
            /*                                 ^                                      */
            /* for a single rule (entered from commandline) only need to eval Rule[0] */
        }

        /* If F.Rules > 1, want to evaluate several rules entered from the Rule file */
        if (F.Rules > 1) {
            /* need to hash the RE names into an array so that we can check against it to see if
               all the nec patterns are loaded. */
            /* set el = RE index, NOT Proto index; need even those REs that DON'T cut (Protos filters
               non-cutters out - see ~line 574) */
            for (rule=0; rule<F.NumRules; rule++) {
                /*  printf(stderr, "\n\nRule  to be eval'ed:\n%s\n", Rules[rule].Name);  */
                NumPosSlidWin = tacg_SlidingWindow(rule, NumREs, 1, SelEnz, Protos, NumProtos);
            }
        }


        /* if want to use the --clone functions */
        if (F.Clone > -1 && ProtoCutters > 1) { /*  and there are REs defined */
            /* only alloc space for Clone[].matching_REs if it looks like it'll be needed */
            for (i=0; i<=F.Clone; i++) {
                if (Clone_S[i].end > seq_len || Clone_S[i].begin > seq_len) {
                    fprintf(stdout, "\n\nSorry - in trying to fulfill the '--clone' request, rule %d specifies \n"
                            "one or both endpoints (%s) greater than the length of the \n"
                            "sequence you've submitted (%ld).  Please correct and try again.\n\n",
                            i, Clone_S[i].range_rule, seq_len);
                    exit(1);
                }
                /* start with an initial size of 5 to force test the re-alloc */
                if ((Clone_S[i].matching_REs = (long *) calloc(CLONE_INITIAL_MRE, sizeof(long))) == NULL) {
                    BadMem("Can't init mem for Clone_S[].matching_REs", 1);
                } else {
                    Clone_S[i].matching_REs[0] = CLONE_INITIAL_MRE; /* current alloc range */
                    Clone_S[i].matching_REs[1] = 2;                 /* 1st free el to use */
                }
            }
            Clone(ProtoCutters, Protos);
        }


        /* if you want a linear map of the sequence */
        if (F.LinMap == 1) {
            LinearMap(seq_len, basesPerLine, NumREs, sequence, max_okline, Codons, codon_Table, Protos);
        }

        /* if want a postscript PLASMID map  */
        if (F.GrafMap != -1) {
            if (SEQinfo->description != NULL) {
                strncpy(s, SEQinfo->description, 255);
                s[255]='\0';
            } else {
                s = "Sequence has no header";
            }
            GrafMap(Degen_Log, Protos, s, datestamp, FrACGT);
        }

        /* if want a postscript LADDER map  */
        if (F.psladder != 0) {
            if (SEQinfo->description != NULL) {
                strncpy(s, SEQinfo->description, 255);
                s[255]='\0';
            } else {
                s = "Sequence has no header";
            }
            tacg_psladder(Protos, s, datestamp);
        }


    } /* This is the End of the scope of the --idonly flag  */

    /* Following allows multiple analyses to be appended without getting them all mixed up*/
    /* Now, HAVE to be careful with vars that get re-used... like ALL of them! What a pain, but it'll work!!  */

    //fprintf(stderr, "pre-free, sequence is at: %x  and is of size %d \n", &sequence, sizeof(sequence)); */

    //free(sequence);  /* give back the mem! */
    //memset(sequence, 0, seq_len); //no need to free  it, just clear it
    F.NBases[0] = F.NBases[1] = F.NBases[2] = F.NBases[3] = 0;


    if (F.Hookey != -1) {  /* Norm Drinkwater's fix  */
        mm = RE[F.Hookey].E_Ncuts;  /* see below  */
    }

    /* and of course all the vars in RE that are sequence-dependent */
    for (i = F.RE_ST; i < NumREs; i++) {
        if (RE[i].proto == i) {
            free(RE[i].Sites);
            free(RE[i].Frags);
            RE[i].EstSites  = 0.0;
            RE[i].E_Ncuts = 0;
        }
    }
    if (F.Hookey != -1) {
        NumREs--; /* bc otherwise, the 'combined' routines will keep incr with each new seq  */
        if (mm > 0) NumREs--; /* ditto */
    }

    /* should we free() DD->Dat or just re-zero it?  What would be more efficient? Probably just rezeroing it,
       since otherwise we have to both free() and re calloc it, BUT it would have to grow and then stay to the largest
       size seen so far...not too bad, really, but do the free/re-calloc 1st cuz it's easier.. */
    free(DD->Dat);
    DD->Dat = (long *) calloc(10000, sizeof(long)); /* init the pointer to 10,000 */
    if (DD->Dat== NULL) BadMem("calloc failed -> DD->Dat mem!", 1);
    DD->Dat[0] = 2;      /* this is Cnt and it has to start at 2 to skip these two unholy counters */
    DD->Dat[1] = 10000;  /* this is Siz */
    /* want to memset with a 0 (zero) rather than a ' ' */
    memset(SelEnz, 0,(MAX_SEL_ENZ*sizeof(*SelEnz)));
    if (F.Verbose >= 1) fprintf(stderr, "Finished sequence # %d\n", SEQNum);

    if ((F.RawSeq == 0) && ((SEQseq = seqfgetseq(SEQfp, &SEQlen, 1)) != NULL)) {
        /* get the sequence via SEQIO and munge it into the standard format with PadEndForEnd(),
           pass it thru the rest of tacg, and then come back for more until there's no more sequences */
        /* entry can be either thru seqio or via raw sequence input  */
        /* SEQNum = 0; */ /* why was this ever functional?? */
        SEQinfo = seqfinfo(SEQfp, 0);         /* this one doesn't */
//               fprintf(stderr,"Got another seq from seqfgetseq\nSEQinfo->description: %s\nLen=%d", SEQinfo->description, SEQlen);
        seq_len = (long)SEQlen;
        free(sequence);
        sequence = 0x0;
        sequence = GetSequence2(SEQseq, &tot_seq_Cnt, &seq_len, &Degen_Log);
        free(SEQseq); /* as soon as it's moved into sequence, don't need the original anymore */
        MoreSequence = 1;
    } else MoreSequence = 0;
}  /* ##################   End of the SEQIO while() loop ######################### */
// if we need to generate the whole page,close it out, so
if (HTML == 0) fprintf(stdout, "<br><pre></body></html>\n");

// Free the unfreed
// if modify to do multiple rules, have to do a recursive free() down the number of rules.
// only .RuleString is allocated, Name and Window are static
if (Rules!= NULL) free(Rules[0].RuleString);
if (Protos != NULL)free(Protos);
if (DD->Dat != NULL)free(DD->Dat);
if (DD != NULL)free(DD);
if (SelEnz != NULL)free(SelEnz);
//if (Fbuf != NULL)free(Fbuf); // no context anymore for Fbuf?
if (rulestring != NULL)free(rulestring);
if (Degen_Log != NULL)free(Degen_Log);
if (SEQfp != NULL) free(SEQfp);

/* free all the hashtable mem - not nec req for proper function, but if tacg eventually moves to a fn().. */
if (F.Matrix !=0 && F.Regex !=0) {
    for (i=0; i<4096; i++) free(hashtable[i]);
}
/* and free() the mem for that particular orf or will leak like a sieve */
i=0;
while (F.ORFs2Do[i] != 0) {
    j=0;
    while (ORFs[F.ORFs2Do[i]-1][j].orflength != 0) {    /* ORFs[i][j] */
        free(ORFs[F.ORFs2Do[i]-1][j++].orf);  /* gotta do this !! */
    }
    free(ORFs[i++]);
}

if (F.PDF) {
    gs_bin = "/usr/bin/gs";
    putenv("TMPDIR=/tmp"); /* is this more portable than setenv? - solaris didn't like setenv */
    /*      setenv("TMPDIR", "/tmp", 1); */
    if (stat(gs_bin, Fbuf) == -1) { /* if failure, die politely */
        gs_bin = getenv("TACG_GS_BIN");
        fprintf(stderr, "\n\nNo ghostscript/gs at '/usr/bin/gs'. Trying environment variable TACG_GS_BIN (%s)...\n",  gs_bin);
        if (stat(gs_bin, Fbuf) == -1) {
            fprintf(stderr, "\n\nNo ghostscript/gs at %s either; No gs, no plasmid maps. Sorry.", gs_bin);
            fprintf(stderr, "\nTry setting the TACG_GS_BIN environment variable to the fully qualified \n"
                    "path name of your ghostscript binary.\n");
            exit(1);
        }
    }
    fprintf(stderr, "\nF.tmppath = %s\n", F.tmppath);
    //compose the cmd
    sprintf(gs_cmd, " %s -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -q -sOutputFile=%s/tacg_Map.pdf  %s/tacg_Map.ps" , gs_bin, F.tmppath,  F.tmppath);
    fprintf(stderr, "%s\n", gs_cmd);

    fprintf(stderr, "\n\nPDF Output: Don't forget - tacg APPENDS the newest page to the END of the \n"
            "postscript & PDF files!\n"
            "The files are named 'tacg_Map.pdf' and 'tacg_Map.ps' in the current directory.\n\n");
    fflush(stdout); /* make sure that tacg and gs output don't get mixed up! */
    sys_err = system (gs_cmd);//      system("exit(0)");
    sleep(1);
    if (sys_err == -1) {
        fprintf(stderr, "\n\nsystem() call to  Ghostscript failed.\n\n");
    }
}
/* free the bits that should be freed  */

for (i=F.RE_ST; i<NumREs; i++) {
    free(RE[i].Orient);
    if (RE[i].E_wsit[0] != NULL) free(RE[i].E_wsit[0]);
    free(RE[i].E_wsit[1]);
}
if (sequence != NULL) free(sequence);
if (RE != NULL)free(RE);  //but would have to free all the bits of RE as well to free it all.

/* the following free() sequence is excessive for a single pass app, but is formally a good thing */
if (hashtable[0] != NULL) {
    for (i=0; i<4096; i++) {
        free(hashtable[i]);
    }
}
/* The SEQIO loop should go down to here - would really like it to cite the analyses of
    multiple sequences.  If it slows things down too much, I'll just have it done once at the end...
*/

return 0;
} /* End main() */


/* tacg_SetF does just that - sets the F struct values for all the flags in an isolated
   way so that it doesn;t contribute code lines to tacg per se.  */
void tacg_SetF(void)
{
    int i;

    /* longs */
    F.Seq_Len=0;
    F.SeqBeg = 1, F.SeqEnd = 0;
    F.NumStart=0;
    for (i=0; i<15; i++) F.NBases[i] = 0; /* rather than FrA, FrC, FrG, FrT */

    /* chars */
    F.RuleFile = '\0'; /* not really nec to set this; it's a char *, but put something there JIC  */
    if ((F.tmppath = calloc(128, sizeof(char))) == NULL) {
        BadMem("tacg_SetF:tmppath", 1);
    }
    strcpy(F.tmppath,".");
    /*shorts */
    F.Topo=1, F.Overhang=1, F.Overhanglen = 0, F.Mag=3, F.Ladder=0, F.Clone=-1,
      F.psladder = 0, F.Orfmap = 0, F.Timeops=0;
    F.XLinFr=0, F.Lab1or3=0, F.ORFs=-1, F.Rev = 0, F.Comp = 0, F.RevComp = 0, F.RevNum = 0,
      F.Version = -1, F.Help = -1, F.Degen = 1,  F.Xplicit = 0, /* -X */
        F.AltPatFile = -2,  F.CodTab = 0, F.Frags = -1, F.LinMap = 0, F.GrafMap = -1;
    F.LogDegens = 0, F.PDF = 0, F.Sites = 0, F.HTML = -1, F.StriderSort = 0, F.Verbose = 0,
      F.Summary = 0, F.Pat = 0, F.Prox = 0;
    /* remember - F.ORFs2Do: 0 -> 5, NOT 1 -> 6 */
    for (i=0; i<6; i++) F.ORFs2Do[i] = 0;   /* ORFs2Do[0-5] = flags[31-36] */
    F.Err = 0, F.GrfBBin = 0, F.GrfFmt = 0, F.Output = 0, F.DB = 0,  /* not implemented yet or rather de-implemented */
      F.Hookey = -1, F.Matrix = 0, F.AllInc = 0, F.Regex = 0,
        F.Xtract = 0, F.XtBeg = 0, F.XtEnd = 0,
          F.Silent = 0, F.IDonly = 1, F.Example = 0, F.OrfXtra = 0,
            F.Rules = 0,  /* 0 = not set, 1 = use the --rules entered string, 2 = use Rules from Rules File */
              F.NumRules = 0, /* NumRules is the number of rules to be evaluated from a file of rules */
                F.SelEnzUsed = 0, F.RE_ST = 1;
    F.dam = 0, F.dcm = 0;
    F.Strands = 2;
    F.Tics = 1;
    F.RawSeq = 0;
    F.GrafMap = -1;
    F.nmer = 0;
    /* ints */
    F.Min = 1, F.Max = 32000,
      F.SlidWin = 0;
    F.Width = 60,
      F.GelLO = 0, F.GelHI = 0;
    F.NumSubPatts = 0; /* number of subpatterns found in the Rules file */
    F.Cost = 0,
      F.PID = 0;
}
