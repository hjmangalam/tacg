/*  $Id: tacg.h,v 1.4 2005/01/26 22:44:08 mangalam Exp $

...And The Small Print Taketh Away... 

COPYRIGHT NOTICE
----------------
 Copyright (c) 1994-2005 by Harry Mangalam, tacg Informatics
 (hjm@tacgi.com, 949 856 2847)


    Permission to use, copy, modify, and distribute this software and its
    documentation is hereby granted, subject to the following restrictions
    and understandings:

      1) Any copy of this software or any copy of software derived
         from it must include this copyright notice in full.

      2) All materials or software developed as a consequence of the
         use of this software or software derived from it must duly
         acknowledge such use, in accordance with the usual standards
         of acknowledging credit in academic research.

         The use of tacg may be cited as:
         Mangalam, HJ. (2002) tacg, a grep for DNA. BMC Bioinformatics. 3:8
         http://www.biomedcentral.com/1471-2105/3/8
         tacg Informatics. hjm@tacgi.com 949 856 2847

      3) The software may be used by anyone for any purpose, subject to the
         Free Software Foundation's General Public License, a copy of which
         should have been included with this distribution (see the file
         COPYING).
         If it was not included with this distribution, a copy may be obtained
         from: http://www.gnu.org/copyleft/gpl.html .

         Those persons or organizations wanting to use code derived from tacg
         in proprietary applications so as not to have their code 'infected'
         by the FSF GPL have the option to license the code directly from from
         the author under a separate license.  This option takes advantage of
         the  so-called Dual Licensing policy, which is approved by the FSF.

      4) This software is not meant to be used as a sole tool to design
         recombinant molecules for theraputic use.  While it is accurate to
         the best of my knowledge about the results it gives, it has not
         been validated to this stringency yet.  Hence the release of the
         source code.

      5) This software is provided AS IS with NO WARRANTIES OF ANY KIND.  The
         author shall have no liability with respect to the infringement of
         copyrights, trade secrets or any patents by this software or any part
         thereof.  In no event will the author be liable for any lost revenue
         or profits or other special, indirect and consequential damages.

		Inclusion of other software
      ===========================

      This software uses code from other authors and packages that are
      copyright  of the individual authors.  Their release policies are not in
      conflict with  those of the GPL, as far as I can tell.  The externally
      developed code is described below.

      The use of the included code does not imply an endorsement of the tacg
      package by the individual authors.

      This software uses the SEQIO package for reading and writing sequences.
      Copyright (c) 1996 by James Knight at Univ. of California, Davis.  See
      the file 'seqio.c' for his release policy.

      This software uses the 'strsep' routines, Copyright (c) 1998 by Gray
      Watson. See the file 'strsep.c' for his release policy.
      NB: strsep() is now included in many stdlibs.

*/

#include "seqio.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

/* this stanza sets up PCRE to direct which bits to execute, BUT it still requires
   conditional compilation to skip the REGEX fn()s so the makefile doesn't choke
   on missing libs and that hasn't been done yet.  */
#ifdef HAVE_PCRE
#define PCRE 1
#elif defined HAVE_LIBPCRE /* configure defines dif vars on dif platforms.. */
#define PCRE 1
#else
#define PCRE 0
#endif


/* try to use a more portable ghostscript-finding code
#ifdef GHOSTSCRIPT
#define GS_BIN GHOSTSCRIPT
#endif
*/

/* Defines */
/* VERSION is defined in configure.in -> configure, to change following one to
	something different */
#define TACG_VERSION VERSION
#define TACG_BUILD_DATE BUILD_DATE
#define TACG_BUILD_PLATFORM UNAME
#define TACG_GCC_VER GCC_VER

#define MAXSITE 200 	/* max length of site - may not need if do use pointer math intelligently */
#define MAX_NUM_RES 1000000  /* Max number of REs that we'll accept for text output routines - will later
                       			be taking some of these from the command line - may be made archaic by dynamic mem*/
#define INIT_RES 250   	/* Initial # of REs for use in dynamic mem allocation */
#define RE_INCR 200     /* the size by which RE is increased at each step */

#define BASE_OVERLAP 50 /* # of bases past the end of a block in order to 'trap' REs whose cut offset put them
                        	into block - was 20, see if it works at 30 - yes!  Also used in topology routines as the
                        	amount that needs to be added to beginning and end of sequence to allow 'trapping'
                        	as above, if DNA is circular.
                            NB: due to various cuations about padding the end of the buffer with a spare byte, this
                            value results in padding BASE_OVERLAP+1 bases front and back for almost all cases.  This doesn't
                            change anything and the numbering is accurate, but it's a consideration for others to note if
                            they're going to code in more functions using it. */
#define INIT_PPs 20            /* Initial number of Proximity structs to allocate */
#define INIT_CLONE_S 20        /* Initial number of Clone structs to allocate */
#define O_SEQ_LINE 190         /* # of line in MAX_OUTPUT_LINES where the sequence gets writ; should be about MAX_OUTPUT_LINES - 10 */
#define MAX_OUTPUT_LINES 200   /* max # of lines per block - also influences how many blank lines go after
                              	the data at each block; was 25, now incr to 80 to allow 6 frames of translation and all the REs you could imagine */
#define MAX_PAT_NAME_LEN 30   /* Max pattern name length for labelling */
#define O_LMARG 10            /* left margin in text output */
#define O_RMARG 10            /* right margin in text output */
#define MAX_ORGS 17           /* number (+1) of different Codon Usage tables */
#define N_CODONS 64           /* number of codons - this would ever change???? */
#define PREFERED_NUM_DIVS 10  /* # of major divisions in the ladder map */
#define TIC_REPEAT 4          /* after this many lines should the tics be repeated in the ladder map */
#define NUMBERLINE_REPEAT 20  /* after this many lines should the number header be repeated in the ladder map */
#define SUMMARY_CUTS 0        /* this many # cuts is the cutoff for the summary map - if set to 0, then it will
                               	be calculated so that >10000 bases, you'll get log(seqlen)-2 cuts reported
                               	if set to a # > 0, it will stay statically at that number  */
/* the following 2 constants should always be changed together and in fact
    one should be set to the other at some point */
#define MAX_LADGEL_WIDTH   1024   /* # chars to print out the ladder/gel map - obv only use if use a v. compr font */
#define MAX_BASES_PER_LINE 1024   /* the max allowable output width; placeholder until can grab basesPerLine
												from F. to set wanted output width; need it to be able to decl and
												init many arrays this can be bypassed only by specifying '-w 1' which
												sets F.Width to 10,000,000 allowing the output to be in one line -
												easier to handle for post-processing by other apps, all the fn()s that
												depend on a width setting use 'reps' which is calculated based on
												F.Width; O_txt[][] is based directly on MAX_BASES_PER_LINE
												so it WON'T create an array of 10Mx*MAX_OUTPUT_LINES */
#define RE_NAME_TAB_SIZ 3331   /* the (prime) size of the hashtable for handling the RE Names - */
#define MAX_SEL_ENZ 50        /* max number of enzymes that can be selected from the commandline */
#define MAX_DEGEN 256         /* max degeneracy that can be in the keying hexamer for the degenerate cutting fn() */
#define ALRDY 100             /* size of already[] in Cutting() and HorribleAccounting() */
#define MAX_PP 10             /* max # of proximity pattern relationships */
#define MAX_ERR 4             /* max # errors allowed in a sequence for -p, -P - 5th term in a REBASE entry */
#define MAX_RAW_MATRIX_LEN 50 /* max length of the RAW matrix length (in Read Matrix() */
#define INIT_DEGEN_LOG 1000   /* initial size of the Degen_Log */
#define INCR_DEGEN_LOG  500		/* calloc increment size for Degen_Log */

#define MAX_NUM_CLONE_RANGES 15  /* Max # of ranges that can be specified for --clone option */
#define CLONE_INITIAL_MRE 5   /* initial alloc for Clone[].matching_REs */
#define CLONE_INCR_MRE    5   /* incremental alloc for Clone[].matching_REs */

#define INIT_NUM_RULES   15   /* initial # of rules */
#define RULES_INCR        5   /* amount Rules[] should be incremented at a shot */
#define MAX_NUM_RULES  100    /* initial # of rules */
#define MAX_LEN_RULESTRING 2000  /* yup */

#define NUM_PLASMID_LABELS 180 /* number of slots in the label array HP (Holding Pen) */

#define _XOPEN_SOURCE 500


/* ------------------ Global declarations ------------------------ */

extern struct RE_struct *RE;
extern struct ORF_struct *ORFs[6];
extern struct Digest_Data *DD;
//extern struct Prox_struct PP[10];
extern struct Prox_struct *PP;
extern struct PW_struct *PosWin;
//extern struct Clone_struct Clone_S[MAX_NUM_CLONE_RANGES];
extern struct Clone_struct *Clone_S;
extern struct Rule_struct *Rules;
extern struct flag_struct F;
extern struct SE_struct *SelEnz[];
extern char *sequence;
extern long *Degen_Log;
extern char *optarg;
extern int optind, opterr, optopt;


/*  --------------------- Struct Declarations --------------------------  */


struct NextHit_struct {
   double alpha, /* incremental (middle) angle to the next hit site, assuming an upper right quadrant angle */
   		 delta, /* lower angle  */
          beta, /* upper angle */
          radius,
          xc, yc,  /* the origin / center of the plasmid circle */
         xet, yet, /* end of tic coords */
          x1, y1,  /* the coordinates of the last hit */
          xn, yn;  /* the coordinates of the next hit - what is 'returned' from the function */
};

struct SE_struct {
   char  PName[MAX_PAT_NAME_LEN+1];
/*    char  PName[11];  */
   short match,  /* whether the Name gathered in SetFlags matched a pattern
                    name in a REBASE, regex, or matrix file */
   hookey;       /* whether this name is a hookey enzyme */
   int   min, Max;  /* the per pattern min, Max values, extracted from --rules
                       and entered in the Read....c files */
 /*   int REi;  the index to the RE[] when doing a 'reverse lookup'  */
};

/* struct SE_struct SelEnz[MAX_SEL_ENZ]; */

/* this is the struct that keeps track of all the Flags for all the options.  Many of the options
   currently allocated as 'short' could also be allocated as bit flags with the bit masking option */
struct flag_struct {
   char *RuleFile, *tmppath, *infile;
   short Topo, Overhang, Overhanglen, Mag, Ladder, Clone, psladder,
      XLinFr, Lab1or3, ORFs, MaxPatNameLen, Orfmap, Timeops,
      Rev, Comp, RevComp, RevNum,
      Version, Help,
      Degen,
      Xplicit,
      AltPatFile,
      CodTab,
      Frags,
      LinMap,
      LogDegens,
      Sites,
      HTML,
      StriderSort,
      Verbose,
      Summary,
      Pat,
      Prox,
      /* remember - F.ORFs2Do: 0 -> 5, NOT 1 -> 6 */
      ORFs2Do[6],  /* ORFs2Do[0-5] = flags[31-36] */
      Err,
      GrfBBin, GrfFmt,
      Output,
      DB,  /* not implemented yet or rather de-implemented */
      Hookey,
      Matrix,
      AllInc,
      Regex,
      Xtract, XtBeg, XtEnd,
      Silent,
      IDonly,
      Example,
      Rules, NumRules,
      OrfXtra,
      SelEnzUsed,
      RE_ST,
      dam,
      dcm,
      Strands,
      Tics,
      RawSeq,
      GrafMap,
      PDF,
      nmer;

   int   Min, Max, PID,
      SlidWin,
      Width,
      GelLO, GelHI, NumSubPatts, Cost;
   long  Seq_Len, SeqBeg, SeqEnd, NumStart,
                    /*             0 1 2 3   4 5 6 7 8 9  10 11 12 13 14  */
      NBases[15];   /* for each of a c g t   y r w s k m   b  d  h  v  n (rather than FrA, FrC, FrG, FrT) */
};

/* Rule_struct keeps track of the Rules that have to be parsed and evaluated
	in the course of satisfying the --rulesfile flag. */

struct Rule_struct {
	char Name[MAX_PAT_NAME_LEN+1];      /* the name of the rule, as in RE */
   char *RuleString;   /* the verified rulestring, concatenated into one long packed
                          line */
   long Window;        /* the size of the window over which the rule must be applied
                          if not specified, then it's the length of the whole seq
                          (indicated as zero) */
};


/* very simple struct to keep track of ranges for the --clone option; since it's so small,
   could make this static for ~15 ranges, which should be about as many as any sane person
   would want ... but the world is full of insane people... */
struct Clone_struct {
  long  begin,          /* begin of this range */
        end,            /* end of range  */
        *matching_REs;  /* index of RE that matches this crit; [0] = alloc range,
                                                               [1] = 1st free el to use */
  short to_cut;         /* whether cuts are allowed in this entry range */
  char range_rule[25];  /* character string def of the rule itself */
  /* Clone[MAX_NUM_CLONE_RANGES].to_cut is set to the number of ranges we have to consider */
};


/* this is the basis for keeping track of the patterns (RE[].etc in main()) - IS dynamically allocated &
   checked so it should never run out of mem.  It has grown a bit large, but it handles all the
   patterns pretty well. */
struct RE_struct {
 /*These vars DON't change with each iteration */
   char E_nam[MAX_PAT_NAME_LEN+1],      /* RE name */
        *E_raw_sit,     /* RE name with the ' and _ and all n's included for a label, if needed.
                           also the Consensus direct from the matrix database */
        *E_wsit[2],     /* RE whole recognition site, minus _ ' trailing n's */
          E_hex[2][7];  /* the hexamer(+1 for terminator) extracted from E_wsit, used to generate the hash key */
   int E_tcut[2],    /* the cut site on the top strand, relative to the besthexindex calc'ed in BestHexWWPal()
                   		the val from REBASE is read into this var when the line is read, but is then
                     	adjusted to the correct val in BestHexWWPal() a bit later.  When the best
                    		hex was simply the 1st 6mer it was exactly the raw value read in from the
                     	REBASE, but since the best hex is calculated (chiefly for non-RE patterns),
                     	it also has to be calculated.
                        In the matrix matching it's just the lenth of the matrix/2, so it could be
                        offset considerably from the Best Hex. */
        E_olap,      /* the overlap if any in the cut - dif between the top cut and bottom cut
                        - this could be checked to provide selection on the basis of the size of
                        overhang - sometimes useful for cloning..; combine with -o? */
        E_len,       /* the length of the recognition sequence */
        E_pal,       /* whether or not the recognition sequence is a palindrome */
        E_dgen,      /* degeneracy in the recognition sequence */
        E_mag,       /* the magnitude of the degeneracy - acgt=1; n=0; yrmkws=1/2; b,d,h,v = 1/4, etc */
        WW[2],       /* Which Way - # that indicates which way Degen_Cmp has to go to do a complete comparison
                         this can also be used to indicate the orientation of the site  (?) */
        proto,       /* the index of the prototype RE entry for the list of degenerated names for ERRORs */
        Err,         /* Errors allowed in the sequence, entered from cmdline or from mod to REBASE entry */
        E_nam_l,   	/* the length of the name (E_nam, above) for calcing the  name positions in output */
        Cost;       /* units per $  - gives a # >> 1 for most, so can specify as an int if needed */
   short dcm, dam;   /* vars for dam/dcm, if needed */

/* These vars CHANGE with each iteration using a new sequence, so they have to be re-zeroed/free()'ed */
   int  E_Ncuts;    	/* the # of cuts that the RE causes in the sequence */
   int  min, Max;    /* min/Max are the PER PATTERN equiv of the -m -M flags (which refer to the TOTAL
                        number of hits), and can refer to the whole sequence or to the window defined
                        by the -W (aka --slidwin) flag if specified */
   int  E_Nfrags;  	/* the # of fragments generated by the analysis - NOT always = to Sites !! */
   long *Sites,		/* For Sites, the 0th el points to 1st free spot, so #ing starts at 1 */
        	           	/* list of hit sites for the DNA for each enzyme - was separate (Dig_Sits) */
                   	/* but for Frags, data starts in el [0].  The last el pointed to (by the 0th el)
                     	points to the last el checked in DamDcmOK(). */
        *Frags;     	/* list of fragments for the DNA for each enzyme  */
   char *Orient;   	/* will be size of (Sites/8) +1, as it's a bitarray that holds the data allows
                     	for the orientation of each pattern in parallel with Sites.  BitArray()
                      	for reasonably easy (if slightly slow) access to the info on a bit by bit basis */
   float EstSites; 	/* estimated # of sites, based on the distribution of ACGT */

/* These vars are for the Matrix Matching bit so keep them together for now */
	char *DE;         	/* Description line from ReadMatrix() */
	int *DNAMatrix[2][4]; /* To support Matrix matching: the [0] array is for the
   								FORWARD matrix; the [1] is for the REVERSE matrix - both have to be done
   								at once, of course. the [0][1][2][3] els are for the a c g t numbers.
   								Should be able to support multiple matrix searches at
   								once in a similar manner */
	int CutOff[2];		/* the cutoff value for scoring a sequence as a positive vs the matrix:
  								= the maximum value possible from the matrix multiplied by the CutOff value
  								supplied with the '-#' flag.  Optimally, this should be settable for every matrix,
  								but will use the single value for now */
	int MaxScore[2];	/* Max Possible Score for the matrix if each position gets the top score */

	int *ScoreMustBe[2];	/* array calloc'ed and calc'ed in ReadMatrix that keeps track of what the score must
  									exceed at each position so that it has a chance of succeeding by the end.  If a score
  									falls below this, there's no point in continuing. */
	int M_Acc;       	/* Matrix accession number */
};  /* 2.14.1999: sizeof(*RE) = 176 so for the standard REBASE (250 REs), it'll cost a base of 32K.
 			Populated, it'll easily double or triple that */


/* ORFs is a 2D array of type struct that keeps the orfs found in each fram separate.  The problem in
   doing it in this way (as opposed to one continuous stream of orfs) is that I have to keep track of
   which frame and keep mallocing mem as I find more orfs.  */

struct ORF_struct {
   int frame;           /* what frame it is */
   long  B_offsetAA;    /* Beginning and ending offset from aa 0 in same frame in aas */
   long  E_offsetAA;    /* Beginning and ending offset from aa 0 in same frame in aas */
   long  B_offsetBP;    /* Beginning and ending offset from aa 0 in same frame in bps */
   long  E_offsetBP;    /* Beginning and ending offset from aa 0 in same frame in bps */
   int   orflength;     /* duh... */
   float pI;            /* yup, pI */
   int sum_AAs[28];     /* array to carry the sum of each AA to print out and to calc % compo */
                        /* index by (decimal 'char' - 65) 26 + < + >  */
   char  *orf;          /* the orf itself in UPPER case, single character form */
   float MolWt;         /* in KD */
};  /* use a pointer to it and then keep mallocing more mem as needed?  I think so ...  this
                makes it 6 wide (for sep frames and then indiv ORFs can be appended to the array by
                mallocing more mem */

struct Digest_Data {   /* Replacement for Dig_Dat and associated vars, the easier to pass to functions */
    long BOS;  /* used to be BOS, now just incorporated into this stuct for compactness */
    long EOS;  /* used to be EOS, ditto */
    long *Dat; /* What used to be Dig_Dat - holds all the raw cutting/matching data.
                  It's particularly evil as the 1st 2 elements now hold Current Element Count [0] and
                  Allocated Size [1], so all writes and reads for actual data have to start at el [2]
                  [even] = Real World Coordinates (== the base # after which the DNA is cut, or the ~middle of
                           the pattern if it's not an RE, or the beginning of a regex)
                  [odd]  = RE # (== pattern index, which IS NEGATED if the pattern is a
                                NONPAL on the other strand)
               */
};

 struct Prox_struct {
    short upstr,    	/* whether N1 is upstream or downstream of N2 */
            gt;     	/* whether the dist is greater than or less than Dlo */
    int REi[2], 		/* the RE index corresponding to the names */
            matchlen; /* how long the below *matches is, for easier use by succeeding fn() */
    long    range,  	/* diff between Dlo and Dhi */
            *matches,/* variable sized array that holds the matches between N1, N2 */
            Dlo,    	/* Distance between N1,N2; if there is a Dhi, then Dlo is the lower
                     	bound on the distance */
            Dhi;    	/* If present, places an upper limit on the distance; also obviates '.gt' */
    char    st,     	/* = search type - one of 9, depending on the diff combos of flags, using it as pseudo-int */
            *optline, /* the orig option line for the -P flag that generated the results */
            *N[2];  	/* holder for Names of the factors */
};

/* PW_struct (for describing Positive Windows) is also dyn. realloc'ed to the correct size as we go */
struct PW_struct {
   int  Lo, Hi;       /* the DD->Dat[] indices that mark the SlidWin */
   long SeqLo, SeqHi; /* the RWC of the sequence that is included in the SlidWin */
   long **PatHits;     /* a dynalloc'ed array (of size NumREs) that holds the hits themselves  */
/*                                ie with 3 patterns
                       -----this is the 2nd dim in SlidWin ----->
                                 0   1    2    3    4    5
        |                    0  [ ]  where the last pos window was for each RE ([0] el is unused)*
        |                    1  [5][111][126][189][198][202]
this is 1st dim in SlidWin   2  [2][103][144]
        |                    3  [4][166][188][199][227]
        |         RE index---+   | [=======================]
        v                        |   |
                                 |   +-- where the hits actually were in Real World Coords
                                 +- how many hits there were in this window

    * the 0th row is only used for this in the 0th el of PatHits: PosWin[0].PatHits[0][0] [1] [2] ...
    	none of the succeeding els are used in this way.
*/
};


/***********************  Function Prototype Declarations  ************************/

/* fill_out_sum duplicates the degeneracy the correct # of times in the array 'sum[]' for submission
   to the hash function */
void fill_out_sum(int degen, int N_degen, int dgn_sits[]);

/* palindrome func returns 1 if the sequence is a pal, 0 if it's not  */
int palindrome(char *site, int length);

/* hash() is a central function to this program - feel free to improve it.  It takes a
(possibly degenerate) n-mer (most often a hexamer) sequence and calculates the hash value
from (for a hexamer) 0 (=aaaaaa) to 4095 (=tttttt) */
int hash(char *nmer, int dgn_sits[], int num_bases);

/* Degen_Cmp compares a degenerate DNA seq with a non degenerate sequence;
   returns 1 if they're compatible, 0 if not */
int Degen_Cmp(char *pattern, char *target, int length, int mode, int WW);

/* Rev_Compl returns the Reverse Complement of the original sequence
	if original = ggatcatttc, reverse complement = cctagtaaag */
void Rev_Compl(char *orig, char *rev, long len);

/* Anti_Par returns the anti parallel sequence of the original sequence
	if original = ggatcatttc, anti parallel = gaaatgatcc */
void Anti_Par(char *ori, char *anti, long len);

/*	Function Reverse; straight from K+R (p62) - reverses a string s in place
	if original = ggatcatttc, reverse = ctttactagg */
void Reverse(char *s);

/* Function Triplet_Reverse reverses a string triplet by triplet ie:
   ArgTrpPheAsnCys ==> CysAsnPheTrpArg  so as to make the 6 frame translations readable in the oppo
   orientation */
void Triplet_Reverse(char *str);

/* Translate translates (nondegenerate, for now) DNA sequence into protein sequence, using
    one of 8 codon preferences */
void Translate(char *DNA_in, char *Prot_out, int len, int n_letters, char Codons[8][64][7], int organism);

/* Read_Codon_Prefs reads in the file that has the codon preferences in human-readable form
   then hashes the codons and inserts them into the right place in the Codon arary.  The file
	also has a description line for each table that is read into a char array for labelling
	purposes. */
/* ARCHAIC - replaced with Read_NCBI_Codon_Data() below, but will continue to carry it for a release
	in case anyone wants to replace it for some reason */
/* void Read_Codon_Prefs(char *Codons[8][64], char *Codon_Labels[]); */

/* Read_NCBI_Codon_Data reads a slightly modified version of NCBI's long form
   of codon data:
      http://www.ncbi.nlm.nih.gov/htbin-post/Taxonomy/wprintgc?mode=t#SG4
   into a form compatible with my previous version (which was almost certainly
   derived from an earlier form of that document).  This allows frequent updating
   of the data with a minimal of editing.  */
void Read_NCBI_Codon_Data (char Codons[MAX_ORGS][N_CODONS][7], char *Codon_Labels[MAX_ORGS]);

/* SetFlags() takes argc/argv and a bunch of variables that hold all the option values and
   implements the routines to read them in and does most of the error checking.  It's related
   functionally to Interactive(), which it calls, and to CommandLine(), which makes use of
   these flags to compose a commandline argument */
char *SetFlags (int argc, char *argv[], struct SE_struct SelEnz[MAX_SEL_ENZ],
					FILE *tmpfp, char *rulestring, char datestamp[80]);

/* Usage() spits out some useful info about how to use the program if it's invoked incorrectly or
   without flags*/
void Usage(void);

/* IndexScreen(), HelpScreen#(), ScreenSep() are the help screens from
	Usage() broken apart for clarity and usability */
void IndexScreen(void);
void HelpScreen1(void);
void HelpScreen2(void);
void HelpScreen3(void);
void HelpScreen4(void);
void HelpScreen5(void);
void HelpScreen6(void);
void HelpScreen7(void);
void ScreenSep(void);

/* compare is a dippy little function that qsort needs to perform its sordid little sort
   abscompare does the same thing but sorts in order of abs magnitude instead  */
int compare(const void *n1, const void *n2);
int abscompare(const void *n1, const void *n2 );

/* GetSequence reads and formats the sequence from **stdin** and returns the address of the read-in
   sequence, as well as length of the sequence, along with the bracketing repeats - needs only a
   few extra variables to do it and it cleans up main considerably */
char *GetSequence(long *tot_seq_Cnt, long *seq_len, long *Degen_Log[]);

/* GetSequence2 reads, filters, and formats the sequence handed to it from the SEQIO pkg fn
	seqfgetseq(), which supposedly takes care of all format conversions, malloc'ing mem, etc.  It may
	also filter the sequence correctly, but I'm not counting on that yet, so there's still quite a bit
	of error-checking in GetSequence2 */
char *GetSequence2(char *SEQseq, long *tot_seq_Cnt, long *SEQ_len, long *Degen_Log[]);

/* PrintSitesFrags prints out the sites and or fragments either sorted or unsorted, depending on
	whether the array has been sorted - it doesn't care, it just prints whatever's in the array
	nicely */
void PrintSitesFrags(int NumProtos, int reps, int *Protos, int sites, long seq_len);

/* PrintGelLadderMap does what it implies - prints a Gel and/or Ladder map sort of like the GCG
	program, but in textmode only for now.  It uses a preset (#defined) output width of ~200 chars, vs
	the regular output which is settable via a flag.  This being the case, it writes to a user-
	selectable file, rather than to stdout.  The output is meant to be processed by an postscript
	wrapper program like genscript, so that the output can be printed landscape in small font to be
	usable.  This is one of my mulitvalent functions that wil process either Sites or Fragments
	depending on a mode switch (gel) and the data fed it - perhaps not such a great idea...  In this
	fn(), if it's called with gel=1, it will want Fragment data (RE[].Frags) If it's called with gel=0
	(ie make a ladder) it will want Sites data */
void PrintGelLadderMap(int NumGREs, long seq_len, int *GREs, int gel);

/* SearchPaths () takes a filename and examines various environment variables to locate that
	filename.  If it can find the file, it returns a pointer to the full path name; if not, it
	returns NULL */

char *SearchPaths(char *InputFileName, char *FileType4Err);

/* realhash() is a real hashing function, straight from Sedgewick (p233) that returns a good hash
	value within the scope of the allocated space, because of the mod() tablesize MUST be prime for
	reliable results */
unsigned realhash(char *string2hash, int tablesize);

/* MatchSelREs() takes RE names and (typically,m if there was a hash collision, determines if it was
a real match or just a random collision; returns a 0 if no match, 1 if a true match */
int MatchSelREs(char *suspect, struct SE_struct SelEnz[MAX_SEL_ENZ], int NumSelREs);

/* DownCase() downcases all the letters in the submitted string  */
char *DownCase(char *string);

/* ReadREs does just that - given a filename (EnzFileName), it tries to open it and if successful, reads all the
   uncommented REs into the global RE[] struct, and then returns the # of entries in RE (NumREs).
   Matching against the SELECTED REs specified on the command line takes place after this functions
   returns.  Besides some random declarations that have to be cleaned out, doesn't deal with the
   GRE arrays or indices - they are restricted to main(), from which this code chunk was was excavated.
   It returns a pointer to hashtable which is an important array for the Coutting routines */
int ReadEnz(char *EnzFileName, FILE *tmpfp, int *NumREs, struct SE_struct SelEnz[MAX_SEL_ENZ],
			int *hashtable[4096], long SeqLen);

/* DegenCut() is a functionization of the cutting routines needed to do the cutting analysis on
	degenerate (gyrwsctkmyycg) sequence - also uses HorribleAccounting() as a function.  This fn()
	includes routines for both 'lite' and 'heavy' degenerate cutting, where lite is defined as ignoring
	any key hexamer which has ANY degeneracy and 'heavy' will force degenerate matching of any
	degeneracies up to a compiled in limit (currently 256 = 4 ns (or equivs) out of the keying hexamer)
*/
void Cutting(char *sequence, int *hashtable[4096]);

/* HorribleAccounting() is a function that became one simply because it contains common code for the
	cutting functions mentioned above.  If it turns out to be about as fast as the inline code, I'll
	functionize it in NonDegenCut() as well - OK , did it */
void HorribleAccounting(char *sequence, unsigned int Key, int mode, int *hashtable[4096], long CSP/* , int already[]*/);

/* ORF_Analysis() does what it implies - analyzes the ORFs found - have to hand it the translated sequence
   (*prot), it's length (protlen), anf the frame (thisframe) and it finds and calculates the stats on all
   the ORFs above the cutoff size (F.ORFs) */
void ORF_Analysis(char *prot, long protlen, short thisframe);

/* ORF_Calcs() does the actual calculations for the analyses; so far just the molwt calc
   'ORFNum' is the current orf number and 'f' is the frame  */
void ORF_Calcs(int ORFNum, int f);

/* little fn() that farts the ORF data to stdout - only reason to keep it separate from ORF_Analysis
	is that I'll probably be adding the ORF Map function to ORF_Analysis and this will keep it a tad
	more modular */
void PrintORFs(int ORFNum, int frame, int MinORF);

/* Degen_Calc() is a sleazy little fn() to calculate the degeneracy of the hexamer (possibly later,
	any string) without the added overhead of the other stuff in hash() that also does this */
int Degen_Calc(char *s);

/* this little sucker returns the 'magnitude' of the string given it as per the same parameters as
	hash() */
float MagCalc(char *s, int len);

/* BadMem() just writes a message and calls exit() */
void BadMem(char* message, int exitcode);

/* Proximity() does all the proximity matching for the patterns detected in the search phase.  Takes
    care of it's own output */
void Proximity(void);

/* BadFlag just warns the user of a badly constructed flag, points to where more info
	can be got, and exits */
void BadFlag(char *flag);

/* CalcEstSites() calculates the probablility of the number of hits of the sites based
	on the base distribution of the sequence that's being analyzed */
float CalcEstSites(char *BindSite, int BSLen, float FrACGT[4]);

/* PadEndForEnd() performs the "pad the beginning sequence with the end and the end with the
	beginning" regardless of whether it's the whole sequence or a subsequence, using the vars 'begin'
	and 'end'.  It is a functionized stanza from GetSequence, but the code in GetSequence has been left
	in place because this funtion has to do some memory manipulations that are not needed in
	GetSequence, due to the way that SEQIO sets up the buffers */
char *PadEndForEnd(char *SeqIn, long *Pad_Seq_Len, long *Bare_Seq_Len);

/* Dumps an error about the inappropriate use of the '-w1' (very wide) flag */
void Errw1(void);

/* Hookey() allows the use of 'tagging' a single RE in a combination of REs (typically a
	6 or 8 cutter, paired with a 4 cutter so that you can pull out those fragments that have
	only an end cut by the tagged RE.  John Hookey suggested this as an assist for some genome
	mapping work he was doing.  Remains to be seen if all this work will lead to anything..but
	I've taken it this far.. */
void Hookey(int H, int AI, long seq_len);

/* DumpDataForPlot() simply dumps data that is stored in internal arrays to stdout so
    that it can be manipulated by external plotting and analytical apps - takes care of all
    output as well as the array manipulations that are required to flip the internal arrays
    on the diagonal or make them long and skinny (as gnuplot prefers)  */
void DumpDataForPlot(long seq_len,  int ProtoCutters,  int *Protos);

/* MatrixMatch(sequence) takes the sequence and matches all the matrices that have been read in
	and fills in DD in the same way that Cutting() does, so that all later functions continue to work
	as before.  More details later */
void MatrixMatch(char *sequence, int NumREs);

/* ReadMatrix() reads in the matrix file in Transfac format, into the same RE
   struct that ReadEnz() does (but filling in a few more things) and allows the
   matching to proceed to MatchMatrix and then reporting much the same as before */
int ReadMatrix(char *EnzFileName, int *NumREs, struct SE_struct SelEnz[MAX_SEL_ENZ]);

/* imax return the int max of 2 values (mainly to save space) - it's another
	1 liner, but should be used relatively infrequently */
int imax(int v1, int v2);

/* BitArray() is a 1st pass approx to doing interesting things with bit strings
   to store bit-level amounts of info.
   bstr = the pre-alloc'ed bit string
   I    = the index to the bit under consideration in normal terms
            ie. if there was a bit array bitarray[53], then I would be 53 and
            that would mean the 53rd bit in that array (or 6 Bytes and 5 bits into
            a char array).
   mode = 0 clears the bit; returns 1 for success, 0 for failure
          1 sets the bit; returns 1 for success, 0 for failure
          2 querys the bit; returns 1 if it's set, 0 if not  */
int BitArray(char *bstr, int I, int mode);

/* ExtractSeqAtHit() takes the list of patterns passed to teh program via the '-r' flag and extracts
	the sequences around the hits beginning 'b' bases 5' of the pattern midpoint and ending 'e' bases 3'
	of the pattern midpoint, where 'b' and 'e' are the parameters of the --take b,e,revcomp flag.
	('revcomp' is the variable that indicates whether the sequence should be revcomplemented before
	being sent to stdout).  the only thing that it has to be passed is the number of REs to chew thru
*/
void ExtractSeqAtHit(SEQINFO *SEQinfo, int SEQlen, char *sequence, int ProtoCutters, int *Protos);

/* lil nothing fn() to check that DD->Dat isn't going to blow */
int DD_Mem_Check(void);

/* This puppy takes a pseudo regex expression in DNA IUPAC codes and translates it into REAL
	regex form for exapmple: gcnytvmtk -> gc.[ct]t[^t][ca]t[gt]
	the main problem seems to be getting the expression to this function without being creamed by the
	shell or other intermediate routines that are trying rabidly to interpret the characters as
	meaning something.  This is what they're paid to do, but it doesn't make things easy for the user
	or programmer.  This fn() assumes that the entering DNA string has been **downcased and
	pre-filtered** to remove any NON-IUAC characters from the string - maybe a faulty assumption.  */
void DNA_IUPACtoRegex(char *Regex, char *IUPAC);

/* ReadRegex reads a file of regular expressions in the form:
comments
comments
..
;this line is commented out
RegexName     gg(gc|tty)rma{6,7}ccgtackckkyav   ! If there's a '!', stuff after it is ignored for now
OtherName     gwtcg{2,4}ta(aa|tg)tayng          ! After the '!' various other fields
;Name1        gaan{5,7}ttyy(am|ma)tggtyk        !

It does some primitive munging and error checking and then stuffs what it can into RE[].
If the opt string is 'FILE' alone, it assumes the file name 'regex.data'.  */
int ReadRegex(char *EnzFileName, FILE *tmpfp, int *NumREs, struct SE_struct SelEnz[MAX_SEL_ENZ]);

/* RegexMatch() takes a compiled regex pattern and searches both strands for it (actually searches
	the same strand for the regex, then Anti_Par()s the sequence and searches it again (which I have to
	make a function to do..)  */
int RegexMatch(long SEQlen, char *sequence, int NumREs);

/* ReverseTranslate() is another eponymous fn(), taking a string of single letter amino acids, it
   munges them in to reverse translated DNA, initially towards contributing to '--silent' (the
   fn() that looks for SILENT SITES, RE sites that can be introduced by mutation that will not
   change the coding sequence.
   WARNING: In doing this, there are some translational ambiguities produced.  It depends on which
   Codon table is being used, but for the Standard one, Arg, Leu, and Ser are reverse translated to
   mgn, ytn, and wsn respectively, which will not forward translate unambiguously, so it remains
   to the user to be careful about interpreting the results.  If I can come up with a plan to
   verify the sequence once this is done, I'll implement it. */
char *ReverseTranslate(char *Prot, char Codons[MAX_ORGS][N_CODONS][7], int OrgTable);

/* BestHexWWPal() does much of the calculation to set the RE values related to which way the extended check
   should look (WW), what the best hexamer is, what the center of the site is, etc.  */
void BestHexWWPal(char *RE_rawsite, int len, int Cur);

/* tacg_SlidingWindow() searches thru DD->Dat, and verifies whether the
   per-pattern min/Max limits are being met within that window.  It  returns the
   number of Positive Sliding Windows over all the patterns searched.   The
   struct that holds all the data is global, so it can be accessed from anywhere
   It also handles it's own printing if requested, by setting PrintIf1 = 1  */
int tacg_SlidingWindow(int rule, int NumREs, int PrintIf1, struct SE_struct SelEnz[MAX_SEL_ENZ], int *Protos, int NumProtos);


/* tacg_StoreHits() functionizes some repeated accounting code, copying the RWC of the hits
   found in the sliding window of DD->Dat[] to the PosWin[PWCnt].PatHits[re][] struct for later
   passing to whatever code wants them.  I could have just stored pointers TO the orig DD->Dat
   location, but this wouldn't have saved any space and would have made later referencing the
   data that much arcane */
void tacg_StoreHits(int Hi, int Lo, int re, int PWCnt);

/* tacg_SetF (usually stored in tacg.c) just sets all the default vars in the flags
   struct F */
void tacg_SetF(void);

/* if either the --dam or --dcm flags have been used, DamDcmOK() checks to see
	if the REi under consideration should be printed - has it been masked or is
   it OK to print it.  What it does: passed the RE index of the enzyme, it
   checks to see if that RE is possibly affected and if so, whether the 'raw'
   hit recorded in DD.Dat is present in RE[].Sites.  The index to the
   last-checked position in RE[].Sites is kept in:
   RE[].Sites[RE[].Sites[0]]       as 'RE[].Sites[0]' points to the end of the array.
   In this way, only have to incr by 1 to get to the next check point.
   RWC = Real World Coordinates
   REi = RE index
*/
int DamDcmOK(long RWC, int REi);

/* Add_dam and Add_dcm just fill the RE[1] and RE[2] entries for these special
   methylation enzymes. */

void Add_dam(int NumREs);
void Add_dcm(int NumREs);

/* this just genrates the HTML TOC entries, based on the flags for major output
	sections that it reads from F.___  */
void GenerateTOC(char datestamp[80]);

/* Quoted direct from Gray Watson, who wrote this function.  It's defined here
   to prevent compilers from complaining about the fn() def being missing:
   This is a function which should be in libc in every Unix.  Grumble.
   It basically replaces the strtok function because it is reentrant.
   This tokenizes a string by returning the next token in a string and
   punching a \0 on the first delimiter character past the token.  The
   difference from strtok is that you pass in the address of a string
   pointer which will be shifted allong the buffer being processed.
   With strtok you passed in a 0L for subsequant calls.  Yeach.
   (See strsep.c for more info and more functions).
   hmm - it is now part of the Linux stdlib - defined in <string.h>, so the
   below declaration is no longer needed, but it still generates multiple
   errors of the sort:
   warning: assignment makes pointer from integer without a cast
   (under gcc 2.95 anyway).

*/
/* char *strsep(char **string_p, char *delim);  */


/* Following functions that begin with 'Flags_case' are functionized
	case stanzas from SetFlags() to take large chunks out of the main part of the
   sub.  None of the ones remaining should be more than a screen's worth of code
*/
char *Flags_case_r(char *optarg, FILE *tmpfp, char datestamp[80]);
void Flags_case_p(char *optarg, FILE *tmpfp, char datestamp[80]);
void Flags_case_15(char *optarg, char *rulestring, struct SE_struct SelEnz[MAX_SEL_ENZ], int *SECnt, char datestamp[80]);
void Flags_case_P(char *optarg, struct SE_struct SelEnz[MAX_SEL_ENZ], int *SECnt);
void Flags_case_x(char *optarg, struct SE_struct SelEnz[MAX_SEL_ENZ], int *SECnt);
void Flags_case_T(char *optarg);
void Flags_case_X(char *optarg);
void Flags_case_O(char *optarg);

/* GrafMap is the 1st take on the long-delayed graphical plasmid/linear maps function */
/* on calling, runs thru RE and makes a 1st pass guestimate of how to map all the hits. If there
	are too many hits, it should truncate at a level which does not leave a pattern partially
   mapped.   Drawing ORFs other features will be in the 2nd pass after I have some graphics
   working and have a better idea of what kinds of computation will be nec.  Needs to open a
   file to write postscript to, may need to call for fork ghostscript as a system level call
*/

void GrafMap(long *Degen_Log, int *Protos, char *Description, char datestamp[80], float FrACGT[4]);


/* LinearMap() is another functionization of a chunk of code that is modular
	and functional enough to makeinto its own fn().  It takes the DD.Dat struct
   and parses the hit data out of it and composes a Linear Map out of it. */
void LinearMap(long seq_len, int basesPerLine, int NumREs, char *sequence,
 					int max_okline, char Codons[MAX_ORGS][N_CODONS][7], int codon_Table, int *Protos);

/* Clone() goes thru RE[] and figures out which REs match the conditions of the option flag, where
   --clone '#_#,#x#,#x#,#_#' (up to MAX_NUM_CLONE_RANGES).  '_' indicates that the range cannot be cut;
   'x' indicates that it can be cut.
*/
void Clone(int ProtoCutters, int *Protos);

/* realloc_Clone_matching_REs() does what it says - just saves a bit of code.  */
void realloc_Clone_matching_REs(long current, long index);

/* RE_Filter returns 1 if the RE indexed by Proto[REi] matches the tests for
	Magnitude
   Number of Cuts
   Cost
   Overlap
   etc that are can be used to fine-tune the selection of REs */
int RE_Filters_OK(int REi);

/* tacg_Eval_Rules_from_File() reads the rulestrings from a file and evaluates
	them, using the tacg_SlidingWindow() that was originally designed to test
   one rule at a time */
/* returns an int (SEi); the number of different subpatterns that are represented in the ENTIRE file
   this is still going to have to use the realhash() to check to see which are the same, tho...  */

int tacg_Eval_Rules_from_File(struct SE_struct SelEnz[MAX_SEL_ENZ], char datestamp[80] );
                             /* long ProtoNameHash[], int NumREs,    */


/* stripwhitespace does what it says.  it strips whitespace (' ', \t, \n, \r) from the passed string and
	operates ON THAT STRING, (ie IT IS DESTRUCTIVE) from the requested begin point to the requested end point
   testing the begin and end points to see if they are within logical limits.  It will NOT detect an already
   bad string.  If end is 0 it will go to the end of the string. It returns the number of whitespace characters
   replaced. */

long stripwhitespace(char *string, long begin, long end);

/* tacg_Fill_SelEnz_fr_RuleString() takes the SelEnz sruct and the rulestring (and a few other housekeeping vars
   and parses the rulestring passed to it to populate SelEnz and do a namecheck to make sure that the names being
   parsed out of the rulestring haven't been seen before and aren't otherwise objectionable
   takes as args:
       SelEnz,
       the starting SEi (as an int * as it's changed in the fn())
       the rulestring (copied in the function so that the original rulestring is untouched)
       NameCheck (has to survive fn() to check for multiple identicals over multiple rulestrings,
          as well as to check for multiple identicals in the same rulestring (A&B|A&C)
       ctemp1 - Name of the Rule (for error messages)
    and back comes
       the (more) populated SelEnz
       and the ending SEi.  */

void tacg_Fill_SelEnz_fr_RuleString(struct SE_struct SelEnz[MAX_SEL_ENZ], int *SEi,
            char *rulestring, short NameCheck[RE_NAME_TAB_SIZ],  char *ctemp1);

/* passed the pointer of the populated struct, it calculates the x,y coords of the next hit and passes them
   back as struct elements xn,yn */
void Calc_Next_Hit_Coords (struct NextHit_struct *NH);

   /*tacg_Draw_Seq_Tics() draws tics on the plasmid map
 		int tic_size - the length of the tics in pts
      float line_width - the width of the tic line
      long seq_len - the sequence length in base pairs
      struct NextHit_struct *NH - the Next Hit struct
      int drawlabels - whether to print labels or just write the tics
      long repeat_period - the tic period (1000 -> 1000, 2000, 3000, etc)
      int skip_period - what period should be skipped.  If already done repeat_period -> 1000 and
            doing minor tics at 100, don't want to overwrite the 1000 label. set to -1 to ignore skip_period
      int font_string - the font characteristics as Postscript likes it:

      FILE *PS - the file pointer to the file we're writing to
   */
void tacg_Draw_Seq_Tics(int tic_size, float line_width, long seq_len, struct NextHit_struct *NH, int drawlabels,
                        long repeat_period, int skip_period, char *font_string, FILE *PS );

/********************** Function Complement_In_Place  ********************************
*   Complement_In_Place takes a char pointer to the sequence to be converted and converts it
*   in-place, for the distance passed by 'len'.  Thus you can operate on sequences and
*   subsequences (by passing in the pointer to an offset and the length you want to
*   operate on) - if original = aggctgctggat,  Complement = tccgacgaccta
**************************************************************************************/
void Complement_In_Place(char *ori, long len);

/*
 * This pitiful fn() just jams the psprolog header into a string and returns it to
 * whatever calls it so that we can skip all the file manipulation and problems when the
 * file isn't in the right place  when we expect it.
 */
char *tacg_PS_Prolog(void);

void tacg_psladder_header(FILE *PS, long flabels[15], double ftics[15], double pt_for_seq, char *Description, char datestamp[80]);

void tacg_psladder(int *Protos, char *Description,  char datestamp[80]);

/* NmerAnalysis() takes the sequence (and the antipar sequence?) and returns the number of
    degens enountered.  Or it could return a pointer to a struct or array that actually holds
    the data.
*/
long NmerAnalysis(char *sequence, long primes[]);

/*************************  Function nmer_hash  **************************************************
*   nmer_hash() takes a pointer to the n-mer string and generates the                            *
*   integer equivalent.  If there are degens in the string, it will exit with an                 *
*   error code of 0; otherwise it returns the long equiv of the string.  If the verbose option   *
*   is set it will tell you what it bonked on.                                                   *
*************************************************************************************************/
unsigned long int nmer_hash(char *nmer);

// nmer_is_not_clean(nmer) does a scan of the nmer for any degen
// returns 1 if nmer contains a degen, 0 if it's clean, so true if its dirty
short nmer_is_not_clean(char *nmer, int nmer_len);


//Prints both an ORF map and a MET / STOP map
void PrintORFMap(long seq_len, int *metstop[]);

// tacg_TimeOps just runs thru a lot of ops to time them out for interested / bored developers
double tacg_TimeOps (void);

// the debug message emitter tha should have been written a long time ago.
void tacg_debug(int dbg_lvel, char *fle_nm, int lin_nbr, char *err_msg);
