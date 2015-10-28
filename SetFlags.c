/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: SetFlags.c,v 1.5 2005/01/26 22:44:08 mangalam Exp $  */

/* The use of this software (except that by James Knight, which is described in
   'seqio.c', and by Gray Watson, which is described in strsep.c) is bound by
   the notice that appears in the file 'tacg.h' which should accompany this
   file.  In the event that 'tacg.h' is not bundled with this file, please
   contact the author.
*/

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include "tacg.h" /* contains all the defines, includes, function prototypes for both main() and functions */

/* SetFlags() takes argc/argv and a pointer to SelEnz (struct that holds the names/other data of the
	Enz's selected from the command line and implements the routines to read them in and do the
	error checking.  It returns a pointer to the name of the Enzyme file name.  Yes, I know there
	are far too many options....  */

char *SetFlags (int argc, char *argv[], struct SE_struct SelEnz[MAX_SEL_ENZ],
                FILE *tmpfp, char *rulestring, char datestamp[80])
{

    /* the compiled_regex is generated in RegexMatch, based on the values stored in RE,
    	either from here (from single or multiple invocations of '-x') or from reading the
    	file of regex patterns in ReadRegex, which means that it's once more thing we DON'T
    	have to keep passing around - it ONLY needs to be declared in RegexMatch */

    char ch='.';
    char *tok, *tstr, *SeqFile;
    char *temp = "             "; /* temp holder for RE name in the selection string */
    char *clonerange = "                         ";
    char *RebFile = "\0\0";    /* holder for alternative Rebase file string */
    char *RuleFileName;        /* holder for Rule file string */
    long kk;
    int   i, j, k, len, badranges=0,
                        SECnt=0, /* Selected Enzyme Counter for SelEnz */
//         pCnt=0,   for patterns
                        clone_rule=0, /* counter for the number of ranges for the Clone[] struct */
                        option_index = 0, /* var for getopt routines */
                        c;    /* holder for the 1st char past the '-' (flag indicator) - has to be of type int because of getopt */

    /* following struct is for the getopt_long options that will have to be used, now that I'm starting
    	to exhaust the single letter flags */
    static struct option long_options[] = {
        /*              + has an arg (0-no arg, 1-required, 2-opt)
            argname     |  + specifies  how  results  are  returned
               |        |  |   + the value to return, or to load into  the variable to index the switch    */
        {"idonly",    1, 0, 'i'},  /*  --idonly (-i) flag */
        {"extract",   1, 0, 'X'},  /*  --extract b,e,revcomp X=extract*/
        {"help",      0, 0, 'h'},  /*  --help  for gnu compliance */
        {"overlap",   1, 0, 'o'},  /*  --overlap = -o  */
        {"regex",     1, 0, 'r'},  /*  --regex = -r    */
        {"silent",    0, 0,  1 },  /*  --silent flag (SILENT SITES)*/
        {"strider",   0, 0, 'c'},  /*  --strider sorting = -c */
        {"HTML",      2, 0, 'H'},  /*  --HTML = -H {0 | 1} */
        {"summary",   0, 0, 's'},  /*  --summary = -s  */
        {"version",   0, 0, 'v'},  /*  --version  for gnu compliance */
        {"example",   1, 0, 12 },  /*  --example for the 'example' function to show how to write one */
        {"slidwin",   1, 0, 'W'},  /*  --slidwin #bases to define the sliding window over which certain stats are gathered */
        {"rule",      1, 0, 15 },  /*  --rule - formal rules for motif logic over whole seq or slidwin */
        {"dam",       0, 0,  4 },  /*  --dam -> simulation with dam methylation (GmATC) */
        {"dcm",       0, 0,  5 },  /*  --dcm -> simulation with dcm methylation (CmC(A/T)GG) */
        {"cost",      1, 0,  6 },  /*  --cost -> in Units/$ */
        {"strands",   1, 0,  7 },  /*  --strands {1|2} -> 1 or 2 strands in linear map (F.Strands)*/
        {"notics",    0, 0,  8 },  /*  --notics  -> no tics in linear map => F.Tics = 0 */
        {"raw",       0, 0,  9 },  /*  --raw  -> take sequence in as raw sequence, as in version 2; use old GetSequence() */
        {"ps",        0, 0, 16 },  /*  --ps   -> takes topology from -f  */
        {"numstart",  1, 0, 11 },  /*  --numstart # -> what number to start numbering raw sequence at */
        {"logdegens", 0, 0,  2 },  /*  --logdegens -> forces logging of alldegeneracies for graphics output */
        {"pdf",       0, 0,  3 },  /*  --pdf -> forks gs to convert the ps file to pdf */
        {"tmppath",   1, 0, 17 },  /*  --temppath used to pass a temp path to tacg from a cgi  */
        {"clone",     1, 0, 13 },  /*  --clone '#_#,#x#,...' ranges separated by '_' cannot be cut; 'x' can be cut */
        {"rulefile",  1, 0, 14 },  /*  --rulefile /path/to/rulefile */
        /* following use optional args format specifier '2' bc no other way to separate the numbers in the optstring */
        {"rev",       2, 0, 18 },  /*  --rev reverses sequence in-place before running the operation */
        {"comp",      2, 0, 19 },  /*  --comp complements sequence in-place before running the operation */
        {"revcomp",   2, 0, 20 },  /*  --revcomp reverse-complements sequence in-place before running the operation */
        {"revnum",    2, 0, 22 },  /*  --revnum causes reverse numbering of the sequence in output,
                                        usually used with --revcomp - option handling is in place, but
                                        the actual code is not currently implemented - would go
                                        into LinearMap() at about line 1208 */
        {"psladder",  2, 0, 23 },  /*  --psladder generates the postscript version of the ladder map.
                                        has optional flag 'p' (portrait) or l (landscape) for the map */
        {"infile",    1, 0, 24},   /*  --infile specifies a file of input sequence to be processed.
                                    if there's no --infile specified, tacg assumes that input is
                                    going to be passed in via stdin, same as it ever was */
        {"nmer",      1, 0, 25},   /*  --nmer # sets the size of the hash word, calls the nmer analysis */
        {"orfmap",    0, 0,  26},  /*  --orfmap - a ladder map for ORFs; uses frames from -O. */
        {"timeops",   0, 0,  27},  /*  --timeops - time different CPU ops. */

        {0,           0, 0,  0 }   /* this last element HAS to be filled with 0's */
    };

    /* in args below, single letters indicate single letter opts; those followed by colons have
    	further args some of which have to be parsed by strsep() and stuffed into the appro array.
    IMPORTANT!!  the optstring below also needs to agree with the "static struct option
    long_options[]" above. see `man 3 getopt` for descriptions to the fields
    EXAMPLE NB: note that the example option is followed by a ':', indicating that it needs an
    option string.  If it doesn't get one, getopt_long() will emit a warning */

    while ((c = getopt_long(argc, argv, "1234589chlLsS::v6:7:10:11:12:13:14:15:24:25:16::17:18::19::20::22::b:23::26::C:D:e:i:f:F:g:G:H::m:M:n:o:O:p:P:t:T:r:R:V:w:W:x:X:#:",
                            long_options, &option_index)) != EOF) {

        switch (c) { /* switch based on the letter after the '-' */

        default:
            fprintf(stderr, "\n\"%c\" - Unrecognized or Unimplemented Flag!\n", c);
            break;

        case 27:
            if (optarg != 0) BadFlag("--timeops"); /* die on bad flag */
            F.Timeops = 1;
            F.Output = 1;
            break;

        case 26:
            if (optarg != 0) BadFlag("--orfmap"); /* die on bad flag */
            F.Output = 1; /* set F.Output to indicate that some outout has been requested */
            F.Orfmap = 1;

            break;
        case 'l': /* include the ladder map; currently a text/graphic, set with -w, same as everything else.  */
            if (F.Width > MAX_BASES_PER_LINE) Errw1();	/* check for -w1 condition */
            F.Output = 1; /* set F.Output to indicate that some outout has been requested */
            F.Ladder = 1;
            break;


        case 25: /* --nmer # to call the nmer analysis & set the hash word size. */
            if (optarg == 0) BadFlag("--nmer"); /* die on bad flag */
            j = atoi(strsep(&optarg, ",:."));
            if (j>1 && j<21) {
                F.nmer = j;
                F.Output = 1; //to allow analysis to go forward.
            } else  {
                fprintf (stderr, "\n'--nmer' flag must be > 1 and < 21.  Try again\n");
                exit(1);
            }
            break;



        case 24:  /* --infile to specify the input seq file. */

            if ((F.infile = (char *) calloc((size_t)(strlen(optarg)+2), sizeof(char))) == NULL) {
                BadMem("F.infile", 1);
            }  // get some appro amount of mem
            if (strlen(optarg) > 128) BadFlag("--infile "); /* die on bad flag */
            if ((SeqFile = SearchPaths(optarg, "--input seq file")) == NULL) {
                fprintf(stderr,"Can't find the alternative input seq file!! - Check spelling, ownership, existance, etc - Bye!!\n");
                exit(1);
            } else {
                strcpy(F.infile, SeqFile);
            }
            break;

        case 23: /* --psladder {p|l} */
            j = 1;
            if (optarg == 0) { /* if no flag, assume portrait format */
                fprintf(stderr, "here we are at case 23\n");
                j = 0;
                F.psladder = 1; /* portrait */
            } else { /* there is a flag - what is it? */
                if (strlen(optarg) > 1) BadFlag("--psladder "); /* die on bad flag */
                else {
                    ch = optarg[0];     /* strncpy(ch, optarg, 1); */
                }

                if (ch == 'p') F.psladder =  1; /* portrait */
                else if (ch == 'l') F.psladder =  -1; /* landscape */
                else BadFlag("--psladder");
            }
            break;

            /* opts to do sequence flipping, flapping, folding, forking  */
        case 18: /* --rev */
            if (optarg != 0) BadFlag("--rev "); /* die on bad flag */
            F.Rev = 1;	/* just set the flag */
            break;

        case 19: /* --comp */
            if (optarg != 0) BadFlag("--comp "); /* die on bad flag */
            F.Comp = 1;	/* just set the flag */
            break;

        case 20: /* --revcomp */
            if (optarg != 0) BadFlag("--revcomp "); /* die on bad flag */
            F.RevComp = 1;	/* just set the flag */
            break;

        case 22: /* --revnum  to reverse the order of sequence numbering on SOME output formats. */
            fprintf(stderr, "\n--revnum has not been entirely implemented yet - sorry.  Write me if you want it\n\n");
            exit(1);
            /* rest of --revnum is non-operational for now */
            if (optarg != 0) BadFlag("--revnum"); /* die on bad flag */
            F.RevNum = 1;	/* just set the flag */
            break;

            /* --tmppath '/path/to/tmp/directory' */
        case 17:
            if (optarg == 0) BadFlag("--tmppath"); /* die on bad flag */
            /* assume path is OK for now, but have to stat it eventually */
            strcpy(F.tmppath, optarg);
            break;

            /* case 3 = fork ghostscript to convert the ps file to a pdf file. */
        case 3:
            F.PDF = 1;
            F.GrafMap = 0;
            F.Output = 1;
            F.Topo = 0; /* plasmid map implies a circular plasmid, no? */
            F.LogDegens = 1;
            break;

            /* case 14 = --rulefile specify the rulefile to pull predefined rules to search for.
            If not specified, assume it's called rules.data and is in the std places.
            --rulefile has to do the whole thing - verify the file is there or find out where it is,
            read the file completely, stash all the bits in the Rules struct and at the same time,
            put all the subpattern names into SelEnz so that it can be passed to ReadEnz() to
            filter the names of the subpatterns coming in.  When all that is done, all that's left is
            to run thru the Rules struct, iteratively passing each set to tacg_SlidingWindow */

        case 14:
            if (F.SelEnzUsed != 1) { /* don't tromp on other fn()s that might want to use SelEnz */
                if (optarg == 0) BadFlag("--rulefile"); /* die on bad flag */
                if ((RuleFileName = SearchPaths(optarg, "Rule File")) == NULL) {
                    fprintf(stderr,"\nCan't find the alternative Rule file you've specified!!\n"
                            "Check spelling, ownership, existance, etc.  As a last effort, searching for\n"
                            "'rules.data' in the std places\n");
                    if ((RuleFileName = SearchPaths("rules.data", "Rule File")) == NULL) {
                        fprintf(stderr,"\nCan't find the standard Rule file either - check more carefully - BYE!\n");
                        exit(1);
                    }
                } else {   /* stuff it in F.RuleFile */
                    F.RuleFile =calloc(strlen(RuleFileName)+1, sizeof(char));
                    strcpy(F.RuleFile, RuleFileName); /* pointer assigned to another pointer to carry */
                    F.Rules = 2; /* notify via F that we want to examine Rules analysis via rules from file (1=from commandline) */
                }

                SECnt = F.NumSubPatts = tacg_Eval_Rules_from_File(SelEnz,  datestamp);
                /* this will have to set some warning flags to warn off anything else that uses SelEnz  */
                F.SelEnzUsed = 1; /*  to warn off other fn()s that might want to use SelEnz */
                F.Output = 1;
            } else {
                fprintf(stderr,"\nYou requested the '--rulefile' option but you've also requested a conflicting option\n"
                        "such as -x, -P, -p, and of course --rule. Please edit the command & try again.\n");
            }
            break;


            /* case 13 -  the range specs should be changed from:
                #x# = CAN hit here
                   to
               #x# = MUST hit here
               #_# = CANNOT hit here (as before)
               DNA not mentioned in the range spec CAN be hit or not, but there should be a pref for fewer hits..? */

        case 13: /* --clone '#_#,#x#,...' ranges separated by '_' cannot be cut; 'x' can be cut */
            F.Clone = 1;
            F.Output = 1; /* allow this to qualify as 'Output' */
            if (optarg == 0) BadFlag("--clone"); /* die on bad flag */
            /* break optstring into ranges of #s */

            /* copy 'optarg => rulestring' to be passed back */
            strcpy(rulestring, optarg);
            if (strlen(optarg) < 3) { /* check if there's something in optarg */
                fprintf(stderr, "!! '--clone' flag option is too short or missing.  The '--clone' flag usage is:\n"
                        "'--clone '#x#,#_#,#_#..' (Maximum of %d ranges), where:\n"
                        "#_# = a range in which there CAN'T BE cuts\n"
                        "#x# = a range in which there MUST BE cuts\n"
                        "where the '#' is an integer indicating the position in the sequence from the beginning\n",
                        MAX_NUM_CLONE_RANGES);
                exit(1);
            }
            // Have to alloc the space for Clone_S
            if ((Clone_S = calloc(INIT_CLONE_S,sizeof(*Clone_S))) == NULL) {
                BadMem("Clone_S fails malloc", 1);
            }
            clonerange = strsep(&optarg, ","); 	/* strsep replaces strtok, returns 1st range into clonerange */
            while (clonerange != 0 && clone_rule < MAX_NUM_CLONE_RANGES) { /* break optarg  into ranges & then load into Clone_S struct */
                // alloc Clone_S.matching_REs at each step - give the matching_REs a starting value of 20
                if ((Clone_S[clone_rule].matching_REs = calloc(20,sizeof(long))) == NULL) {
                    BadMem("Clone_S.matchingREs fails malloc", 1);
                }
                strcpy(Clone_S[clone_rule].range_rule, clonerange);
                Clone_S[clone_rule].to_cut = -1;
                /* figure out if the range is to_cut or NOT_to_cut */
                j = (int)strlen(clonerange);
                for (i=0; i<j; i++) {
                    if      (clonerange[i] == '_') {
                        Clone_S[clone_rule].to_cut = 0;
                    } else if (clonerange[i] == 'x') {
                        Clone_S[clone_rule].to_cut = 1;
                    }
                }
                if (Clone_S[clone_rule].to_cut == -1) {  /* if it still = -1, then the range doesn't have a valid separator */
                    fprintf(stderr, "\n - invalid range separator in --clone option\n"
                            "please correct the option (%s) and re-run.\n", optarg);
                    exit(1);
                }
                Clone_S[clone_rule].begin = (long)atoi(strsep(&clonerange, "_x"));
                Clone_S[clone_rule].end   = (long)atoi(clonerange);

                /* check whether begin < end; switch them if not */
                if (Clone_S[clone_rule].begin > Clone_S[clone_rule].end) {
                    kk = Clone_S[clone_rule].begin;
                    Clone_S[clone_rule].begin = Clone_S[clone_rule].end;
                    Clone_S[clone_rule].end = kk;
                }
                /* this should end with a check to see if it's at the end of the range string */
                clonerange = strsep(&optarg, ","); 	/* and grab the next range */
                /* clonerange[0] will == \0 when it finishes */
                /* this will also be good for multiple passes (if use multiple --clone options) */

                clone_rule++;  /* incr on 1st time thru */
            }
            /* now check for range conflicts for all preceding clone range sets */

            /* Also, have to check for conflicting rules in different ranges - ie '333_444,366x777'  would be conflicting
            	as the allowed range overlaps with the forbidden range.  Check this by checking for each range, is there
               another range that overlaps..?  if range[0].begin > range[1].begin &&  range[0].begin < range[1].end
               && range[0] is '_' && range [1] is 'x') there's a conflict (also check where range[0].end relative to range[1]
            */
            for (i=1; i<clone_rule; i++) {
                for (j=0; j<i; j++) {
                    if (F.Verbose > 1) {
                        fprintf(stdout, "\nTesting rule %d (%s) vs rule %d (%s) \n", i, Clone_S[i].range_rule, j, Clone_S[j].range_rule);
                    }

                    /* only check sets if they are conflicting - obviously same-type ranges won't conflict */
                    if (Clone_S[i].to_cut != Clone_S[j].to_cut) {
                        if (((Clone_S[i].begin  <  Clone_S[j].begin)  &&
                                (Clone_S[i].end    <  Clone_S[j].begin)) ||
                                (Clone_S[i].begin  >  Clone_S[j].end)) {
                            /*  no need to do anything but print out the commented debugging lines below */
                            /* fprintf(stdout, "\nIn '--clone' range NO CONFLICT between clone ranges %d (%s) and %d (%s)\n\n",
                            i, Clone_S[i].range_rule, j, Clone_S[j].range_rule);
                            fprintf(stdout, "\nClone_S[%d].to_cut [%d] != Clone_S[%d].to_cut [%d]\n", i, Clone_S[i].to_cut, j, Clone_S[j].to_cut);
                            fprintf(stdout,  "&& (((Clone_S[%d].begin [%d] <  Clone_S[%d].begin [%d])\n", i, Clone_S[i].begin, j,  Clone_S[j].begin);
                            fprintf(stdout,  "&& (Clone_S[%d].end [%d]  <  Clone_S[%d].begin [%d]))\n", i, Clone_S[i].end, j,  Clone_S[j].begin);
                            fprintf(stdout,  "|| (Clone_S[%d].begin [%d]  >  Clone_S[%d].end [%d])))\n", i, Clone_S[i].begin, j, Clone_S[j].end);
                            */
                        } else {
                            fprintf(stderr, "\nIn '--clone' range THERE IS A CONFLICT between clone ranges %d (%s) and %d (%s)\n",
                                    i, Clone_S[i].range_rule, j, Clone_S[j].range_rule);
                            if (F.Verbose > 1) {
                                fprintf(stderr, "\nClone_S[%d].to_cut [%d] != Clone_S[%d].to_cut [%d]\n", i, Clone_S[i].to_cut, j, Clone_S[j].to_cut);
                                fprintf(stderr,  "&& (((Clone_S[%d].begin [%ld] <  Clone_S[%d].begin [%ld])\n", i, Clone_S[i].begin, j,  Clone_S[j].begin);
                                fprintf(stderr,  "&& (Clone_S[%d].end [%ld]  <  Clone_S[%d].begin [%ld]))\n", i, Clone_S[i].end, j,  Clone_S[j].begin);
                                fprintf(stderr,  "|| (Clone_S[%d].begin [%ld]  >  Clone_S[%d].end [%ld])))\n", i, Clone_S[i].begin, j, Clone_S[j].end);
                            }
                            badranges = 1;  /* mark it as irretrievably bad  */
                        }
                    }
                }
            }
            if (badranges == 1) {
                fprintf(stdout, "\nPlease correct the bad range(s) and try again.\n");
                exit(1);
            }
            F.Clone = clone_rule - 1;
            break;

        case 2: /*  --logdegens - want to log degeneracies for graphing later */
            F.LogDegens = 1;
            break;

        case 11:  /* what number to start numbering raw sequence at */
            /* don't need any end-checking as the seq can start numbering at any anything, just as long as it
               increases.  Could also make it count DOWN, I guess, but wait for a request to do so */
            if (optarg == 0) BadFlag("-b"); /* die on bad flag */
            j = atoi(optarg);
            if (j != 0) {
                F.NumStart = j - 1;
            } else {
                fprintf (stderr, "'--numstart' flag (%d) can't be 0.\n", j);
                exit(1);
            }

            break;

        case 16: /* graphics plasmid/linear maps - currently no args, just start the plasmid map */
            /* possibly in future,
              0  = postscript
               1 = PNG output (maybe via gs call?
              no args defaults to 0
            */
            if (optarg == 0) { /* then it's alone - no argument to the flag  so we set it to 0  */
                F.GrafMap = 0; /* default is -1; this is analogous to -f (topology)  */

            } else { /* there's something there so we find out what  */
                switch (optarg[0]) {
                case '0':
                case '1':
                    F.GrafMap = atoi(optarg);
                    break;
                default:  {
                    fprintf (stderr, "'--ps' flag requires nothing, '0' (postscript) or '1' (PNG).  Try again \n");
                    exit(1);
                }
                }
            }
            F.Output = 1;
            F.Topo = 0; /* plasmid map implies a circluar plasmid, no? */
            F.LogDegens = 1;
            break;


            /* EXAMPLE flag processing */
        case 12:  /* --example { an integer between 1 and 10 } ex: '--example 4'  */
            /* note that this example requires a numeric option; other kinds of option processing
               are described below, including complex processing to parse out things as complex as
              comma delimited names, etc */
            if (optarg == NULL) BadFlag("--example "); /* die on bad flag */
            j = atoi(optarg);
            F.Output = 1; /* allow this to qualify as 'Output' */
            /* F.Example is initialized at the beginning of tacg.c, at ~line 69 with a call to tacg_SetF(), which
            	sets the default values for all the variables in the F struct (F = Flags).  They can be changed in
               this fn() via processing of the various flags.  */
            if (j > 0 && j < 11) F.Example = j;
            else {
                fprintf (stderr, "'--example' flag must be an integer between 1 and 10.  Try again.\n");
                exit(1);
            }
            /* once you have your flag processed, all you have to do is to test it in tacg.c against the
            	value of F.Example and proceed on the basis of that test. The test is */
            break;

        case 7:
            if (optarg == NULL) BadFlag("--strands "); /* die on bad flag */
            j = atoi(optarg);
            if (j == 1 || j == 2) F.Strands = j;
            else {
                fprintf (stderr, "'--strands' flag must be 1 or 2.  Try again.\n");
                exit(1);
            }
            break;


        case 8:
            F.Tics = 0;
            break;


        case 9:
            F.RawSeq = 1;
            break;


        case 4:     /* --dam -> simulation with dam methylation (GmATC) */
            F.dam = 1;
            break;


        case 5:     /* --dcm -> simulation with dcm methylation (CmC(A/T)GG) */
            F.dcm = 1;
            break;

        case 6:     /* --cost -> in Units/$, so '--cost 300' requests those REs that are >= 300 U/$*/
            if (optarg == NULL) BadFlag("--cost"); /* die on bad flag */
            j = atoi(optarg);
            if (j > 0 ) F.Cost = j;
            else {
                fprintf (stderr, "'--cost' flag must be > 0.  Try again.\n");
                exit(1);
            }
            break;

        case 'W': /* --slidwin #bases  */
            /* does this have to be error-checked for conflicts or silly requests, like
               doing sliding windows for EVERY RE? How will it be used?  How can it be abused? */
            if (optarg == 0) BadFlag("-W (--slidwin)"); /* die on bad flag */
            j = atoi(optarg);
            if (j>2) F.SlidWin = j;
            else  {
                fprintf (stderr, "'-W' flag (--slidwin) value must be > 2.\n");
                Usage(); /* let user see what the values are */
                exit(1); /* and then die */
            }
            break;


        case 'i':  /* --idonly */
            /* '-i' determines whether and how much output will be emitted for those seq's that do not
               have any hits:
               0 - as in older versions - everything gets printed, name and headers regardless if there were
                   any hits or not.
               1 - everything gets printed only if there are hits.  If no hits, only the name gets printed.
               2 - ONLY the name gets printed and ONLY if there are hits.
            */
            j = atoi(optarg);
            if (j>=0 && j<3) {
                F.Output = 1; /* allow this to qualify as 'Output' */
                F.IDonly = j;
            } else {
                fprintf (stderr, "\n'-i' flag requires an integer: 0 - 2.  try again..\n");
                exit(1);
            }
            break;


        case 1:  /* --silent */
            F.Silent = 1;
            F.Degen = 3;   /* and set the Degeneracy-handling flag to process the expected degens*/
            break;


            /* 'g' - include the gel map; currently a text/graphic, width-settable with -w flag  */
            /* flag has 2 arguments: '-g minlog,maxlog', where minlog = some exponenent of 10 and indicates
               where the gel should start - no longer just 10 or 100 - can be ANY exponent of 10 and
               maxlog indicates where the gel should range up to */
        case 'g':
            if (optarg == 0) BadFlag("-g"); /* die on bad flag */
            if (F.Width > MAX_BASES_PER_LINE) Errw1(); /* check for -w1 condition */
            F.Output = 1; /* set F.Output to indicate that some outout has been requested */
            tok = strsep(&optarg, ",:.");
            j = atoi(tok);
            /* z = (double)j; */
            if (j >= 10) {
                F.GelLO = (long)log10((double)j); /* truncates to the integer part */
            } else {
                fprintf (stderr, "The '-g' flag form is: '-g m,M' (where m >= 10, M > m).\n");
                exit (1);
            }
            tok = strsep(&optarg, ",:.");
            if (tok == 0) {
                F.GelHI = 0;  /* this means that F.GelHI will be set to the sequence length when determined  */
            } else {
                k = atoi(tok);
                if (k > j) {
                    F.GelHI = (long)log10((double)k); /* truncates to the integer part */
                } else {
                    fprintf (stderr, "The '-g' flag form is: '-g m,M' (where m >= 10, M > m).\n");
                    exit (1);
                }
            }

            break;


        case 'r':  /* aka --regex='RxName:RxPattern'  or -r 'FILE:FileName'  (NEEDS the single quotes)
                       so obviously you can't name a regex pattern 'FILE' */
            /* this stanza also has to duplicate a lot of the functionality of the '-p' stanza,
            	as the regex's read in here have to be parsed and then written out to a temp file
            	as well as a log file.  (The temp file is then read in as a standard regex file and
            	all the error handling is done by the ReadRegex() - Because of the overlapping
            	variables to the -p flag, they are not compatible, so have to warn each other.  */
            /* there must be a regex file waiting to be read */
            if (PCRE == 1) {
                RebFile = Flags_case_r(optarg, tmpfp, datestamp);    /* see below */
            } else {
                fprintf(stderr, "tacg compiled without regex support - install pcre and re-compile\n");
                exit(1);
            }
            break;


            /* 'R' uses SearchPaths() to search for an alternative rebase file - SearchPaths() checks if
                the name starts with a '/','~', or '.' and uses the full path and name if entered
                without leading path specifiers, else checks thru a series of paths depending on
                environment variables and returns the pathname of whatever it finds.  If SearchPaths()
                finds nothing, it returns NULL and everything dies */
        case 'R':
            if (optarg == NULL) BadFlag("-R"); /* die on bad flag */
            if (F.AltPatFile == 2 || F.Regex != 0) { /* then -p or -x 'option' has already been requested. They're incompatible */
                fprintf(stderr,"!! '-p (or -x with options) and -R flags are incompatible.\n");
                exit(1);
            } else if (F.AltPatFile == -2) { /* if nothing has been set so far (-2 is orig value)... */
                F.AltPatFile = 1; /* mark it so I can reconstruct the command line */
                if ((RebFile = SearchPaths(optarg, "DATABASE")) == NULL) {
                    fprintf(stderr,"Can't find the alternative DATABASE file!! - Check spelling, ownership, existance, etc - Bye!!\n");
                    exit(1);
                }
                F.Pat=2; /* set -p flag so that it won't be used */
            }
            break;


        case 'X':  /* --extract b,e,revcomp */
            Flags_case_X(optarg);
            break;


        case '#':  /* Matrix Matching; arg is the Cutoff in (int) percentage */
            if (optarg == 0) BadFlag("-#"); /* die on bad flag */
            j = atoi(optarg);
            if (j>0 && j<=100) F.Matrix = j;
            else {
                fprintf (stderr, "'-#' (Matrix CutOff) value must be between 0 and 100.  Try again\n");
                exit(1);
            }
            break;


        case 'f':  /* form or topology  - circular or linear */
            if (optarg == 0) BadFlag("-f"); /* die on bad flag */
            switch (optarg[0]) {
            case '0':
            case '1':
                F.Topo = atoi(optarg);
                break;
            default:  {
                fprintf (stderr, "'-f' flag requires '0' (circular) or '1' (linear).  Assuming linear \n");
                F.Topo = 1;
                break;
            }
            }
            break;


        case 'p':  /* include patterns from the command line, sharing some of the vars from the -x flag above */
            Flags_case_p(optarg, tmpfp, datestamp); /* see below */
            break;


            /* process '--rule' names into a char *array and return that to main() as well - see below */
        case 15:
            Flags_case_15(optarg, rulestring, SelEnz, &SECnt, datestamp);
            F.Output = 1;
            break;

            /* -x used to pass names *explicitly* to pull patterns out of a file for pattern-finding */
        case 'x':
            Flags_case_x(optarg, SelEnz, &SECnt);
            break;


            /*            vv-dep--|------|      */
            /* -P N1,[+-][lg]#####[-#####],N2    where ##### are non negative #s  */
            /*      ^|-----  NMess  -----|^     */
        case 'P': /* here we go - the final frontier - pattern proximity matching */
            Flags_case_P(optarg, SelEnz, &SECnt);
            break;


        case 'c':  /* sort REs in output by # cuts; by order of listing in REBASE if not */
            F.StriderSort = 1;
            break;


        case 'F': /* How to present fragment data - include or exclude 0-3;
         				it's -1(d) (no frag data) if it wasn't called*/
            if (optarg == 0) BadFlag("-F"); /* die on bad flag */
            F.Output = 1; /* set F.Output to indicate that some outout has been requested */
            j = atoi(optarg);
            if (j>=0 && j<4) {
                F.Frags = j;
                /*                fprintf(stderr, "Frags set to %d\n", F.Frags); */
            } else fprintf (stderr, "The '-F' flag value must be: \n"
                                "           0 (don't print fragment data), \n"
                                "           1 (print fragments sorted by order)\n"
                                "           2 (print fragments sorted by size), or \n"
                                "           3 (print fragments sorted both ways), \n"
                                "not \'%s\'.\n", optarg);
            break;


        case 'D':  /* consider input seq degenerate AND do extended mapping even if key hexamer is degenerate */
            if (optarg == 0) BadFlag("-D"); /* die on bad flag */
            j = atoi(optarg);
            if (j>=0 && j<5) F.Degen = j;
            else fprintf (stderr, "The '-D' flag value must be: \n"
                              "  0 - FORCE exclusion of IUPAC and other characters in sequence\n"
                              "  1 - default - cut as nondegenerate unless IUPAC characters found;\n"
                              "        then cut as '-D4'\n"
                              "  2 - allow IUPAC characters; ignore them in the KEY hexamer,\n"
                              "        but match outside of KEY)\n"
                              "  3 - allow IUPAC characters; find only exact matches\n"
                              "  4 - allow IUPAC characters; find ALL POSSIBLE matches\n"
                              "not '%s' as you entered.\n", optarg);
            break;


        case 's': /* print summary stats of enz's that don't map and how many times they do if any */
            F.Output = 1; /* set F.Output to indicate that some outout has been requested */
            F.Summary = 1;
            break;


        case 'H':  /* include HTML tags if present, exclude if not */
            /* 0 = FULL HTML page with TOC
               1 = HTML content only with TOC, no headers
              no args defaults to 0
            */
            if (optarg == 0) { /* then it's alone - no argument to the flag and we want -H 0, so we set it to 0  */
                F.HTML = 0;
            } else { /* there's something there so we find out what  */
                switch (optarg[0]) {
                case '0':
                case '1':
                    F.HTML = atoi(optarg);
                    break;
                default:  {
                    fprintf (stderr, "'-H' flag requires '0' (Full HTML) or '1' (partial HTML w/ TOC).  Try again \n");
                    exit(1);
                }
                }
            }
            break;


        case 'O': /* sets which ORFs to analyze and minimum ORF length */
            Flags_case_O(optarg);
            break;


        case 'n':  /* magnitude of RE recognition site - 4 cutter, 5 cutter, etc */
            if (optarg == NULL) {
                BadFlag("-n"); /* die on bad flag */
            }
            j= atoi(optarg);
            if (j>2 && j<11) F.Mag = j;
            else {
                fprintf (stderr, "'-n' flag requires an integer between 3 and 10.  Assuming '-n6'..\n");
                F.Mag = 6;
            }
            break;


        case 'o':  /* Type of overhang to consider - omit flag to consider all */
            if (optarg == 0) BadFlag("-o"); /* die on bad flag */
            j = atoi(strsep(&optarg, ",:."));
            if (j==5 || j==3 || j==1 || j==0) F.Overhang = j;
            else  {
                fprintf (stderr, "\n'-o' flag must be 5, 3, 1, or 0; omit to include all overhangs...Try again\n");
                exit(1);
            }
            /* to filter on overhang length as well */
            if (j != 0) { /* only if there's an overhang */
                tstr = strsep(&optarg, ",:."); /* do it again to get (optional) overhang length */
                if (tstr != NULL) {
                    j = atoi(tstr); /* do it again to get (optional) overhang length */
                    if (j > 0 && j < 7) {
                        F.Overhanglen = j;
                    } else {
                        fprintf (stderr, "\nTo filter on overhang length as well as overhang, the '-o' flag must be \n"
                                 "entered as -o5,#, where '#' is the length of the overhang wanted between \n"
                                 "1-6 inclusive. Try again.\n");
                        exit(1);
                    }
                }
            }
            break;  /* above assumption doesn't require changing anything */


        case 'G':  /* Graphics data - argument is the number of bases per bin,format */
            if (optarg == 0) BadFlag("-G"); /* die on bad flag */
            F.Output = 1; /* set F.Output to indicate that some outout has been requested */
            temp = DownCase(strsep(&optarg, ",:.")); /* strsep replaces strtok() */
            len = atoi(temp); /* len = Nbins */
            if (len > 0) F.GrfBBin = len;
            else {
                fprintf (stderr, "'-G' flag option (part1) must be 0 or greater.\n");
                exit(1);
            }
            temp = (strsep(&optarg, ",:.")); /* next call gets second argument */
            len = strlen(temp);
            if (len == 1) {
                temp = DownCase(temp);
                ch = temp[0];
                if (ch == 'x') F.GrfFmt = 1;
                else if (ch == 'y') F.GrfFmt = 2;
                else if (ch == 'l') F.GrfFmt = 3;
                else {
                    fprintf(stderr, "Bad option for '-G'\n");
                    exit(1);
                }
            } else {
                fprintf(stderr, "-G flag option (part 2) must be in form '-G Bin_size,X|Y|L'\n");
                exit(1);
            }
            break;


        case 'm':  /* Minimum # of cuts to be considered */
            if (optarg == 0) BadFlag("-m"); /* die on bad flag */
            j = atoi(optarg);
            if (j>=0) F.Min = j;
            else fprintf (stderr, "'-m' flag must be 0 or greater.\n");
            if (F.Max != 0 && F.Max -j < 0) {
                fprintf (stderr, "'-m' flag (%d) must be less than or equal to '-M' flag (%d)."
                         " Assuming you want all cuts.\n", j, F.Max);  /* Doesn't require changing anything */
            }
            break;


        case 'M':  /* Maximum # of cuts to be considered */
            if (optarg == 0) BadFlag("-M"); /* die on bad flag */
            j = atoi(optarg);
            if (j>=0) F.Max = j;
            else fprintf (stderr, "'-M' flag must be 0 or greater.\n");
            if (F.Min != 0 && j - F.Min < 0) {
                fprintf (stderr, "'M' flag (%d) must be greater than or equal to 'm' flag (%d). "
                         "Assuming you want all cuts.\n", j, F.Min);
            }
            break;


        case 'b':  /* Beginning of subsequence */
            /* !!! check that a -b that requests a begin AFTER the total length of the sequence doesn't crash it */
            if (optarg == 0) BadFlag("-b"); /* die on bad flag */
            j = atoi(optarg);
            if (j>0) F.SeqBeg = j;
            else  {
                fprintf (stderr, "'-b' flag (%d) must be >0.\n", j);
                exit(1);
            }
            /* check if b < e and that sequence is long enough (>4) */
            if ((F.SeqEnd != 0) && (F.SeqEnd - j < 3)) {
                fprintf (stderr, "'-b' flag (%d) must be less than value of '-e' flag (%ld)\n "
                         "and subsequence must be at least 4 bps", j, F.SeqEnd);
                exit(1);
            }
            break;


        case 'e':  /* End of subsequence */
            if (optarg == 0) BadFlag("-e"); /* die on bad flag */
            j= atoi(optarg);
            if (j>=0) F.SeqEnd = j;
            else  {
                fprintf (stderr, "'-e' flag must be >= 0.\n");
                exit(1);
            }
            /* check if e > b and that sequence is long enough (>4) */
            if (F.SeqBeg != 1 && F.SeqEnd != 0 && j - F.SeqBeg < 4) {
                fprintf (stderr, "'-b' flag (%ld) must be less than value of '-e' flag (%ld)\n "
                         "and subsequence must be at least 4 bps", F.SeqBeg, F.SeqEnd);
                exit(1);
            }
            break;


        case 'T': /* Translate in requested frames in 1 (W) or 3 letter (Trp) code  */
            Flags_case_T(optarg);
            break;


        case 'h':
            Usage();
            exit(1);
            break;


        case 'v':
            fprintf (stdout, "\ntacg Version %s\nWritten by Harry Mangalam (hjm@tacgi.com)\n"
                     "\nCopyright  1996-2011 Harry Mangalam, tacg Informatics"
                     "\nThis software is licensed under the Free Software Fndn GPL 3.0"
                     "\nSee the source for copying conditions.  There is NO warranty;"
                     "\nnot even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n",
                     TACG_VERSION);
            fprintf (stdout, " Executable built: %s \n on kernel: %s\n GCC version: %s\n", TACG_BUILD_DATE, TACG_BUILD_PLATFORM, TACG_GCC_VER);
            exit(1);
            break;


        case 'V': /* Verbose - print error messages to stdout; otherwise shaddup */
            if (optarg == 0) BadFlag("-V"); /* die on bad flag */
            j = atoi(optarg);
            if (j>0 && j<4) F.Verbose = j;
            else  {
                fprintf (stderr, "'-V' flag value must be 1-3, not %s.\n", optarg);
                Usage(); /* let user see what the values are */
                exit(1); /* and then die */
            }
            break;


        case 'C': /* Codon usage table to use - 0(d) to 7 */
            if (optarg == 0) BadFlag("-C"); /* die on bad flag */
            j = atoi(optarg);
            if (j>=0 && j<MAX_ORGS) F.CodTab = j;
            else  {
                fprintf (stderr, "'-C' flag value must be > 0 and < %d, not %s.\n", MAX_ORGS, optarg );
                exit(1); /* and then die */
            }
            break;


        case 'L': /* whether to print linear map; -L prints the linear map; omitting it doesn't print it */
            if (F.Width > MAX_BASES_PER_LINE) Errw1();  /* check for -w1 condition */
            F.Output = 1; /* set F.Output to indicate that some outout has been requested */
            F.LinMap = 1;
            break;


        case 'S':  /* whether to print cutsites; if specified, print cutsites; otherwise don't */
            F.Output = 1; /* set F.Output to indicate that some outout has been requested */
            if (optarg == 0) {
                F.Sites = 1; /* std Sites output */
            } else {
                j = atoi(optarg);
                if (j < 1 || j > 2) BadFlag("-S");
                else if (j == 1) F.Sites = 1; /* std Sites output */
                else if (j == 2) F.Sites = 2; /* extended Sites output ('-' for nonpals in bottom strand */
            }
            break;


        case 'w': /* change the output width */
            if (optarg == 0) BadFlag("-w"); /* die on bad flag */
            j = atoi(optarg);
            if ((j >= 60 && j <= MAX_BASES_PER_LINE) || j == 1)  {
                if (j == 1) {
                    /*        -L               -g               -l               -P        */
                    if (F.LinMap == 1 || F.GelLO != 0 || F.Ladder == 1 || F.Prox != 0 ) {
                        Errw1(); /* check if they've been set already */
                    }
                    F.Width = 10000000; /* for single lines */
                } else F.Width = (long) j / 15*15; /* chew it down to mod 15 */
            } else {
                fprintf (stderr, "'-w' value must be 1 or between 60 and MAX_BASES_PER_LINE (usually 210);\n "
                         "program truncates to a number exactly divisible by 15 (178 truncates to 165).\n");
                exit (1); /* and die cleanly */
            }
            /* and so it's passed back to main(); effectively replaces the BASES_PER_LINE #define */
            break;


        }   /* end of switch() statement */
    }  /*    while (--argc > 0 && (*++argv)[0] ... */

    /* DELETE TO ### following stanza when verify */
    /* only if it's been opened can we close it, dummy */
    /*
       if (pCnt != 0){
          fclose(fppat);
       }
     */
    /* close the pattern log; temp file pointer stays open until end of ReadEnz() */
    /* ### */

    /* Do some final error checking and advice to the confused */
    if (F.LinMap == 0 && F.XLinFr != 0) {
        fprintf (stderr, "\nINPUT ERROR: '-T' (Translate) cannot produce output without '-L' (Linear Map). \n"
                 "Do you want the '-O' flag? Type 'tacg -h' for usage; 'man tacg' for man page.\n");
    }

    if (RebFile[0] == '\0') {  /* then it hasn't been set already, so it needs the default */
        /* if F.AltPatFile == -2, it hasn't been touched, so we want rebase.data */
        if (F.AltPatFile == -2) {  /* this should be redundant - if something changed [17], it should have named a file */
            RebFile = SearchPaths("rebase.data", "REBASE"); /* default for REBASE filename */
        }
    }
    /* below means that -R hasn't been called yet, so we need to set the default 'matrix.data' */
    /* if user HAS set '-#' and HASN'T called for an alt Matrix file ('-R'), then look for the default 'matrix.data'
    	and if we can't find it, then issue a warning */
    if (F.Matrix != 0 && F.AltPatFile == -2) {
        if ((RebFile = SearchPaths("matrix.data", "MATRIX")) == NULL) {
            fprintf(stderr,"Can't find 'matrix.data' - you need to specify an alternative with the '-R' flag!\n");
            exit(1);
        }
    }

    /*   Clone[MAX_NUM_CLONE_RANGES].to_cut = clone_rule;  return the number of ranges to consider in the last el */
    /*   F.Clone = clone_rule;  store the number of ranges to consider in the Flag struct */
    return RebFile;
}    /* end of  SetFlags */

/*  ============================ Additional subs to modularize SetFlags =====================================  */

/* ----------------------------------------------------------------------------------------------------------- */
/* sets which ORFs to analyze and minimum ORF length */
void Flags_case_O(char *optarg)
{
    int i, j, l, OK, m;
    char s[80];

    if (optarg == 0) BadFlag("-O"); /* die on bad flag */
    F.Output = 1; /* set F.Output to indicate that some outout has been requested */
    memset(s,' ',80);
    if (strlen(optarg) > 2) { /* check if there's something useful in optarg */
        /* break optarg into separate args and plunk into frames[] */
        i=0;
        j=0; /* j starts at 0 to make the frames to be ORFed fit in F correctly */

        while (optarg[i] != ',' && optarg[i] != '\0') {
            /* have to make the individual char a string 1st */
            if (optarg[i] == 'x' || optarg[i] == 'X') {
                F.OrfXtra = 1;
            } else {
                s[0] = optarg[i];
                s[1] = '\0';
                m = atoi(s);
                l=0;
                OK = 1;
                while (m != F.ORFs2Do[l] && l <= j) {
                    l++;
                }
                for (l=0; l<j; l++) {
                    if (m == F.ORFs2Do[l]) {
                        OK = 0;    /* if someone entered -O11123,##, instead of -O123,# */
                    }
                }
                if (m > 0 && m < 7 && j < 6 && OK) { /* as long as m values are OK and we haven't gone up too far or dup'ed #s*/
                    F.ORFs2Do[j++] = m; /* F.ORFs2Do[] passes frame-2b-ORFed indicators (NOT nec in order) */
                } else {
                    fprintf(stderr, "'-O': frames must be between 1 and 6 inclusive & entered only 1X; try again\n");
                    exit(1);
                }
            }
            i++;
        }
        if (optarg[i++] == ',') {
            strcpy(s,optarg+i);
            m = atoi(s); /* this should now be a properly terminated string */
            if (m > 0) F.ORFs = m; /* store the value in the '-O' position in flags */
        } else { /* we got here cuz we hit '\0' before getting the min ORF value */
            fprintf(stderr, "'-O' flag is missing the Min ORF length value - try again\n");
            exit(1);
        }
    } else {
        fprintf(stderr, "'-O' flag usage is '-O frameframeframe,MinOrfLength';\n"
                "ie. '-O 135,55' (data from frames 1,3,5 with a minimum ORF length of 55 aas - try again\n");
    }
}

/* ----------------------------------------------------------------------------------------------------------- */
void Flags_case_X(char *optarg)     /* --extract b,e,revcomp */
{
    int e;
    if (optarg == 0) BadFlag("-X (--extract)"); /* die on bad flag */
    if (strlen(optarg) < 5) { /* check if there's something in optarg */
        fprintf(stderr, "!! '--extract' flag option is too short or missing.  The '--extract' flag usage is:\n"
                "'-X (or --extract) b,e,revcomp.\nwhere b=beginning, e=end, revcomp=0|1  indicating:\n"
                "b - the seq taken begins 'b' bases before the approx pattern center\n"
                "e - the seq taken ends 'e' bases after the approx pattern center\n"
                "revcomp - whether you want (1) or don't want (0) the reverse\n"
                "complement taken if the pattern is in the reverse orientation\n");
        exit(1);
    }
    F.XtBeg = (long)atoi(strsep(&optarg, ",:.")); 	/* strsep replaces bizarro strtok */
    if (F.XtBeg < 0) { // was 1 - hjm519
        fprintf(stderr, "in '--extract b,e,revcomp', 'b' must be 0 or greater.\n");
        exit(1);
    }
    F.XtEnd = (long)atoi(strsep(&optarg, ",:.")); 	/* 2nd strsep is called same as 1st */
    if (F.XtEnd < 0) { //was 1 - hjm519
        fprintf(stderr, "in '--extract b,e,revcomp', 'e' must be >0.\n");
        exit(1);
    }
    e = atoi(strsep(&optarg, ",:."));
    if (e == 1) { /* 3rd call continues the tokenizing */
        F.XtBeg = -1; /* to mark it using only 2 vars */
        //fprintf(stdout, "in '--extract b,e,revcomp', F.XtBeg = [%d]\n", F.XtBeg);
    } else if (e != 0) {
        fprintf(stderr, "in '--extract b,e,revcomp', 'revcomp' must = 0 or 1.\n");
        exit(1);
    }
    F.Xtract = 1;
    F.Output = 1;  /* indicate that some output has been requested */
}

/* ----------------------------------------------------------------------------------------------------------- */
/* Translate in requested frames in 1 (W) or 3 letter (Trp) code  */
/* following stanza replaces the -t/-T silliness of previous versions, but keeps the
   choices of 1, 3, or 6 frames vs single frame choices as in the '-O1345,55' case above */

void Flags_case_T(char *optarg)
{

    char	s[80];
    int	i,m;

    if (optarg == 0) BadFlag("-T"); /* die on bad flag */
    memset(s,' ',80);
    if (strlen(optarg) > 2) { /* check if there's something useful in optarg */
        /* break optarg into separate args and plunk into frames[] */
        i=0;
        /* should only be 1 pass thru the following 'while' - use it like this in case want to
           expand or change it later to allow the 123456 approach */
        while (optarg[i] != ',' && optarg[i] != '\0') {
            /* have to make the individual char a string 1st */
            s[0] = optarg[i];
            s[1] = '\0';
            m = atoi(s);
            if (m==0 || m==1 || m==3 || m==6) { /* as long as m values are OK and we haven't gone up too far */
                F.XLinFr = m; /* F.XLinFr holds the frames indicator */
                F.LinMap = 1; /* and set the Linear Map flag to 1 as well just in case */
            } else {
                fprintf(stderr, "'-T': frames have to be '1', '3', or '6'; try again\n");
                exit(1);
            }
            i++;
        }
        if (optarg[i++] == ',') {
            strcpy(s,optarg+i);
            m = atoi(s); /* this should now a properly terminated string */
            if (m == 1 || m == 3) F.Lab1or3 = m; /* store the value in the '-T' position in flags */
        } else { /* we got here becuz we hit '\0' before getting the min ORF value */
            fprintf(stderr, "in '-T' flag, the Label code value is not 1 or 3 - try again\n");
            exit(1);
        }
    } else {
        fprintf(stderr, "'-T' flag usage is: -T {[0|1|3|6],[1|3]}\n"
                "(Translates frames 1, 1-3, or 1-6) w/ Linear Map w/ 1 or 3 letter codes.\n");
    }
}


/* ----------------------------------------------------------------------------------------------------------- */
/* -x used to pass names *explicitly* to pull patterns out of a file for pattern-finding */
void Flags_case_x(char *optarg, struct SE_struct SelEnz[MAX_SEL_ENZ], int *SECnt)
{
    char *temp = "             "; /* temp holder for RE name in the selection string */
    int max, i, LastNamei=0;

    if (optarg == 0) BadFlag("-x"); /* die on bad flag */
    if (strlen(optarg) > 1 && F.Prox == 0) { /* check if there's something in optarg and !using   */
        /* break optarg into separate chars and plunk into *SelEnz[10] */
        temp = DownCase(strsep(&optarg, ",:.")); /* instead of strtok() */
        /* SECnt=0  - Selected Enz Count; NOT reset at each pass thru the case 'x' stanza
        	it keeps track of ALL the REs requested from all the '-x' requests */
        max = MAX_SEL_ENZ;     /* shouldn't have to do this, but it failed otherwise */
        F.SelEnzUsed = 1;  /* and be explicit about this  */

        while (temp != NULL && *SECnt < max) {
            temp = DownCase(temp);  /* lclint says this is not kosher, but tracing it looks ok */

            if (strcmp(temp,"c") == 0) {  /* test for the 'Combined' tag */
                F.AllInc = 1; /* set the 'Combined' or 'All' flag (noted as AI in the flags description */
                /* fprintf(stderr, "\nF.AllInc set to %d\n", F.AllInc); */
                break; /* since ',c' has to be the last entry, break here */
            }

            i=0;    /* check if that RE name has already been entered into SelEnz[] */
            while (i < *SECnt && strcmp(SelEnz[i].PName, temp) != 0) i++;
            if (i == *SECnt) { /* and if not, maybe copy it in (if it's not a '=')  */
                /* at this pt, have to check whether the name is a '='; if it is, then the PREVIOUS RE
                   (whose index is LastNamei) is the one that we want to be keyed on in Hookey(); as
                   coded below, *SECnt is used as the index - only useful in referring to SelEnz[]  but no
                   real error checking being done here yet - check for 2 '='s, etc.  SelEnz is passed back
                   and forth and flags is global, so can extract it in main() */
                if (strcmp(temp,"=") == 0) {
                    F.Hookey = LastNamei; /* set the 'Hookey' flag  */
                    F.AllInc = 1; /* if want the Hookey() need to generate the combined data as well */
                    F.Output = 1; /* set F.Output to indicate that some outout has been requested */
                } else { /* it's just a regular name, so process it as normal */
                    LastNamei  = *SECnt;  /*  marking for Hookey just in case */
                    strcpy(SelEnz[(*SECnt)++].PName, temp); /* and copy it into SelEnz to be extracted in main() */
                }
            }
            /* new bit to allow specification of COMBINED cuts, as sometimes you don't want to
            	as when you're looking for several TF sites (which don't cut) rather than REs
            	This means of course, that if your pattern is called 'c', it will have unexpected
            	results */
            temp = (strsep(&optarg, ","));  /* and get the next token/RE Name */
        }

        if (*SECnt >= max) fprintf(stderr, "'-x' flag has too many names; using only the 1st %d\n", max);
        /*     F.Xplicit++;  incr to count the # of iterations of -x */
        F.Xplicit = 1; /*  incr to count the # of iterations of -x */
    } else { /* whine about a badly composed enz name string */
        fprintf(stderr, "The '-x' flag usage is '-x enz,enz,enz...'; try again\n");
        exit(1);
    }
}

/* ----------------------------------------------------------------------------------------------------------- */
/*            vv-dep--|------|      */
/* -P N1,[+-][lg]#####[-#####],N2    where ##### are non negative #s  */
/*      ^|-----  NMess  -----|^     */

void Flags_case_P(char *optarg, struct SE_struct SelEnz[MAX_SEL_ENZ], int *SECnt)
{

    int PPCnt=-1, optlen, i, j, e, f, len;
    char *temp = "             "; /* temp holder for RE name in the selection string */
    char *NMess = "                              ";

    if (optarg == 0) BadFlag("-P"); /* die on bad flag */
    if (F.Width > MAX_BASES_PER_LINE) Errw1(); /* check for -w1 condition */
    F.Output = 1; /* set F.Output to indicate that some outout has been requested */
    F.SelEnzUsed = 1;  /* and be explicit about this  */

    optlen = strlen(optarg);
    if (optlen > 5 && PPCnt < MAX_PP && F.Xplicit <= 0 && F.Rules == 0) {
        /* check if there's something in optarg
          and haven't done too many and don't conflict with '-x', '--rule' */
        e=f=0;
        PPCnt++;
        F.Xplicit = -1;
        F.Prox++;/* reset e,f; incr the P, F.Prox counter */
        /* break optarg into separate elements - the 2 names and then the NMess */
        /* going to have to use the '-x' procedure, if not all the data structures to match the REBASE names */
        /* Some initializations for those vars that need it */
        // and calloc some mem for PP.
        if ((PP = calloc(INIT_PPs,sizeof(*PP))) == NULL) {
            BadMem("PP fails malloc", 1);
        }
        PP[PPCnt].upstr=0;
        PP[PPCnt].gt =-1; 		/* 'less than' is the default */
        PP[PPCnt].range =-1;  	/* NO range is the default */

        /* grab the command line string for attaching to the output */
        PP[PPCnt].optline = (char *) calloc(optlen+1, sizeof(char));
        if (PP[PPCnt].optline==NULL) BadMem("calloc fails at optlen", 1);
        strcpy(PP[PPCnt].optline, optarg);

        /* co-opt the '-x' SelEnz[] for this one as well */
        temp = DownCase(strsep(&optarg, ",:.")); /* strsep replaces strtok */
        len = strlen(temp);
        if (len > 10) {
            fprintf(stderr, "In '-P' option, name %s is too long, truncating to 10 chars.\n", temp);
            temp[10] = '\0';
            len = 10;
        }
        PP[PPCnt].N[0] = (char *)calloc(len+1, sizeof(char));
        if (PP[PPCnt].N[0] == NULL) BadMem("!calloc fails\n", 1);
        strcpy(PP[PPCnt].N[0],temp);
        strcpy(SelEnz[(*SECnt)++].PName, PP[PPCnt].N[0]);

        NMess = strsep(&optarg, ",:."); /* gets the Number Mess */
        if (NMess==NULL) fprintf(stderr,"Bad format in '-P' flag.\n");

        temp = DownCase(strsep(&optarg, ",:.")); /* gets the 2nd name */
        len = strlen(temp);
        if ( len > 10) {
            fprintf(stderr, "In '-P' option, name %s is too long, truncating to 10 chars.\n", temp);
            temp[10] = '\0';
            len = 10;
        }
        PP[PPCnt].N[1] = (char *)calloc(len+1, sizeof(char));
        if (PP[PPCnt].N[1] == NULL) BadMem("!calloc fails\n", 1);
        strcpy(PP[PPCnt].N[1],temp);
        strcpy(SelEnz[(*SECnt)++].PName, PP[PPCnt].N[1]);

        /* so now have optarg broken into usable tokens */
        len = strlen(NMess); /* how much pain is there? */
        j = 0;
        for (i=0; i<len; i++) {
            switch (NMess[i]) {
            case '+':
                if (i== 0) PP[PPCnt].upstr = -1; /* can only be used once at start */
                else fprintf(stderr,"'+' out of place in '-P' option\n");
                break;

            case '-':
                if (i==0) PP[PPCnt].upstr = 1; /* can only be used in this context once at start */
                else { /* indicates the end of Dlo, so process it as a decent # */
                    temp[j] = '\0';
                    PP[PPCnt].Dlo = abs(atoi(temp));
                    f = -1;
                    memset(temp, 0, 13);
                    j = 0; /* re/set the related vars */
                    if (PP[PPCnt].gt != 0) {
                        /* fprintf(stderr,"in '-P' option 'l or g' illegal when specifying a range\n"
                                   	"unsetting it and continuing...\n"); */
                        PP[PPCnt].gt = 0;
                    }
                }
                break;

            case 'l':
                PP[PPCnt].gt = -1;
                if (i>1) fprintf(stderr,"'l' out of place in '-P' option, but continuing...\n");
                break;

            case 'g':
                PP[PPCnt].gt = 1;
                if (i>1) fprintf(stderr,"'g' out of place in '-P' option, but continuing...\n");
                break;

            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                temp[j++] = NMess[i]; /* must be part of one of the distances, so add it */
                break;

            default:
                fprintf(stderr,"Bad char (%c) in '-P' option.\n", NMess[i]);
                break;
            }
        }
        temp[j] = '\0';  /* regardless */
        if (f == -1) {
            PP[PPCnt].Dhi = abs(atoi(temp)); /* if f= -1, this is the *2nd* number.. */
            if (PP[PPCnt].Dhi < PP[PPCnt].Dlo) {  /* swap em - maybe user made a mistake  */
                e = PP[PPCnt].Dhi;
                PP[PPCnt].Dhi = PP[PPCnt].Dlo;
                PP[PPCnt].Dlo = e;
                if (F.Verbose > 0) fprintf(stderr, "in '-P' option, Lo > Hi.  Mistake?  Swapping and continuing.\n");
            }
            PP[PPCnt].range = abs(PP[PPCnt].Dhi - PP[PPCnt].Dlo); /* the range of the distances */
        } else {
            PP[PPCnt].Dlo = atoi(temp);   /* else it's the 1st */
            PP[PPCnt].range = -1; /* -1 indicates NO range */
        }
        /*set the searchtype here, rather than decide on the fly later */
        if (PP[PPCnt].upstr==0) {
            if (PP[PPCnt].gt==-1) {
                if (PP[PPCnt].range==-1) PP[PPCnt].st = '0';   /* upstr:0 gt:-1 range:-1 - default*/
            } else if (PP[PPCnt].gt==1) PP[PPCnt].st = '1';   /* upstr:0 gt: 1 range:-1 */
            else  						 PP[PPCnt].st = '2';   /* upstr:0 gt: 0 range:calc */
        } else if (PP[PPCnt].upstr==1) {
            if (PP[PPCnt].gt==-1) {
                if (PP[PPCnt].range==-1) PP[PPCnt].st = '3';   /* upstr:1 gt:-1 range:-1 */
            } else if (PP[PPCnt].gt==1) PP[PPCnt].st = '4';   /* upstr:1 gt: 1 range:-1 */
            else  						 PP[PPCnt].st = '5';   /* upstr:1 gt: 0 range:calc */
        } else if (PP[PPCnt].gt==-1) {
            if (PP[PPCnt].range==-1) PP[PPCnt].st = '6';   /* upstr:-1 gt:-1 range:-1 */
        } else if (PP[PPCnt].gt==1) PP[PPCnt].st = '7';   /* upstr:-1 gt: 1 range:-1 */
        else 					 		 PP[PPCnt].st = '8';   /* upstr:-1 gt: 0 range:calc */

        if (strcmp(PP[PPCnt].N[0], PP[PPCnt].N[1]) == 0) PP[PPCnt].upstr = 0; /* if same, there is no 'upstream' */
    } /* if (strlen(optarg) > 1 && PPCnt < MAX_PP)  ...  */
}

/* --------------------   Self mumblings for case 15   ----------------------------
   --rule 'RuleName,((LabA:m:M&LabB:m:M)|(LabC:m:M&LabD:m:M))|((LabE:m:M)&LabF:m:M),window'

   The SelEnz struct also used by -P to filter enz's from REBASE; ergo -x, --rule (aka case
   15)  and -P are incompatible.


   using strsep() to tokenize the parts. If there's not 3 parts to the 'Name:min:Max'
   triplet, the one numeric is the min flag (name:3).   If there's 3 parts but the middle
   is NULL, the last is Max (name::45).  There has to be a name, obviously.

   6.1.99 - changing this again as it's getting WAY too complex for a single flag.
   Bringing back the old '-x, which can handle names, the Hookey tag, and the 'Combined'=C
   tag, and breaking the new options (T/F logic, per-pattern min:Max, and paren grouping
   into its own case (case 15) = '--rule' tokenize 1st on '|&', strip parens off the
   front, then off the back, then tokenize on ':' and stuff the  Name, min, Max into the
   SelEnz struct.  However, the WHOLE optarg for  this flag has to be passed back as well
   to be parsed into TRUTH values in SlidWin() and  perhaps elsewhere.
   8.3.99 - WHY does all this need to return 'rulestring'?  We're already parsing the
   commandline down to names, so couldn't these names be included in SelEnz??  The only
   reason for doing this  before was that the 2 options were separate.  As it stands now,
   they will be exclusively separate, so they can both use SelEnz, no?  Yes.
   8.9.99 --rule has to invoke -W / --slidwin, so that all the rule-evaluating logic gets
   included. If it's not set explicitly to something < the sequence length, implicitly set
   to the entire sequence. This can be figured out at run-time by values of F.Xplicit and
   F.SlidWin.   If F.Xplicit is set and F.SlidWin is not, then F.SlidWin gets set to the
   seq length.
*/

/* process '--rule' names into a char *array and return that to main() as well */

void Flags_case_15(char* optarg, char *rulestring, struct SE_struct SelEnz[MAX_SEL_ENZ],
                   int *SECnt, char datestamp[80])
{
    long window=0;
    int   i, j, max;
    char *temp; /* temp holder for RE name in the selection string */
    char *temprule, *temprule_addr;
    char *rulelogic;
    char *RuleName="           ";
    char *str, *tok; /* another temp string */
    short tokCnt;

    FILE *fppat;
    struct stat statrec;

    /* and give some mem to the vars we need (MAX_LEN_RULESTRING=2000 at this time) */
    if ((temp = calloc(MAX_PAT_NAME_LEN+10, sizeof(char))) == NULL) {
        BadMem("Can't init mem for temp", 1);
    }
    if ((temprule = calloc(MAX_LEN_RULESTRING, sizeof(char))) == NULL) {
        BadMem("Can't init mem for temprule", 1);
    }
    temprule_addr = temprule; // have to track this separately due to strsep.
    if ((rulelogic = calloc(MAX_LEN_RULESTRING, sizeof(char))) == NULL) {
        BadMem("Can't init mem for rulelogic", 1);
    }

    if ((Rules = calloc(2, sizeof(*Rules))) == NULL) {
        BadMem("SetFlags:Rules", 1);
    }


    if (optarg == 0) BadFlag("--rule"); /* die on bad flag */
    if (F.SelEnzUsed == 1) {
        fprintf(stderr, "\nYou requested the option '--rule' but it conflicts with another flag you\n"
                "requested, probably one of -x, -p, -P, --rulefile or a previous --rule.\n"
                "Please edit the command and try again.\n");
        exit(1);
    } else {
        /* check if there's something in optarg and !using -P && !using -x */
        if (strlen(optarg) > 3 && F.Prox == 0 && F.Xplicit == 0) {
            F.Rules = 1; /* register that --rule is being used */
            F.Xplicit = 1; /* register that '--rule' is using it */
            F.SelEnzUsed = 1;  /* and be explicit about this  */

            strcpy(rulestring, optarg); /* copy optarg => rulestring to be written to the log */

            /* break optarg into separate bits and plunk subPattern names into *SelEnz[].PNames */
            /* 1st break on ',' to sep RuleName from rulestring from window  */
            RuleName = strsep(&optarg, ",");
            if (strlen(RuleName) > 10) {
                fprintf(stderr, "\nErr: in '--rule', Rule Name (%s) too long. Must be <= 10 chars\n", RuleName);
                exit(1);
            }
            /* there may need to be better format checking here to catch all the vile stuff people can do */
            /* temprule is the char * that gets chewed to bits below in place of optarg. */
            temprule = strsep(&optarg, ",");
            strcpy(rulelogic, temprule); /* copy temprule(only the rule logic) => rulestring to be passed back */
            window = atoi(strsep(&optarg, ","));
            F.SlidWin = window;
            /* at this point the rule logic hasn't been validated; only have 3 tokens to copy to Rule[] and pass on */


            /* AND..at this point, we don't know how big the sequence is going to be, so we can't check it here. */
            /* then, use logical separators (&|^) to tokenize */
            max = (int)MAX_SEL_ENZ;     /* shouldn't have to do this, but it failed otherwise */
            strcpy(temp, temprule);
            temp = DownCase(strsep(&temp, "&|^")); /* get 1st major token group */
            /* *SECnt - Selected Enz Count; NOT reset at each pass thru the case '3' stanza
                it keeps track of ALL the REs requested from all the '--rule' requests */
            while ((temp != NULL) && (*SECnt < max)) {
                /* chew off leading and lagging '()'s */
                i = 0;
                while (temp[i] == '(' || temp[i] == ')') i++;
                j = strlen(temp);
                while ((temp[j-1] == '(') || (temp[j-1] == ')')) j--;
                /* now left with i->start, j->end of 'Name:m:M' set */
                /* but need to copy it into another string if going to pass it to strsep() */
                if ((str = (char *)calloc(strlen(temp)+1, sizeof(char))) == NULL) {
                    BadMem("No mem for str", 1);
                }

                strncpy(str,temp+i,j-i);
                //free(temp); // let the orig temp go
                temp = str; /* and flip it */

                /* now break it into bits to stuff into SelEnz struct */

                /* brute-force way to check if that RE name has already been entered into SelEnz[]
                	as opposed to using realhash to check it.. whatever..*/
                i=0;

                while (i < *SECnt && strcmp(SelEnz[i].PName, temp) != 0) i++;
                if (i == *SECnt) { /* and if not, copy it in */
                    /* following maze of if / else / if's is to parse down the variants of 'Name:min:Max'
                       and all the possible combos: the std:
                       Name || Name:min:Max as well as
                       Name:min:Max,C Name=:min:Max  Name::Max  Name:min  NameWayTooLong:min:Max */

                    /* So at this pt, the token may contain the 'min' and 'Max' components (but not the
                       Hookey indicator ('=' suffixed to the name), as it's available only in the '-x'
                       flag),  so have to check if min:Max are included in the current token.
                       As coded below, *SECnt is used as the index - only useful in referring to
                       SelEnz[] but no real error checking being done here yet.  SelEnz is passed back
                       and forth and flags is global, so can extract it in main() */

                    tokCnt = 1;

                    while (tokCnt > 0 && tokCnt < 4) { /* this loop parses the token into Name, min, Max */
                        tok = strsep(&temp, ":");  /* now break the 'Name1:min:Max' into subtokens on ':' */
                        if (tok == NULL) {  /* if it's the end, this will be true */
                            if (tokCnt == 3) { /* if we're looking for the Max value and there's no token for it.. */
                                SelEnz[*SECnt].Max = 32000;
                            }
                            tokCnt = -1;    /* and the initial test will fail */
                        } else if (tokCnt == 1) { /* if it's the Name part, cp it into the struct */
                            if (strlen(tok) <= 10) {
                                strcpy(SelEnz[*SECnt].PName, tok); /* this is what we need for validation against the ReadXXX
                                                               routines - should be no need for rulestring */
                            } else {
                                fprintf(stderr, "\nThe 'Name' for one of the '--rule' options is too long.  Try again.\n");
                            }
                        } else if (tokCnt == 2) { /* looking for the 'min' value */
                            if (tok[0] == '\0') { /* but if there's no min value (ie Name1::Max), then it's the default */
                                SelEnz[*SECnt].min = 1;
                            } else if ((i = atoi(tok)) >= 0) {
                                SelEnz[*SECnt].min = i;
                            } else {
                                fprintf(stderr, "\nThe 'min' value (%d) for one of the '--rule' options is bad.  Try again.\n", i);
                                exit(1);
                            }
                        } else if (tokCnt == 3) { /* looking for the Max value */
                            if (tok[0] == '\0') { /* but if there's no Max value (ie Name1:min:), then it's the default */
                                SelEnz[*SECnt].Max = 32000;
                            } else if ((i = atoi(tok)) > 0 && i >= SelEnz[*SECnt].min) { /* allows Max = min */
                                SelEnz[*SECnt].Max = i;
                            } else {
                                fprintf(stderr, "\nThe 'Max' value (%d) for one of the '--rule' options is either < 0 or \n"
                                        "> the min value (%d).  Try again.\n", i, SelEnz[*SECnt].min);
                                exit (1);
                            }
                        }
                        tokCnt++;
                    }
                }
                (*SECnt)++;
                temp = (strsep(&temprule, "&|^"));  /* and get the next token/RE Name */
            }
            if (*SECnt >= max) {
                fprintf(stderr, "'--rule' flag has too many names; using only the 1st %d\n", max);
                exit (1);
            }

            /* and finally write the valid rulestring to the tacg.patterns file for later perusal/reuse */

            if (stat("tacg.patterns", &statrec) != 0) { /* if pat file doesn't exist, write header */
                fppat = fopen("tacg.patterns","a+"); /* open/create file to log patterns created this way */
                if (fppat == NULL) {
                    fprintf(stderr,"!!tacg.pattern open failed - que??\n");
                    exit(1);
                }
                fprintf(fppat, "Created by <tacg -p name,pattern> or <tacg -r 'name:regex_patterns'>\n"
                        "This file is a log of the patterns entered from the command line.\n"
                        "They are in REBASE format and can be edited into other such files if you wish.\n"
                        "BUT, some are in REBASE format, some are in REGEX format, and some are\n"
                        "in Rules format - careful..\n\n");
            } else { /* then the file exists, so just open the pre-existing file for appending */
                fppat = fopen("tacg.patterns","a+b");
                if (fppat == NULL) {
                    fprintf(stderr,"!!tacg.pattern open failed - que??\n");
                    exit(1);
                }
            }

            fprintf(fppat, "%s   ! Pattern added @ %s\n", rulestring, datestamp); /* and to the pattern log */


        } else {  /* whine about a badly composed enz name string */
            fprintf(stderr, "The '--rule' flag is used like this:\n"
                    "--rule 'RuleName,(NameA:m:M^NameB:m:M)|(NameC:m:M&NameD:m:M),Window'. 'man tacg' for details.\n"
                    "You may also have used the '-P' or the '-p' flags with '--rule' which isn't allowed.\n");
            exit(1);
        }
    }
    /* Finally, since everything else has been traversed OK, copy to Rule[0] for execution later on */
    strcpy(Rules[0].Name, RuleName);
    /*    if (Rules[0].RuleString = (char *)(calloc((size_t)(strlen(rulelogic))+1, sizeof(char))) == NULL) {BadMem("SetFlags:Rules[0].RuleString", 1); } */
    if ((Rules[0].RuleString = calloc(MAX_LEN_RULESTRING, sizeof(char))) == NULL) {
        BadMem("SetFlags:Rules[0].RuleString", 1);
    }
    fprintf(stdout,"ruleslogic = %s\n", rulelogic);
    strcpy(Rules[0].RuleString, rulelogic);
    Rules[0].Window = window;

//   free(temp_addr);
    free(temprule_addr);
    free(rulelogic);
}  // end of Flags_case_15()


/* ------------------------------------------------------------------------------------------ */

/* -p flag needs 3 storage parts: eg:  -p Name,gattgcnnnrygts,1
1) the name (same length as the regular names (10+\0).  As written, it creates or appends
	a REBASE style entry to the file name 'tacg.patterns', as well as a temp file that is processed
	exactly as rebase.data - this both prevents more contention for SelEnz and decreases the
	amount of supporting code and creates a record of the things you've been looking for.
2) the pattern - written to the files as described above
3) the number of ERRORs allowed in searching for these patterns.  This will cause another var to
	be writ to these files - the # of errors allowed and requires another step to detect it
	in ReadEnz(), but it should still allow most of the pre-existing code to run unmodified.
	This also allows the REBASE entries to be modified by appending this ERROR number to the data
	part of the line ie:
		hjm	   2 gcgggtswnnnnnt  0 ! anything after the '!' is ignored
		becomes
		hjm	   2 gcgggtswnnnnnt  0 2 ! where the 2 is the max number of errors tolerated.

	in ReadEnz(), the sequences corresponding to the errors are generated on the fly, checked
	against each other, and if novel, added to RE, with the mod that the (new) RE entry 'proto'
	refers all hits by this sequence to the prototype  */

void Flags_case_p(char *optarg, FILE *tmpfp, char datestamp[80])
{
    char *temp = "             "; /* temp holder for RE name in the selection string */
    char *errstr = "           "; /*  ditto for error string */
    int pCnt=0, err=0, mid;
    char *patt = "                                                                                                    "; /* for the pattern */

    FILE *fppat;
    struct stat statrec;

    if (optarg == 0) BadFlag("-p"); /* die on bad flag */
    if (F.Pat != 2 && F.Regex != 1  && F.Rules == 0) { /* go thru this only if -R, -x hasn't been set */

        if (strlen(optarg) < 3) { /* check if there's something in optarg */
            fprintf(stderr, "!! '-p' flag option is too short or missing.  The '-p' flag usage is:\n"
                    "'-p name,pattern{,error}'; repeat as needed for more patterns.\n");
            exit(1);
        }
        F.AltPatFile=2; /* warning to -R F. element to prevent possible contention */

        if (stat("tacg.patterns", &statrec) != 0) { /* if pat file doesn't exist, write header */
            fppat = fopen("tacg.patterns","a+"); /* open/create file to log patterns created this way */
            if (fppat == NULL) {
                fprintf(stderr,"!!tacg.pattern open failed - que??\n");
                exit(1);
            }
            fprintf(fppat, "Created by <tacg -p name,pattern> or <tacg -r 'name:regex_patterns'>\n"
                    "This file is a log of the patterns entered from the command line.\n"
                    "They are in REBASE format and can be edited into other such files if you wish.\n"
                    "BUT, some are in REBASE format, some are in REGEX format, and some are\n"
                    "in Rules format - careful..\n\n");
            pCnt++;
        } else if (pCnt++ == 0) {
            /* then the file exists but this is the 1st time thru, so just open the
                                        pre-existing file for appending */
            fppat = fopen("tacg.patterns","a+b");
            if (fppat == NULL) {
                fprintf(stderr,"!!tacg.pattern open failed - que??\n");
                exit(1);
            }
        } /* and if pCnt is > 0, then this is the 2nd or more time thru, so just keep adding to the file */

        /* break optarg into separate chars  */
        temp = DownCase(strsep(&optarg, ",:.")); 	/* strsep() replaces strtok() */
        patt = DownCase(strsep(&optarg, ",:.")); 	/* and works same way, regardless of # of invocations */
        if (strlen(patt)>BASE_OVERLAP) {
            fprintf(stderr, "In '-p' option, pattern must be < %d bases - truncating + continuing.\n", BASE_OVERLAP);
            patt[BASE_OVERLAP] = '\0';
        }
        mid = strlen(patt)/2;
        errstr = strsep(&optarg, ",:.");
        if (errstr == 0L) { /* then user forgot to add error term; probably wants ',0' */
            err = 0;
        } else err  = atoi(errstr); 		/* otherwise set it the usual way */
        if (err<0) {
            fprintf(stderr, "In '-p' option, error term must be >= 0 - assuming 0 and continuing.\n");
            err = 0;
        }
        if (err>MAX_ERR) {
            fprintf(stderr, "In '-p' option, error term must be < %d - assuming %d and continuing\nbut you might want to use the regex options (-r) for this magnitude of degeneracy.\n", MAX_ERR, MAX_ERR);
            err = MAX_ERR;
        }
        if (err >= strlen(patt)-1) {
            fprintf(stderr, "In '-p' option, error term must be less than (pattern length - 1).\nIf this was really intended, you're probably using the wrong application for that kind of degeneracy)\n");
            exit(0);
        }
        /* compose and print out the strings to both the temp and pattern files */
        /* date composition code went here */
        fprintf(tmpfp, "%s  %d  %s  0  %d ! Pattern added @ %s\n", temp, mid, patt, err, datestamp); /* print to temp file */
        fprintf(fppat, "%s  %d  %s  0  %d ! Pattern added @ %s\n", temp, mid, patt, err, datestamp); /* and to the pattern log */
        F.Pat = 1; /* and set the F element */
        F.Err = err; /* use F.Err to carry the error var as well */
        //F.Sites = 1; /* special case for acting like 'grep'  if you request -p, you probably want some output */
        //F.Output = 1; /* ..so set F.Output to indicate that some outout has been requested */
    } else {
        fprintf(stderr, "!! '-p' flag incompatible with already set '-R' flag; \n"
                "continuing with '-R' value.\n");
    }
}

/* ------------------------------------------------------------------------------------------ */

char * Flags_case_r(char *optarg, FILE *tmpfp, char datestamp[80])
{
    /* aka --regex='RxName:RxPattern'  or -r 'FILE:FileName' (NEEDS the single quotes)
      	so obviously you can't name a regex pattern 'FILE' */
    /* this stanza also has to duplicate a lot of the functionality of the '-p' stanza,
       as the regex's read in here have to be parsed and then written out to a temp file
       as well as a log file.  (The temp file is then read in as a standard regex file and
       all the error handling is done by the Readmatrix() - Because of the overlapping
       variables to the -p flag, they are not compatible, so have to warn each other.  */

    char *RxName, *RebFile = "\0\0", *Regex;
    char *patt = "                                                                                                    "; /* for the pattern */
    FILE *fppat;
    int pCnt=0;

    /*  struct tm *date;
    	time_t now;
    */
    struct stat statrec;     /* for file stat functions */


    if (strcmp((RxName=strsep(&optarg, ":")),"FILE") == 0) {
        RebFile = strsep(&optarg, ":");   /* next call to strsep should finish it by getting the filename */
        /* could also be that it's a mistake, but have to delay til we've exhausted the args */
        if (RebFile == NULL) {
            RebFile = "regex.data";
        }
        F.Regex = 1;  /* if patterns are read from a named file, set to 1 for -p, -R to check */
        if (stat(RebFile, &statrec) != 0) { /* if pat file doesn't exist, there's a problem. */
            fprintf(stderr, "Can't find the regex file you mentioned: %s;\n please check it and try again!\n", RebFile);
            exit(1);
        }  else { /* it at least exists, so set the file params for the pointer to be passed back */
            tmpfp = fopen(RebFile,"r"); /* open the file for reading */
            if (tmpfp == NULL) {
                fprintf(stderr,"!! %s open failed - check permissions!\n", RebFile);
                exit(1);
            }
            /*          else { fprintf(stderr, "The regex file you entered is %s;\n\n", RebFile);} */
        }

    } else if (F.Pat == 0 && F.AltPatFile == -2) {
        /*    do only if -p                -R   haven't been set */

        /* if the mode is to be used to enter regex's, the patterns are written to the temp
           and log file (the same log file as regular patterns entered with the '-p' flag.  */

        /* there's a chance that someone will try to use BOTH FILE and a named regex -should do
           error checking for this unlikely eventuality after all the rest is taken care of */

        /* warn off -p flag */
        F.Regex = -1; /* if patterns are read from the temp file, set to -1 for -p, -R to check
                         so if F.Regex == 0, don't worry about it */
        if ((Regex = calloc(200, sizeof(char))) == NULL) {
            fprintf(stderr, "calloc fails @ SetFlags:Regex\n");
            exit(1);
        }
        /* to handle naming of regex's have to pull apart optarg and tokenize at ':'
           as it'll be entered as:  -x 'RegexName:RegexPattern'
                   used as a label in RE  ------^ ^----- pattern written to tmp/log file  */
        if (optarg == NULL) {
            fprintf(stderr, "\nCheck the flag for -r or --regex.  It looks like you forgot to\n"
                    "give the regex a name; the format is:  -r 'RE_name:regex_pattern'.\n"
                    "The RE_name and the ':' ARE REQUIRED!\n");
            exit(1);
        }
        patt = DownCase(strsep(&optarg, ":"));  /* this call to strsep gets the regex pattern */

        /* OK - just write the regex's to the files */
        if (stat("tacg.patterns", &statrec) != 0) { /* if pat file doesn't exist, write header */
            fppat = fopen("tacg.patterns","a+"); /* open/create file to log patterns created this way */
            if (fppat == NULL) {
                fprintf(stderr,"!!tacg.pattern open failed - que??\n");
                exit(1);
            }
            fprintf(fppat, "Created by <tacg -p name,pattern> or <tacg -x 'name:regex_patterns'>\n"
                    "This file is a log of the patterns entered from the command line.\n"
                    "They are in REBASE format and can be edited into other such files if you wish.\n"
                    "BUT, some are in REBASE format and some are in REGEX format - careful..\n"
                    "..\n");
            pCnt++;
        } else if (pCnt++ == 0) {
            /* then the file exists but this is the 1st time thru, so just
                                     open the pre-existing file for appending */
            fppat = fopen("tacg.patterns","a+b");
            if (fppat == NULL) {
                fprintf(stderr,"!!tacg.pattern open failed - que??\n");
                exit(1);
            }
        }  /* and if pCnt is > 0, then this is the 2nd or more time thru, so just keep adding
            to the file */
        /* compose and print out the strings to both the temp and pattern files */
        /*      time(&now);                    Get the current calendar time.   */
        /*      date = localtime(&now);        Convert it to a readable date .   */
        /*      strftime(s, 80, "%c", date);    Convert that to a string.     */
        fprintf(tmpfp, "%s    %s     ! Pattern added @ %s\n", RxName, patt, datestamp); /* print to temp file */
        fprintf(fppat, "%s    %s     ! Pattern added @ %s\n", RxName, patt, datestamp); /* and to the pattern log */
    } else {
        fprintf(stderr, "!! '-p' flag incompatible with already set '-R' flag; \n"
                "continuing with '-R' value.\n");
    }
    return RebFile;
}

void Errw1(void)
{
    fprintf(stderr, "\nYou've selected 1-line output (-w 1) with a Linear (-L), ladder (-l), gel (-g),\n"
            "or proximity (-P) map, which doesn't make sense; '-w1' should only be used with\n"
            "the '-S' or '-F' flags.\n\n");
    exit(1);
}


void BadFlag(char *flag)
{
    fprintf(stderr, "\nThe %s flag either REQUIRES a value after it or can be suffixed\n"
            "with an optional value.  It therefore can't be embedded with\n"
            "single letter flags. Type 'tacg -h' for a brief description of\n"
            "the flags, or 'man tacg' for more info.\n", flag);
    exit(1);
}

long stripwhitespace(char *string, long begin, long end)
{
    long i, j, len=0;
    len = strlen(string);
    if (end == 0) {
        end = len;
    }
    if (end > len) {
        fprintf(stderr, "\nstripwhitespace: requested end (%ld) appears to be > length of string (%ld)\n"
                "will process to the end of the string\n", end, len);
        end = len;
    }
    if (begin < 0) {
        fprintf(stderr, "\nstripwhitespace: requested begin (%ld) < 0; starting from 0\n", begin);
        begin = 0;
    }
    j = i = begin;
    while (string[i] != '\0' && i <= end) {
        if (string[i] == ' ' || string[i] == '\t' || string[i] == '\n' || string[i] == '\r') {
            i++;
        } else  {
            string[j++] = string[i++];
        }
    }
    string[j] = '\0'; /* terminate */
    return (i-j);
}
