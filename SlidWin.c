/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com | 949 856 2847) */

/* $Id: SlidWin.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

/* The use of this software (except that by Harald T. Alvestrand, which is described in 'udping.c')
   and that by James Knight, which is described in 'seqio.c' is bound by the notice that appears
   in the file 'tacg.h' which should accompany this file.  In the event that 'tacg.h' is not bundled
   with this file, please contact the author.
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <regex.h>

/* #include "seqio.h" */

#include "tacg.h"

/***********************************  tacg_SlidingWindows  ************************************/
/* The following code searches thru DD->Dat, and verifies whether the per-pattern min/Max limits
   are being met within that window.  It returns the number of Positive Sliding Windows over all
   the patterns searched.   The struct that holds all the data is global, so it can be accessed
   from anywhere. It also handles it's own printing if requested, by setting PrintIf1 = 1 */

/* where rule = the rule to evaluate; allows fn() to evaluate a single rule at a time,
	so it can be fed indices of arbitrary rules (to allow it to be configured to eval
	only selected rules from a file of rules) */

/* this should be changed to handle the --rule as a special (single) case of the --rulefile flag;
 	then we wouldnt need to pass it rulestring - it would already have it (named) in the Rule struct. */
/* and ProtoNameHash[] shouldn;t be required either, as we don't have to do this kind of backchecking - ReadEnz
	will only read those REs that were represented in the Rule files. */
/* shouldn't need the skip hack as all the patterns need to be printed */

int tacg_SlidingWindow(int rule, int NumREs, int PrintIf1, struct SE_struct SelEnz[MAX_SEL_ENZ], int *Protos, int NumProtos)
{

    short NameCheck[RE_NAME_TAB_SIZ]; /* indices from realhash(PatName) to check against multiple identicals  */
    int *Counts=0x0, // callocable array that is indexed on # of patterns and holds the counts per pattern
         PWIncSiz = 20,  // PosWin Incr/Ending Size
         PWEndSiz = 20,
         PWCnt=0, p, EOLA, Pspace=0, r, rl, i, ii, j, k, m,
         LP, RP, SLP, Op, R, ok, L, lo, loopCnt=0,
                                        /* if gonna allow multiple rules to be checked at once,
                                         gonna have to make the 2 things following
                                         into 2D arrays ie: **LogicArr, **TruthArr */
                                        *LogicArr=0x0, *TruthArr=0x0,
                                         SWHits,	/* counter for the Sliding Window hits  */
                                         SWHBIncSiz = 200, /* the initial and incr size of SWHitBuf (carries the hits in the sliding windows) */
                                         Lparen[10] = {0,0,0,0,0,0,0,0,0,0},
                                                 Lo=2,  /* the Low end pointer to DD->Dat (bottom of the SlidWin) */
                                                 Hi=4,  /* the High end pointer to DD->Dat (top of the SlidWin) */
                                                 namelen, NameOpen, Matched, //vars for extracting the Name fr the Name:#:# string
                                                 lpi, CLVp, CTruth, RVp, COp, CRV, amt2mov, SEi, lparen_open=0;

    long seq_len, SW=0, *SWHitBuf=0x0;	/* follows DD->Dat, where [0] is current pointer; [1] is the last alloc'ed el */

    char *Name, rulematch, *tmp_str, *lcREName=0x0, *rulestring; //, *rulecopy;

    if ((rulestring = calloc(MAX_LEN_RULESTRING+2,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'rulestring'", 1);
    }
    if ((Name = calloc(MAX_PAT_NAME_LEN+2,sizeof(char))) == NULL) {
        BadMem("SlidWin.c: Name", 1);
    }
    if ((tmp_str = calloc(MAX_PAT_NAME_LEN+2,sizeof(char))) == NULL) {
        BadMem("SlidWin.c: tmp_str", 1);
    }
    seq_len = F.Seq_Len;

    /* copy in the bits we need from the Rule struct and proceed */
    strcpy(rulestring, Rules[rule].RuleString);
//   fprintf(stdout, "rulestring = %s\n", rulestring);
    strcpy(Name, Rules[rule].Name);
    if (F.SlidWin > 0) SW=F.SlidWin; // abbriev. for the Flag
    else if (Rules[rule].Window != 0) SW = Rules[rule].Window; // take it from the appro Rule struct value
    else SW = seq_len; // else it defaults to the entire sequence length */

    /* Following stanza (or subpart thereof) must be wrapped in a loop that moves the window up DD->Dat
    by hits, testing at each correctly sized window, all of the conditions.  Seems that after the
    initial error checking there should be a way to avoid going thru the whole stanza repeatedly */

    if (SW <= seq_len) {  /* 'global' if that spans the file - SW has to be smaller than whole sequence*/
        /* PW_struct is external, so don't have to return it, but only its size */
        if ((SWHitBuf = calloc(SWHBIncSiz, sizeof(long))) == NULL) BadMem("No mem for SWHitBuf", 1);
        SWHitBuf[0] = 2;
        SWHitBuf[1] = SWHBIncSiz; /* set up the initial values  */

        /* alloc PosWin to a suitable starting size */
        if ((PosWin = calloc(PWIncSiz, sizeof(*PosWin))) == NULL) BadMem("No mem for PosWin", 1);

        /* alloc PosWin[].PatHits[] (1st dimension) for the initial group */
        for (i=0; i<PWIncSiz; i++) {
            if ((PosWin[i].PatHits = calloc(NumREs, sizeof(long*))) == NULL) BadMem("No mem for PosWin[]PatHits", 1);
            /* and now calloc mem for the 0th el which keeps track of the base of the lst sliding win so we don't
            	have to count up to it each time we want to figure to how many hits there have been in the window */
            //if ((PosWin[i].PatHits[0] = calloc(NumREs, sizeof(long))) == NULL) BadMem("No mem for PosWin[].PatHits[0]", 1);

            for (j=0; j< NumREs; j++) {
                if ((PosWin[i].PatHits[j] = calloc(NumREs, sizeof(long))) == NULL) BadMem("No mem for PosWin[].PatHits[0]", 1);
            }
        }
        /* this is only for the 0th PatHits el of the 0th el of PosWin: PosWin[0].PatHits[0][...] */
        for (j=0; j< NumREs; j++) {  /* and set them all to 1 (where RE[].Sites begins */
            PosWin[0].PatHits[0][j] = 1;
        }

        rl = strlen(rulestring);

        if ((lcREName = (char *)calloc(MAX_PAT_NAME_LEN+1, sizeof(char))) == NULL) {
            BadMem("No mem for lcREName", 1);
        }
        if ((TruthArr = (int *)calloc(rl, sizeof(int))) == NULL) {
            BadMem("No mem for TruthArr", 1);
        }
        if ((LogicArr = (int *)calloc(rl, sizeof(int))) == NULL) {
            BadMem("No mem for LogicArr", 1);
        }

        /* ************************ calculate LogicArr ************************* */
        /* i = counter for rulestring
           j = counter for LogicArr  (which will be shorter) */
        /*  'bani:2:4&hinfi:3:5'  */
        i = 0;
        j = 1;
        namelen = 0;
        NameOpen = 1;
        Matched = 0;
        while (i < rl) {
            if (NameOpen == 0 && Matched == 0) { /* the name was closed by some other concluding character, so now processs the name */
                /* by running thru RE to try to match the names, unless it's already been processsed */
                k = F.RE_ST;
                m = 0;  /* k has to start @ F.RE_ST, not 1, not 0 */
                while (k < NumREs && m == 0) { /* this placement wastes 1 set of loops; maybe move it? */
                    strcpy(lcREName,RE[k].E_nam);
                    lcREName = DownCase(lcREName);
                    strcpy(tmp_str,Name);
                    tmp_str = DownCase(tmp_str);
                    if (strcmp(lcREName, tmp_str) == 0) {
                        m = 1;
                        LogicArr[j++] = k;
                    } else {
                    }
                    k++;
                }
                namelen = 0;
                Matched = 1; /* reset the Name params to allow a new name to be started */
            }

            switch (rulestring[i]) {
            default:  /* name chars - anything other than    : ( ) | &   or a  #   */
                NameOpen = 1;
                Matched = 0; /* Name is now open and hasn't been matched yet */
                /* process it as normal */
                if (namelen < MAX_PAT_NAME_LEN)  {    /* there IS a Name being composed, so if it's < Max */
                    Name[namelen++] = rulestring[i];   /* copy it over as well */
                } else {
                    fprintf(stderr, "FATAL: In '--rules', a name was too long (must be < %d chars).\n", MAX_PAT_NAME_LEN);
                    exit(1);
                }
                break;

            case '|':
                LogicArr[j++] = -1;
                Name[namelen] = '\0';
                NameOpen = 0;
                break;
            case '&':
                LogicArr[j++] = -2;
                Name[namelen] = '\0';
                NameOpen = 0;
                break;
            case '^':
                LogicArr[j++] = -7;
                Name[namelen] = '\0';
                NameOpen = 0;
                break;
            case '(':
                LogicArr[j++] = -3;
                Name[namelen] = '\0';
                NameOpen = 1;
                lparen_open=j-1;
                break;
            case ')':
                if (j == lparen_open) {
                    fprintf(stderr, "FATAL: Detected an empty paren set at ~ rulestring [%d].  Please correct!\n", i);
                    exit(1);
                }
                LogicArr[j++] = -4;
                Name[namelen] = '\0';
                NameOpen = 1;
                break;
            case ':':
                Name[namelen] = '\0';
                NameOpen = 0;
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
                /* skip it, unless NameOpen == 1 (in which case, it's part of the Name) */
                if (NameOpen == 1) {
                    Name[namelen++] = rulestring[i]; /* copy it over as well */
                }
                break;
            }
            i++;
        }
        LogicArr[0] = j; /* sets EOLA */

        /* init and zero Counts[] - will only use 1->NumREs */
        if ((Counts = calloc(NumREs+1, sizeof(int))) == NULL) BadMem("Can't init mem for Counts", 1);

        if (DD->Dat[0] > 3) { /* if > 2 hits (if [0] > 3, there has to be a paired [5]) */
            Lo = Hi = 2;
        } else {
            fprintf(stderr, "FATAL: Sliding Window: Only 1 match in whole sequence; not enough to calc sliding windows.\n");
            /* maybe need some other notes printed as well */
            exit(1);
        }
        /* in following, the Current pattern (aka RE) will be 'DD->Dat[Hi+1]' or 'DD->Dat[Lo+1]',
           a long way to represent it, but unfortunately about as short as possible */

        /* initialize the Hi pointer to just below the SW length, counting hits as we go.
           following will check for mult patterns that map to *same* site, but will NOT check for
           extremely jittered sites, if searching for extremely degenerate sites that might map
           all over a small area, so hits from degenerate patterns that occur near a boundary
           condition might not be completely accurate */

        /*      Lo = Hi = 2;  initial conditions */
        Counts[abs(DD->Dat[Lo+1])]++; /* count the 1st hit */

        while (Lo > Hi) { /* if Lo has somehow passed Hi (as would happen in crossing a BS) */
            Hi += 2; /* keep incr'g the Hi pointer until they're = */
            Counts[abs(DD->Dat[Hi+1])]++; /* & incr the counts for the new pattern */
        }

        while (((DD->Dat[Hi+2] - DD->Dat[Lo]) <= SW) && (DD->Dat[Hi+2] != -22222)) {
            Hi += 2;
            Counts[abs(DD->Dat[Hi+1])]++; /* incr the counts for the new pattern */
        } /* now both Lo and Hi are set correctly, bracketing the max Win that's <= the SW */

        /* this is where the loop for real sliding windows starts.  The LogicArray has already
        	been set up, memory has been alloc'ed for many of the structures that are needed.  Just need to set
           up the loop tests and keep looping until the end of the seq has been reached. */

        /* Next stanza is the BIG LOOP to do until the end of the array */

        /* ################################################################################## */
        SEi = 0; /* for this instance, want to restart it for each rule */
        memset(NameCheck,0, sizeof(short)*RE_NAME_TAB_SIZ); /* make sure everything is set to 0 */
        tacg_Fill_SelEnz_fr_RuleString(SelEnz, &SEi, rulestring, NameCheck, Name);

        while ((Hi <= DD->Dat[0] && DD->Dat[Hi] != -22222)) { 	/* BIG LOOP to do until the end of the array */
            /*  && (DD->Dat[Hi] - DD->Dat[Lo] <= SW) */
            loopCnt++;
            if (F.Verbose >= 2) fprintf(stderr, "Loop # %d : Window = %d : Hi = %d : Lo = %d : Diff = %ld\n", loopCnt, PWCnt, Hi, Lo, (DD->Dat[Hi] - DD->Dat[Lo]));

            /* 1st test if the window (set to be of OK size) is POSITIVE, based on rules */
            /* this will have to be redesigned to parse formal logical conditions:
               ie: ((1 & 2) | ((3 & 4)) & ((5 & 6) | 7) but it will have to do for now
               rules: (1    |     2)    &    (3    | 4)  etc  */

            /* This code eliminates the need to test for only ANDs or only ORs, as it covers both -
               I assume that the phrase contains both | and & so need to use the formal logic rules */
            /* check! will this handle flat complex strings?  ie: 'T|F&T&T|T|T&F&T|F&F&F&T&F|T&T&F' ? I think so. */
            /* 1st, reduce and serialize all the conditions, changing this:
               A&((B|C)|(D&E)&(K&G)&H)|(I&J)  (where the letter names subst for 'NameX:min:Max') to
               1&((2|3)|(4&5)&(6&7)&8)|(9&10)   done in SetFlags() and stored in passed-in 'int LogicArr[]' var
               and EOLA (EndOfLogicArray) stored as as LogicArr[0]
               EOLA points to .....................................................................v
               0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 [array index]
               1 -2 -3 -3  2 -1  3 -4 -1 -3  4 -2  5 -4 -2 -3  6 -2  7 -4 -2  8 -4 -1 -3  9 -2 10 -4 [array values]
               1  &  (  (  2  |  3  )  |  (  4  &  5  )  &  (  6  &  7  )  &  8  )  |  (  9  & 10  ) [equiv chars]
               where: positive #s = RE indices and negative #s are defined as:
                     -1 = '|' = OR     ,    -3 = '(' = LP    -5 = TRUE     -7 = '^' = XOR (exclusive OR)
                     -2 = '&' = AND    ,    -4 = ')' = RP    -6 = FALSE    -8 = reserved for later
               so all I do here is set:
               T&((F|F)|(T&T)&(T&F)&T)|(F&T)  (their equiv truth values in this routine)
               - sort of already done -> names and min:Max values are registered as RE entries in SetFlags.
               what should also be done in SetFlags is the logic string should be re-written using
               the RE indices: as above.
               so, given LogicArr as above, calculate the equiv TruthArr,
            */
            /*        p = 1;  cuz LogicArr[0] carries EOLA */
            EOLA = LogicArr[0]; /* EOLA = EndOfLogicArray, set above */

            /* Check the Min:Max values */
            /* Following sets the 'truth value' for each pattern by going thru LogicArr, checking each
               pos#  vs the RE[pos#].min/Max  */
            if (F.Verbose >= 2) {
                fprintf(stderr, "\n%s:%d:LogicArr[]:\n",__FILE__, __LINE__);
                for (i=0; i<rl; i++) {
                    fprintf(stderr, "%3d ", LogicArr[i]);
                }
                fprintf(stderr, "\n");
            }
            for (p=1; p<=EOLA; p++) {
                if (LogicArr[p] > 0) { /* for all pos #s (=RE indices) */
                    r = LogicArr[p]; /* shorten */
                    if (Counts[r] <= RE[r].Max && Counts[r] >= RE[r].min) {
                        TruthArr[p] = -5;
                    } else  TruthArr[p] = -6;
                } else { /* this should have an 'else' to allow copying over neg values that indicate logic symbols */
                    TruthArr[p] = LogicArr[p];
                }
            }
            /* so now the LogicArr has effectively been changed to TruthArr
            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29[array index]
            1  &  (  (  2  |  3  )  |  (  4  &  5  )  &  (  6  &  7  )  &  8  )  |  (  9  & 10  ) [equiv chars]
            1 -2 -3 -3  2 -1  3 -4 -1 -3  4 -2  5 -4 -2 -3  6 -2  7 -4 -2  8 -4 -1 -3  9 -2 10 -4 [logic values]
            -5 -2 -3 -3 -6 -1 -6 -4 -1 -3 -5 -2 -5 -4 -2 -3 -6 -2 -6 -4 -2 -5 -4 -1 -3 -6 -2 -5 -4 [truth values]
            T  &  (  (  F  |  F  )  |  (  T  &  T  )  &  (  F  &  F  )  &  T  )  |  (  F  &  T  ) [equiv chars]
            could these just be written right back into LogicArr? YES!, but leave for now */

            /* now need to 'flatten' the array to wipe out the parens, substituting a single TRUTH value for
               any logical expr; anything preceding a '(' is left where it is for now */

            p=1;
            LP = RP = lpi = 0;
            if (F.Verbose >= 1) {
                fprintf(stderr, "%s:%d:TruthArr:\n", __FILE__, __LINE__);
                for (i=0; i<rl; i++) {
                    fprintf(stderr, "%3d ", TruthArr[i]);
                }
                fprintf(stderr, "\n");
            }

            while (p < EOLA) {        /* (p <= EOLA) */
                r = TruthArr[p];
                switch (r) {
                case 0: /* reached at end of sweep thru the array - sweep again until only a single Truth value remains */
                    /* so test to see how much Truth has yet to be extracted , and if there is some, reset EOLA and p and
                       restart the sweep */
                    if (EOLA > 4) p = 0; /* then another sweep is required */
                    else p = EOLA + 10; /* make p > EOLA to break out of the loop */
                    break;

                default: /* there should ONLY be negative numbers here - this should never be hit with a +# */
                    /* fprintf(stderr, "%d left in place at postion %d\n", r, p); */
                    break;
                    /* only paren detection matters at this point; all else is handled inside of the parens */
                case -3: /* (, LP */
                    Lparen[LP++] = p; /* note where it lies in Lparen */
                    break;
                case -4: /* ), RP */
                    RP++; /* incr RP */
                    if (RP > LP) { /* should never be the case - if it is, err */
                        fprintf(stderr, "Too many ')'s detected (%d) in '--rules' string - try again\n", RP);
                        exit (1);
                    } else { /* eval the previously paren'ed subarray to a single TRUTH value */
                        SLP = Lparen[LP-1]; /* Starting Left Paren - where the expr has to be eval'ed FROM*/
                        /* define the Left Paren, LeftVal, Operator, RightVal, Right Paren */
                        L        = SLP; /* where the eval'ed expr starts from */
                        CLVp     = SLP+1; /* CLVp = Current Left Value pointer, set to lpi */
                        CTruth   = TruthArr[CLVp];  /* CTruth (=Current Truth) */

                        if ((p-SLP) == 2) {  /* if the paren pair just encloses a single truth value '..|(hpaii:4:9)' */
                            /* then just collapse it immediately by removing the parens */
                            TruthArr[SLP] = TruthArr[SLP+1]; /* move the Truth value over */
                            amt2mov = EOLA-p-1;
                            if (amt2mov > 0) {
                                memmove(&TruthArr[SLP+1], &TruthArr[p+1], (sizeof(int)*(EOLA-p-1)));
                                Pspace = p - Lparen[LP-1];
                                for (i=EOLA-Pspace; i<=EOLA; i++) TruthArr[i] = 0;  /* and zero the rest */
                            } else { /* just zero the rest of it */
                                for (i=p-1; i<EOLA; i++) TruthArr[i] = 0; /* and zero the rest of the array to clarify for debugging */
                            }
                            p = SLP;
                            EOLA -= Pspace;  /* EOLA = EOLA - p + 1;  correct EOLA */
                            RP--;
                            LP--; /* decr the counters, now that one pair is taken care of */
                        } else {
                            while ((RVp = CLVp+2) < p) { /* do it until all the Ops in the paren pair have been eval'ed, L -> R */
                                /*  RVp      = CLVp+2;  Right Value pointer (Right border of the expr to be eval'ed */
                                COp      = TruthArr[CLVp+1]; /* Current Operator */
                                CRV      = TruthArr[RVp]; /* Current Right Value */
                                /*  ..and CTruth is already set  */

                                if (CTruth == -5) {  /* -5 = T, -6 = F */
                                    if (COp == -1) {  /* -1 = |, -2 = & , -7 = ^  */
                                        CTruth = -5;   /*  T | X -> T */
                                    } else if (CRV == -5) {  /* if COp = '&' or '^', have to eval the last part */
                                        if (COp == -7) {  /* if CRV is T and COp is ^, then expr is F, as XOR only allows 1 truth */
                                            CTruth = -6;
                                        } else { /* if COp is '&' and the CRV is T, then expression is true */
                                            CTruth = -5; /* T & T -> T */
                                        }
                                    } else CTruth = -6; /* T & F -> F */
                                } else {  /* CTruth = F, so have to finish the eval */
                                    if (COp == -2) { /* F & X -> F*/
                                        CTruth = -6;
                                    } else if (CRV == -5) { /* if F [^|] T , have to test for which it is */
                                        if (COp == -7) { /* it's '^'  */
                                            CTruth = -5;  /* F ^ T -> T  (as below, but keep separate for now) */
                                        } else { /* it's '|' */
                                            CTruth = -5; /* F | T -> T */
                                        }
                                    } else {
                                        CTruth = -6; /* F | F -> F */
                                    }
                                    /* now finished with the 1st 2 values, reset the pointers, etc for another round */
                                }
                                CLVp = RVp;
                            } /* when leave loop, last CTruth value is value for the whole parened expr */
                            /* so drop it on top of TruthArr[SLP] (where it all began) and shift the uneval'ed
                               part of the string over to the left */
                            TruthArr[SLP] = CTruth;
                            memmove(&TruthArr[SLP+1], &TruthArr[p+1], (sizeof(int)*(EOLA-p-1)));
                            Pspace = p - Lparen[LP-1];
                            for (i=EOLA-Pspace; i<=EOLA; i++) TruthArr[i] = 0; /* and zero the rest of the array to clarify for debugging */

                            /* where ( was----^;                ^---where ) was; need 1 more to move in the new value */
                            /* EOLA = EOLA - (p+1-SLP);  correct EOLA */
                            EOLA -= Pspace; /* correct EOLA */  /* EOLA = EOLA - p + 1; */
                            p = SLP-1;	/* bc of the 'p++' at the end of the loop */
                            RP--;
                            LP--; /* decr the counters, now that one pair is taken care of */
                        }
                    }
                    break;
                } /* end of switch(r) statement */
                p++;
            }  /* end of 'while (p <= EOLA)' loop */
            /* end of the big 'else' expression that handles complex logic/truth statements */
            /* must check to see that there aren't any hanging '(' left over either.. */
            if (LP > 0) {
                fprintf(stderr, "FATAL: Problem with too many Left Parens '(' - check the rules string carefully!\n");
                exit(1);
            }
            /* now that all the parens have been elim'ed, go L -> R to calc final Truth, in this case just
               setting the Rightmost of the 2 Truth values being tested to the outcome and then starting from there
               on the next comparison */
            L = 1;
            Op = 2;
            R = 3;
            while (R <= EOLA) {  /* yes, this should be fn()ized and I will as soon as it verifies */
                if (F.Verbose >= 1) {
                    fprintf(stderr, "%s:%d:TruthArr:\n", __FILE__,__LINE__);
                    for (i=0; i<rl; i++) fprintf(stderr, "%3d ", TruthArr[i]);
                    fprintf(stderr, "\n");
                }
                if (TruthArr[Op] == -1) { /* if the Op is an '|' */
                    if (TruthArr[L] == -5 || TruthArr[R] == -5) TruthArr[R] = -5;
                    else TruthArr[R] = -6;
                } else { /* if the Op is '&' */
                    if (TruthArr[L] == -6 || TruthArr[R] == -6) TruthArr[R] = -6;
                    else TruthArr[R] = -5;
                }
                /* ie: now have the prev 'T|F' triplet replaced by 'T|T' */
                /* now set up pointers for next cmp */
                L = R;
                Op = L+1;
                R = L+2;
            }

            /* and now test the FINAL TRUTH */
            if (TruthArr[L] == -5) {
                ok = 1;
                tacg_debug(1,__FILE__,__LINE__,"--rules string evals to TRUE");
            } else {
                ok = 0;
                tacg_debug(1,__FILE__,__LINE__,"--rules string evals to FALSE");
            }

            /* and that's THAT */

            /* at end, if ok = 1, AT LEAST 1 condition has been met and Window is Pos */
            /* if in either case ok = 0, then none of the nec conditions have been met and
               the test fails -> NOT a Positive Window */

            if (ok == 1) { /* then mark the Window as Pos and set all the vars that need to be set */
                /* otherwise just store it in the struct whose pointer gets returned */
                PosWin[PWCnt].Lo = Lo;
                PosWin[PWCnt].Hi = Hi;
                PosWin[PWCnt].SeqLo = DD->Dat[Lo];
                PosWin[PWCnt].SeqHi = DD->Dat[Hi];
                /* we have to put something in PosWin[].PatHits if we expect to get anything out, right?  */

                /* Following uses Counts[] as a guide to alloc mem to hold the actual sites */

                /* #### this is why it prints thru the entire RE[] #### */
                /* #### this is now changed to match the new way of doing things #### */

                for (i=1,p=0;  (p<NumProtos && RE[Protos[p]].proto == Protos[p]); i++,p++) {
                    fprintf(stderr, "\nBEGIN: p=%d, Protos[p]=%d, RE[Protos[p]].proto=%d, RE[Protos[p]].E_nam=%s\n", p, Protos[p], RE[Protos[p]].proto, RE[Protos[p]].E_nam);
                    /* 1st, figure out how many els are needed - how many hits were in each window and then
                       alloc mem for each RE/row needed */
                    lo = PosWin[0].PatHits[0][i]; /* lo = index of the appro Sites[] array -  base of last sliding window */

                    while (RE[Protos[p]].Sites[lo] < PosWin[PWCnt].SeqLo) {
                        fprintf(stderr, "%ld:%ld -> ", RE[Protos[p]].Sites[lo], PosWin[PWCnt].SeqLo);
                        lo++;  /* just runs it up to the current base  */
                        fprintf(stderr, "%ld:%ld   ", RE[Protos[p]].Sites[lo], PosWin[PWCnt].SeqLo);
                    }
                    /* and then record this for later use so we don't have to count up from 0 each time - set BEFORE the counting
                    	loop as this base won't change after it - it's the OLD FLOOR from where we start to count up (just above) */
                    PosWin[0].PatHits[0][i] = lo; /* rememb to keep track of where the newly set, but OLD floor will be */

                    m = 2; /* SWHits = 0; */

                    /* this can be re-written to use PosWin[PWCnt].PatHits directly, instead of SWHitBuf - do the same mem jig */
                    /* have to calloc PosWin[PWCnt].PatHits[i] based on Counts[i] (+ 2 for accounting + 3 for a jitter buffer), test for overrun  */
                    /* calloc the mem for how many hits it was.. */
                    if ((PosWin[PWCnt].PatHits[i] = calloc(Counts[i]+5, sizeof(long*))) == NULL) BadMem("No mem for PosWin.PatHits[i]", 1);
                    PosWin[PWCnt].PatHits[i][1] = Counts[i]+5; /* [1] tracks the end of alloc'ed mem */
                    while (RE[Protos[p]].Sites[lo] <= PosWin[PWCnt].SeqHi && lo < RE[Protos[p]].Sites[0]
                            && RE[Protos[p]].E_Ncuts > 0 ) { /* record all the sites in here */
                        PosWin[PWCnt].PatHits[i][m++] = RE[Protos[p]].Sites[lo];
                        PosWin[PWCnt].PatHits[i][0] = m; /* [0] tracks next open position/how many */
                        SWHitBuf[0] = m;  /* add it and keep track of it */
                        lo++; /* SWHits++;   count the hits in the window by stepping lo up to the top of the window */
                        if (PosWin[PWCnt].PatHits[i][1] - PosWin[PWCnt].PatHits[i][0] < 2) { /* check the mem for overrun */
                            /* since we've approximated it pretty accurately, we shouldn't need a much bigger alloc */
                            k = PosWin[PWCnt].PatHits[i][1] + 10;
                            if ((PosWin[PWCnt].PatHits[i] = realloc(PosWin[PWCnt].PatHits[i], sizeof(long)*k)) == NULL) BadMem("Can't realloc mem for SWHitBuf", 1);
                            PosWin[PWCnt].PatHits[i][1] = k;
                        }
                    }
                    /* RE[Protos[p]].Sites[lo] ends being <= to PosWin[PWCnt].SeqHi && lo < RE[Protos[p]].Sites[0]  */
                    /* PosWin[0].PatHits[0][i] (below) holds the RE[].Sites position of the floor (using lo as the index)
                    	so we don't have to count up each time */
                }
                ii = i;
                if (F.Verbose >= 2) {
                    for (p=0; p< PWCnt; p++) {
                        fprintf(stderr,"PWCnt = %d ---\n", p);
                        for (i=0; i<ii; i++) {
                            fprintf(stderr,"%d (%s); ", i, RE[Protos[p]].E_nam);
                            for (m=0; m< PosWin[p].PatHits[i][1]; m++) {
                                fprintf(stderr, "xxx%ld ", PosWin[p].PatHits[i][m]);
                            }
                            fprintf(stderr,"\n");
                        }
                        fprintf(stderr,"\n");
                    }
                    fprintf(stderr,"\n");
                }

                /* and check to see if we're too close to the end of mem for PosWin */
                if ((PWEndSiz - PWCnt) < 2) {
                    m = PWEndSiz; /* now it's the new start site for alloc'ing new memory */
                    PWEndSiz += PWIncSiz;
                    if ((PosWin = realloc(PosWin, sizeof(*PosWin)*PWEndSiz)) == NULL) BadMem("Can't realloc mem for PosWin", 1);
                    /* now set all the new vars and calloc new mem for all of them   */
                    for (i=m; i<PWEndSiz; i++) {
                        if ((PosWin[i].PatHits = calloc(NumREs, sizeof(long*))) == NULL) BadMem("No mem for new PosWin", 1);
                    }
                }
                /* if output should be printed as well as stored in struct, print the header only once */

                if (PrintIf1 == 1 && PWCnt == 0) {

                    fprintf(stdout, "\n== Positive Sliding Windows for %s (Window Size = %ld) \n\n",
                            Rules[rule].Name, Rules[rule].Window);
                    /* and the warning that there may be some funny biz with degenerate sequences */
                    fprintf(stdout, "NB: When using very degenerate sequences (ie: 'wwwwwwwwww') in conjunction\n"
                            "with Sliding Windows,  inaccuracies in the way that the hits are counted\n"
                            "may allow hits to be OVERCOUNTED.  If you suspect that this is happening,\n"
                            "check the Sites table (-S) against the Sliding Windows.  Sites data has\n"
                            "been  filtered to remove this 'jitter'.  \n"
                            "Each stanza includes header info about the window and one line on each pattern:\n\n"
                            "Serial#    Low Border    High Border       Diff\n"
                            "   \"  Name:min:Max   (# hits in window)   Pattern sequence   Actual sites....\n\n");
                }

                /* printing section also needs to check current rule and ONLY print that pattern data for the the patterns
                in that particular rule (not all patterns in the rule file as is currently the case */

                if (PrintIf1 == 1) { /* and if this should be printed as well as stored in struct */
                    /*               Serial #   Low Border   High Border    Difference     */
                    fprintf(stdout, "%d Rule Name = [%s] Req'd Win = [%ld] Lo = [%ld]  Hi = [%ld]  Actual Win = [%ld]  \n%d Logic = %s\n",
                            PWCnt, Rules[rule].Name, Rules[rule].Window, PosWin[PWCnt].SeqLo, PosWin[PWCnt].SeqHi,
                            (PosWin[PWCnt].SeqHi - PosWin[PWCnt].SeqLo), PWCnt, rulestring);


//  This needs to check the protos, not the RE stuct directly!!

                    for (i=1,p=0; p < NumProtos && RE[Protos[p]].proto == Protos[p] ; i++,p++) {
                        fprintf(stdout, "\nLooping i=%d p=%d NumProtos = %d, Protos[p] = %d, RE[Protos[p]].proto=%d, RE[Protos[p]].E_nam = %s\n", i, p, NumProtos, Protos[p], RE[Protos[p]].proto, RE[Protos[p]].E_nam);
                        //P = Protos[p];
//               for (i=F.RE_ST; i<NumREs  && RE[i].proto == i; i++) { /* start @ F.RE_ST, not 1, not 0 */
                        if (PosWin[PWCnt].PatHits[i][0] == 0) SWHits = 0;
                        else SWHits = PosWin[PWCnt].PatHits[i][0]-2;

                        if (NameCheck[realhash(DownCase(RE[Protos[p]].E_nam), RE_NAME_TAB_SIZ)] == 1) {
                            if (SWHits >= RE[Protos[p]].min && SWHits <= RE[Protos[p]].Max) {
                                rulematch = 'x';
                            } else rulematch = ' ';
                            fprintf(stdout, "%d %*s:%d:%d %c (%d)  %s   ",
                                    PWCnt, MAX_PAT_NAME_LEN, RE[Protos[p]].E_nam, RE[Protos[p]].min, RE[Protos[p]].Max, rulematch, SWHits, RE[Protos[p]].E_raw_sit);
                            /* following is suposed to print out the actual hits but it's currently buggered
                            	it may not be a good idea to reiterate this data for every window as that's a lot of
                               repeated data.  Better to just use the Sites output? */
                            for (j=2; j< PosWin[PWCnt].PatHits[i][0]; j++) {
                                fprintf(stdout, " %ld  ", PosWin[PWCnt].PatHits[i][j]);
                            }
                            fprintf(stdout, "\n");
                        }
                    }
                    fprintf(stdout, "\n");
                }
                PWCnt++;
                fprintf(stderr, "PWCnt now %d\n",PWCnt );
            }  /*   if (ok == 1) block */
            /* Now move the SW up.  Lo always moves up one position.  Hi may or may not,
            depending on whether the next hit is too far away.  If Lo becomes > Hi, it
            indicates that we've jumped a BS and Hi has to be run up to be >= to Lo.   */
            Counts[abs(DD->Dat[Lo+1])]--;  /* decr the count for the current (old) hit */
            Lo += 2;                       /* incr Lo to move up DD->Dat */
            while (Lo >= Hi && (DD->Dat[Hi + 2] - DD->Dat[Lo] <= SW)) { /* if Lo has somehow passed Hi (as would happen in crossing a BS) */
                Hi += 2; /* keep incr'g the Hi pointer until they're = */
                Counts[abs(DD->Dat[Hi+1])]++; /* & incr the counts for the new pattern */
            }
            while (((DD->Dat[Hi+2] - DD->Dat[Lo]) <= SW) && (DD->Dat[Hi] != -22222)) {
                //(DD->Dat[Hi+2] != -22222)){ - incredibly stupid - tested wrong el - off by 2
                // retain as warning to self about stupid incrementing.
                Hi += 2;
                Counts[abs(DD->Dat[Hi+1])]++; /* incr the counts for the new pattern */
            } /* now both Lo and Hi are set correctly, bracketing the max Win that's <= the SW */
        }  /* while (Hi <= DD->Dat[0] || (( Hi - Lo == 2) && (DD->Dat[Hi] - DD->Dat[Lo] > SW)))  */
    } else {		/* if (SW <= (int)seq_len)  - file spanning if loop */
        fprintf(stderr, "Sliding Window (%ld) > entire sequence (%ld), skipping Sliding Window analysis.\n", SW, seq_len);
    }
    free(rulestring);
    free(PosWin);
    free(LogicArr);
    free(TruthArr);
    free(Counts);
    free(Name);
//   free(tmp_str_addr);
    free(lcREName);
    free(SWHitBuf);
    for (i=0; i<PWIncSiz; i++) {
        free(PosWin[i].PatHits[0]);
        free(PosWin[i].PatHits);
    }


    return PWCnt; /* PosWin is already available as it's global */
}
