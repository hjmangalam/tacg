/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: ReadRegex.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

/* The use of this software (except that by Harald T. Alvestrand, which is described in 'udping.c')
   and that by James Knight, which is described in 'seqio.c' is bound by the notice that appears
   in the file 'tacg.h' which should accompany this file.  In the event that 'tacg.h' is not bundled
   with this file, please contact the author.
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "tacg.h" /* contains all the defines, includes, function prototypes for both main() and functions */

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

int ReadRegex(char *EnzFileName, FILE *tmpfp, int *NumREs,
              struct SE_struct *SelEnz)
{

    int e=0, i, j, ii, GotH=0,
        EORE,    /* End of RE - marks and tracks the mem used by RE */
        SelREs, NumSelREs=0, RENameHash[RE_NAME_TAB_SIZ], SEi=0;

    char ct, RxName[MAX_PAT_NAME_LEN+1], RxPattern[MAXSITE], RxRest[256], *FormalRegex,
         *line;
    FILE *fpinRE;

    memset(RENameHash,0,sizeof(int)*RE_NAME_TAB_SIZ);

    SelREs = (int)abs((int)F.Xplicit); /* this way because both '-r' and '-P' (puts a -1 in F.Xplicit) use SelEnz[] */
    EORE = INIT_RES;
    if ((line = (char *) calloc(256, sizeof(char))) == NULL) {
        BadMem("No mem for ReadRegex:'line'.\n",1);
    };

    /* If the -r flag has been set, hash the names into the name space */
    if (SelREs == 1) { /* both '-r' and '-P' use SelEnz[] - should act as the correct input filter for both */
        while ((SelEnz[NumSelREs].match != 1) && (NumSelREs < MAX_SEL_ENZ)) {
            RENameHash[realhash(SelEnz[NumSelREs].PName, RE_NAME_TAB_SIZ)]++; /* incr the index */
            NumSelREs++;
        } /* and reset the '-n' and '-o' flags to default, regardless of what they were set to in SetFlags() */
        F.Mag = 3;
        F.Overhang = 1;
    }

    /* open the REbase input file, specified either by flag or by setting above */
    if (F.Regex == -1) { /* either the temp file opened and fed in SetFlags() */
        fpinRE = tmpfp; /* point one to the other */
        rewind(fpinRE); /* and rewind it */
    } else if ((fpinRE=fopen(EnzFileName,"r")) == NULL) {  /* or the standard/optional one from SetFlags() */
        fprintf(stderr,"Cannot open the REGEX file \"%s\" for reading!!\n", EnzFileName); /* DON'T HIDE behind -V */
        exit(1); /* print an error and die gracefully */
    }

    /*  Scan thru comments in rebase file to separator ("..")   */
    line[0] = 'k';   /* make sure the strcmp below fails the first time */
    while ((feof(fpinRE) == 0) && (strncmp(line, "..",2)!=0))  {
        e=fscanf (fpinRE, "%2s", line);
        if (e <= 0) {
            fprintf(stderr, "Readregex.c:74 Failure to read\n");
            exit(1);
        }
        ct = 'm'; /* and read to end of line - this needs to be able to deal */
        while (ct != '\n' && ct != '\r' && (feof(fpinRE) == 0)) {  /* with end of line conditions for Mac, Win, and unix - what have I missed? */
            ct = fgetc(fpinRE);    /* I think this is where the the lock-up happed with bad REBASE files */
        }
    }

    /* we're at "..", so go into file and start slurping directly data into struct */
    *NumREs = 1;  /* start the count at 0, to be incr as we go after each line verifies!!  */

    while ((fscanf(fpinRE, "%s %s %s", RxName, RxPattern, RxRest) != EOF)  && (*NumREs < MAX_NUM_RES )) {
        ct = 'm';
        while (ct != '\n' && ct != '\r' && (feof(fpinRE) == 0)) ct = fgetc(fpinRE); /* and read to end of line */
        if (RxName[0] != ';') {  /* if the line !commented out with a ';' */

            if (SelREs == 1) {
                SEi = MatchSelREs(RxName, SelEnz, NumSelREs);
                /* now match the SEi index to the NumREs index, so it can be special-handled in Hookey() */
                if (F.Hookey != -1) {
                    if (F.Hookey == SEi && GotH == 0) {
                        F.Hookey = *NumREs;
                        GotH = 1;
                    }
                }
            }
            /* Now match the names to those used in -P options so Proximity() knows what REs it's supposed to match */
            if (SEi != -1 && F.Prox >0) {  /* if there's a match and we need to match them for -P */
                for (i=0; i<F.Prox; i++) {  /*  check all the entries in PP[] for matches */
                    for (j=0; j<2; j++) {  /* need to check both names of each PP[] entry */
                        if (strcmp(SelEnz[SEi].PName, PP[i].N[j]) == 0) {  /* if they match entries */
                            PP[i].REi[j] = *NumREs;
                        }
                    }
                }
            }

            if (SelREs != 1 || SEi != -1) {
                if ((ii=strlen(RxName)) > 10 ) {
                    ii = 10;
                }
                strncpy(RE[*NumREs].E_nam, RxName, ii);  /* copy the name in, truncating at 10 */
                RE[*NumREs].E_nam_l  = ii; /* have the var, so enter it */
                RxName[ii] = '\0';
                SEi = -1;    /* SE index */

                if ((ii = strlen(RxPattern)) > MAXSITE) {
                    fprintf(stderr, "The REGEX pattern you specified (%s) is too long (%d chars).\n"
                            "Try again.\n", RxPattern, ii);
                    exit(10);
                }
                RE[*NumREs].E_len = ii;    /* have the var, so enter it */
                /* calloc the mem for the next few vars */
                if ((RE[*NumREs].E_raw_sit = (char *) calloc(ii+1, sizeof(char))) == NULL) {
                    BadMem("mem fails at ReadRegex:E_raw_sit\n",1);
                }
                strcpy(RE[*NumREs].E_raw_sit, RxPattern); /* copy the entered pattern into E_raw_sit (what will be used for the label */

                if ((FormalRegex = (char *)calloc((ii = (strlen(RxPattern) * 3)), sizeof(char))) == NULL) {
                    BadMem("mem fails at ReadRegex:FormalRegex\n",1);
                }

//            strcpy(FormalRegex, RxPattern);
                DNA_IUPACtoRegex(FormalRegex, RxPattern);  //now munge the entered regex into shape
		
                if ((RE[*NumREs].E_wsit[0] = (char *) calloc((ii = strlen(FormalRegex)+1) , sizeof(char))) == NULL) {
                    BadMem("mem fails at ReadRegex:E_wsit\n",1);
                }
                strcpy(RE[*NumREs].E_wsit[0], FormalRegex); /* wsit is what will be passed to the RegexMatch */
                RE[*NumREs].E_pal    = 0;
                RE[*NumREs].E_olap   = 0;
                RE[*NumREs].E_dgen   = 1;
                RE[*NumREs].WW[0]    = 0;
                RE[*NumREs].WW[1] = 0;
                RE[*NumREs].proto    = *NumREs;
                RE[*NumREs].Err      = 0;
                if (SEi >= 0) {
                    RE[*NumREs].min = SelEnz[SEi].min;  /* added 8.19.99 for SlidWin */
                    RE[*NumREs].Max = SelEnz[SEi].Max;  /* added 8.19.99 for SlidWin */
                } else {
                    RE[*NumREs].min = 1;
                    RE[*NumREs].Max = 32000;
                }

                RE[*NumREs].E_mag    = 8;   /* too hard to calculate just now */
                RE[*NumREs].EstSites = 0.0; /* and if can't calculate ", can't calc this */
                if (EORE - *NumREs < 12) {
                    if (F.Verbose > 2) fprintf(stderr, "\n\t\t#REs=%d - Need mem for RE..", *NumREs);
                    EORE += RE_INCR;
                    RE = realloc(RE, sizeof(*RE)*EORE);
                    if (RE == NULL) BadMem("Failed to get more mem for RE.\n", 1);
                    if (F.Verbose > 2) fprintf(stderr, "..Got it! RE NOW contains %d elements\n", EORE);
                }
                (*NumREs)++;
            }
        } /* end of if statement that checks for ';' at beginning of RE name */
    } /* End of while loop that reads in  RE enzyme data into struct ... whew! */

    if (SelREs==1) { /* used in both -r and -P flags; adjust the checking above to extend past 1st check */
        for (i=0; i<MAX_SEL_ENZ; i++) {   /* check for RE names used in -r flag but not matched in REBASE used */
            if ((SelEnz[i].match != 1) && (SelEnz[i].PName[0] != '\0')) { /* 1 indicates a match */
                fprintf(stderr, "Missing or Misspelled enzyme ('-x' flag): %s\n", SelEnz[i].PName);
            }
        }
    }
    return *NumREs; /* Nprotos = *NumREs in this case */
}



