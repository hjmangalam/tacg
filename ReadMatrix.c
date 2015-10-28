/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: ReadMatrix.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

/* The use of this software (except that by Harald T. Alvestrand, which is described in 'udping.c')
   and that by James Knight, which is described in 'seqio.c' is bound by the notice that appears
   in the file 'tacg.h' which should accompany this file.  In the event that 'tacg.h' is not bundled
   with this file, please contact the author.
*/

/* ReadMatrix does just that - it takes a TRANSFAC formatted file of matrices and extracts a
   number of infobits, loading those that need to be into the global RE struct.  It also calls
   CreateConsensus() which also does as it suggests, taking the matrix and creating a consensus
   out of IUPAC degeneracies (as opposed to simply the 4 bases, as is given by the Matrix file).
   This consensus is used both for labeling and for searching via tacg's more direct search
   algorithm.  The consensus is associated with a cutoff value to choose the degeneracy that best
   represents the profile at that point.  */

/* RebaseFile/MatrixFile can be either the default or an optional filename - decided in/immed after
   SetFlags() */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "tacg.h"

int ReadMatrix(char *EnzFileName, int *NumREs, struct SE_struct *SelEnz)
{

    int e=0, i, ii, j, jj, k, l, KeepGoing, tmax, Order,
        which, /* some of these counters can prob be re-used */
        line = 0,
        EORE,          /* End of RE - marks and tracks the mem used by RE */
        ConsCutOff,    /* the cutoff for calc'g the Consensus, entered on the command line as the
                        arg to -#, also used for the cutoff for MatrixMatch() */
        len,
        SelREs, NumSelREs=0, RENameHash[RE_NAME_TAB_SIZ], SEi=0;

    float A, C, G, T;
    int Pass, Accession=0,
              DMat[15],               /*   a c g t   w s y r m k    b  d  h  v    n  (Degeneracy Matrix */
              /*   0 1 2 3   4 5 6 7 8 9   10 11 12 13   14  */
              MaxPossibleScore[2],    /* for each matrix */
              MatLen,                 /* Matrix Length */
              NumSites,               /* # of members of the Matrix */
              TmpNumSites,            /* tmp ditto */
              HiScore[MAX_RAW_MATRIX_LEN],     /* highest score possible for that point - used to calc ScoreMustBe */
              RMat[2][4][MAX_RAW_MATRIX_LEN];  /* holds the raw counts, fractions from the Transfac Matrix file */

    char ct,  Name[11],
         *ID,   /* Holds the IDentification word from Matrix file that IDs the data on that line */
         TransfacConsensus[MAX_RAW_MATRIX_LEN], /* holds the Transfac-supplied consensus (which may or may
                                                not be the one calculated by own calculations)  */
         Degeneracies[] = { 'a','c','g','t','w','s','y','r','m','k','b','d','h','v','n' },
                          MyConsensus[MAX_RAW_MATRIX_LEN]; /* the consensus calc'ed here, from the read-in matrix */

    FILE *fpinRE;
    memset(RENameHash,0,sizeof(int)*RE_NAME_TAB_SIZ);

    ConsCutOff = (int)F.Matrix * 10;
    SelREs = (int)abs((int)F.Xplicit); /* this way bc both '-r' and '-P' (puts a -1 in F.Xplicit) use SelEnz[] */
    *NumREs=1;  /* this has to start at 1 rather than 0, because negation is used later to mark the doppels */
    EORE = INIT_RES;
    ID = (char *) calloc(256, sizeof(char));

    /* If the -r flag has been set, hash the names into the name space */
    if (SelREs == 1) { /* both '-r' and '-P' use SelEnz[] - should act as the correct input filter for both */
        while ((SelEnz[NumSelREs].PName[0] != '\0') && (NumSelREs < MAX_SEL_ENZ)) {
            RENameHash[realhash(SelEnz[NumSelREs].PName, RE_NAME_TAB_SIZ)]++; /* incr the index */
            NumSelREs++;
        } /* and reset the '-n' and '-o' flags to default, regardless of what they were set to in SetFlags() */
        F.Mag = 3;
        F.Overhang = 1;
    }

    /* open the Matrix input file, specified either by flag or by setting above */
    if ((fpinRE=fopen(EnzFileName,"r")) == NULL) {  /* or the standard/optional one from SetFlags() */
        fprintf(stderr,"Cannot find or open the Matrix file \"%s\" for reading!!\n", EnzFileName); /* DON'T HIDE behind -V */
        exit(1); /* print an error and die gracefully */
    }

    /*  Scan thru comments in rebase file to separator ("//")   */
    ID[0] = 'k';   /* make sure the strcmp below fails the first time */
    while ((feof(fpinRE) == 0) && (strncmp(ID, "//",2)!=0))  {
        e=fscanf (fpinRE, "%2s", ID);
        if (e <= 0) {
            fprintf(stderr, "ReadMatrix.c:92 Failure to read\n");
            exit(1);
        }

        ct = 'm'; /* and read to end of line - this needs to be able to deal */
        line++; /* for debugging */
        while (ct != '\n' && ct != '\r' && (feof(fpinRE) == 0)) {  /* with end of line conditions for Mac, Win, and unix - what have I missed? */
            ct = fgetc(fpinRE);    /* I think this is where the the lock-up happed with bad REBASE files */
        }
    }

    /* we're at "//", so go into file and start slurping directly data into struct */
    *NumREs = 1;  /* start the count at 1, not 0!!  */

    while ((fgets(ID, 256, fpinRE) != NULL)  && (*NumREs < MAX_NUM_RES )) {
        line++;
        if (ID[0] != ';' && ID[0] != 'X') {  /* if ! commented out or it isn't a separator line (XX) */
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

            if (strstr(ID, "NA") != NULL) {  /* if it's a NA (Name) line  */
                i=4;
                len=0;
                while (ID[i] != '\000' && i < 10 && ID[i] != '\n') {  /* while not at end of string and i<10 (length of array) */
                    if (ID[i] != '_' && ID[i] != '\'' && ID[i] != ' ') { /*As long as it's a valid char */
                        Name[len++] = ID[i++]; /* The struct value gets the whole site */
                    } else i++;  /* incr i only, skipping the _ and ' chars  */
                }
                Name[len] = '\0';
                SEi = -1;    /* SE index */
                if (SelREs == 1) {
                    SEi = MatchSelREs(Name, SelEnz, NumSelREs);
                }
            }  else if (strstr(ID, "AC") != NULL) {   /* Is it the TRANSFAC Accession # ? */
                ID[10] = '\0';
                Accession = atoi(ID+5);     /* copy it into RE */
            }  else if ((strstr(ID, "P0") != NULL) && (SelREs != 1 || SEi != -1)) {  /* header just before the matrix comes in */

                strcpy(RE[*NumREs].E_nam, Name);
                RE[*NumREs].M_Acc = Accession;
                j = 0;
                while ((fgets(ID, 256, fpinRE) != NULL) && (*NumREs < MAX_NUM_RES && ID[0] != 'X')
                        && (j < MAX_RAW_MATRIX_LEN)) {
                    sscanf(ID, "%d %f %f %f %f  %c", &Order, &A, &C, &G, &T, &TransfacConsensus[j]);
                    A *= 100;
                    C *= 100;
                    G *= 100;
                    T *= 100;
                    RMat[0][0][j] = (int)A;
                    RMat[0][1][j] = (int)C;
                    RMat[0][2][j] = (int)G;
                    RMat[0][3][j] = (int)T;
                    j++;
                }
                TransfacConsensus[j] = '\0'; /* terminate it so it doesn't cause catastrophy */
                RE[*NumREs].E_len = MatLen = j;  /* store this in RE as well for other fn()s */
                for (i=0; i<2; i++) RE[*NumREs].E_tcut[i] = (int) MatLen/2;
                RE[*NumREs].E_olap = 0;

                /* Raw Matrix now read in  - now process it down to an RE entry  */
                /* add up the entries to figger how many sites were used to calc it or rather what the
                   sum is across the rows */
                NumSites = 0;
                for (i=0; i<MatLen; i++) {
                    TmpNumSites = 0;
                    for (j=0; j<4; j++) {
                        TmpNumSites += RMat[0][j][i];
                    }
                    if (NumSites < TmpNumSites) NumSites = TmpNumSites; /* should all be =, but if not, take biggest */
                }
                /* convert forward matrix to reverse complement matrix */
                for (i=0, ii=MatLen-1; i<MatLen; i++, ii--) {
                    for (j=0, jj=3; j<4; j++, jj--) {
                        RMat[1][jj][ii] = RMat[0][j][i]; /* does a reverse diagonal copy from [0] to [1]  */
                    }
                } /* this should now be ready to copy into RE in the Matrix part */

                /* ok - now know how many sites were used to generate the cons; now convert to fractions;
                   calculate the REAL consensus; also have to convert to reverse complement fractions */
                /* these numbers will have little relation to the numbers read in because of the manipulation
                   I've done, but the outcome should be OK */
                for (k=0; k<2; k++) {
                    MaxPossibleScore[k] = 0;
                    /* alloc the mem in RE[].ScoreMustBe  at same time, saving a little 'for'*/
                    RE[*NumREs].ScoreMustBe[k] = (int *) calloc(MatLen, sizeof(int)); /* get the mem */
                    if (RE[*NumREs].ScoreMustBe[k] == NULL) BadMem("calloc failed - ScoreMustBe[i] in ReadMatrix.c", 1);
                    for (l=0; l <= MatLen; l++) HiScore[l] = 0;  /* reset HiScore */
                    for (i=MatLen-1; i>=0; i--) {
                        for (j=0; j<4; j++) {
                            RMat[k][j][i] = RMat[k][j][i] * 1000 / NumSites; /* jacking it up to avoid floats */
                            HiScore[i] = ((HiScore[i] > RMat[k][j][i]) ? HiScore[i] : RMat[k][j][i]); /* keep track of the highest score */
                        }
                        MaxPossibleScore[k] += HiScore[i];
                    }

                    /* calculate the RE[].ScoreMustBe[] values here since know everything that's required.
                       This value is the value at each el that must be equalled or exceeded or the match will
                       fail - allows early failure to a pattern.  While I'm thinking about it, could also try to
                       set up an array which will help to fail early by taking the top scoring els and test them
                       1st.  That way, should only have to test the 1st few and it will fail, rather than always
                       delaying the failure to near the end.  Could arrange it with an array of pointers to the
                       els, so just have to walk down the primary array in order which point to the out-of-order
                       top-scoring els until it fails */

                    /* can do it this way because the array RMat[][][] was rev copied above so */
                    /* HiScore must be * by the cutoff score -otherwise only BEST matrices would match */
                    Pass = (int)((float)MaxPossibleScore[k] * F.Matrix / 100);
                    for (i=MatLen-1; i>=0; i--) {
                        RE[*NumREs].ScoreMustBe[k][i] = imax(Pass - HiScore[i] ,0);
                        Pass = RE[*NumREs].ScoreMustBe[k][i];
                    }
                    RE[*NumREs].MaxScore[k] = MaxPossibleScore[k];
                    RE[*NumREs].CutOff[k] = MaxPossibleScore[k] * F.Matrix / 100;
                }
                /* and copy these RMat values into RE[].DNAMatrix after alloc'ing the correct mem */
                /* float *DNAMatrix[2][4] */
                for (i=0; i<2; i++) {
                    for (j=0; j<4; j++) {
                        RE[*NumREs].DNAMatrix[i][j] = (int *) calloc(MatLen, sizeof(int)); /* get the mem */
                        if (RE[*NumREs].DNAMatrix[i][j] == NULL) BadMem("calloc failed - DNAMatrix[i][j] in ReadMatrix.c", 1);
                        RE[*NumREs].E_raw_sit = (char *) calloc(MatLen+1, sizeof(char));
                        if (RE[*NumREs].E_raw_sit == NULL) BadMem("calloc failed - RE[*NumREs].E_raw_sit in ReadMatrix.c", 1);
                        for (k=0; k<MatLen; k++) {
                            RE[*NumREs].DNAMatrix[i][j][k] = RMat[i][j][k];  /* and copy the values in */
                            RE[*NumREs].E_raw_sit[k] = TransfacConsensus[k];
                        }
                    }
                }
                /* now calculate the real consensus, using the user-supplied cutoff from SetFlags()  */
                /*      DegenMatrix aka DMat[15],   a c g t   w s y r m k    b  d  h  v    n  */
                /*                                  0 1 2 3   4 5 6 7 8 9   10 11 12 13   14  */

                for (i=0; i<MatLen; i++) { /* at each char for the length of the matrix */
                    for (l=0; l<15; l++) DMat[l] = 0; /* zero DMat, the degeneracy counter */
                    for (j=0; j<4; j++) {  /* for each base position that has been read in */
                        if (RMat[0][j][i] > 0) {  /* if there was at least 1 base of this type */
                            ii = RMat[0][j][i];    /* for clarity/brevity */
                            if (j == 0) {  /* if it was an 'a' */
                                DMat[0] = ii;
                                DMat[4] += ii;
                                DMat[7] += ii;
                                DMat[8] += ii;
                                DMat[11] += ii;
                                DMat[12] += ii;
                                DMat[13] += ii;
                                DMat[14] += ii;
                            } else if (j == 1) {  /* if it was a 'c' */
                                DMat[1] = ii;
                                DMat[5] += ii;
                                DMat[6] += ii;
                                DMat[8] += ii;
                                DMat[10] += ii;
                                DMat[12] += ii;
                                DMat[13] += ii;
                                DMat[14] += ii;
                            } else if (j == 2) {  /* if it was a 'g' */
                                DMat[2] = ii;
                                DMat[5] += ii;
                                DMat[7] += ii;
                                DMat[9] += ii;
                                DMat[10] += ii;
                                DMat[11] += ii;
                                DMat[13] += ii;
                                DMat[14] += ii;
                            } else if (j == 3) {  /* if it was a 't' */
                                DMat[3] = ii;
                                DMat[4] += ii;
                                DMat[6] += ii;
                                DMat[9] += ii;
                                DMat[10] += ii;
                                DMat[11] += ii;
                                DMat[12] += ii;
                                DMat[14] += ii;
                            }
                        }
                    } /* finished loading the DegenMatrix; now determine what the consensus char should be */

                    /* this should be stratified so that it evals the 3 levels of degeneracy at each step,
                    taking the best estimate at each level of degeneracy; otehrwise it's going to end up
                    possibly taking passing, but suboptimal levels at each point ..  but leave it for now...  */

                    /* check 1x degeneracies - acgt */
                    KeepGoing = 1; /* set the flag to Keep Going */
                    for (k=0; k<4 && KeepGoing == 1; k++) { /* check the acgt's 1st */
                        if (DMat[k] >= ConsCutOff) { /* 1st one > cutoff wins here */
                            MyConsensus[i] = Degeneracies[k];
                            KeepGoing = 0; /* and let the rest know that we got it */
                        }
                    }
                    /* check 2x degeneracies - wsyrmk */
                    if (KeepGoing == 1) { /* if didn't get it on the acgt level */
                        tmax = 0;
                        which = 0;
                        for (k=4; k<10; k++) { /* check ALL the 2x degens - wsyrmk */
                            if (DMat[k] > tmax && DMat[k] >= ConsCutOff) {
                                tmax = DMat[k];
                                which = k;
                            }
                        }

                        if (tmax > 0) { /* and only after checking them ALL ... */
                            MyConsensus[i] = Degeneracies[which]; /* choose the top score if it exceeded the cutoff */
                            KeepGoing = 0; /* and tell that there's no point in going on */
                        }
                    }
                    /* check 3x degeneracies - bdhv */
                    if (KeepGoing == 1) { /* but if still didn't get it... */
                        tmax = 0;
                        which = 0;
                        for (k=10; k<14; k++) {   /* check ALL the 3x degens - bdhv */
                            if (DMat[k] > tmax && DMat[k] >= ConsCutOff) {
                                tmax = DMat[k];
                                which = k;
                            }
                        }
                        if (tmax > 0) {
                            MyConsensus[i] = Degeneracies[which];  /* ditto - choose only TOP score */
                            KeepGoing = 0;
                        }
                    }
                    /* and finally, if everything else failed, take the default 'n' */
                    if (KeepGoing == 1) {
                        MyConsensus[i] = 'n';
                    }
                }

                MyConsensus[i] = '\0'; /* terminate it correctly */
                /* now should have the full self-calc'ed consensus so copy it into RE[].E_wsit. */
                /* so grab mem..*/
                RE[*NumREs].E_wsit[0] = (char *)calloc(MatLen+1, sizeof(char));
                if (RE[*NumREs].E_wsit[0] == NULL) BadMem("calloc - RE[*NumREs].E_wsit[0]", 1);
                RE[*NumREs].E_wsit[1] = (char *)calloc(MatLen+1, sizeof(char));
                if (RE[*NumREs].E_wsit[1] == NULL) BadMem("calloc - RE[*NumREs].E_wsit[1]", 1);
                /* and copy values */
                strcpy(RE[*NumREs].E_wsit[0],MyConsensus);
                Anti_Par(MyConsensus, RE[*NumREs].E_wsit[1], MatLen);
                RE[*NumREs].proto = *NumREs;
                RE[*NumREs].E_mag = MagCalc(RE[*NumREs].E_wsit[0], RE[*NumREs].E_len);
                RE[*NumREs].E_nam_l = strlen(RE[*NumREs].E_nam);
                RE[*NumREs].min = SelEnz[SEi].min; /* added 8.19.99 for SlidWin */
                RE[*NumREs].Max = SelEnz[SEi].Max; /* added 8.19.99 for SlidWin */

                if (F.Verbose > 1) {  /* if you really want to know... */
                    fprintf(stderr, "Matrix for %s (# %d) is(len=%d):\n", RE[*NumREs].E_nam, *NumREs, RE[*NumREs].E_len);
                    for (k=0; k<2; k++) {
                        fprintf(stderr, "FR=%d\n", k);
                        fprintf(stderr, "El#");
                        for (i=0; i<MatLen; i++) {
                            fprintf(stderr, " %4d", i);
                        }
                        fprintf(stderr, "\n");
                        for (j=0; j<4; j++) {
                            fprintf(stderr, "%d  ", j);
                            for (i=0; i<MatLen; i++) {
                                fprintf(stderr, " %4d", RE[*NumREs].DNAMatrix[k][j][i]);
                            }
                            fprintf(stderr, "\n");
                        }
                        fprintf(stderr, "\nMaxPossible=%d, Cutoff=%d   ScoreMustBe:\n   ", RE[*NumREs].MaxScore[k], RE[*NumREs].CutOff[k]);
                        for (i=0; i<MatLen; i++) {
                            fprintf(stderr, " %4d", RE[*NumREs].ScoreMustBe[k][i]);
                        }
                        fprintf(stderr, "\n\n");
                    }
                }

                (*NumREs)++; /* incr the counter */
                /* and check that we're not going to run into th eend of the array */
                if (EORE - *NumREs < 12) {
                    if (F.Verbose > 1) fprintf(stderr, "\nNeed mem for RE (matrices)..");
                    EORE += RE_INCR;
                    RE = realloc(RE, sizeof(*RE)*EORE);
                    if (RE == NULL) BadMem("Failed to get more mem for RE.\n", 1);
                    if (F.Verbose > 1) fprintf(stderr, "..Got it! RE NOW contains %d elements\n", EORE);
                }

                /* this next bit should be functionized - called 2 times so far */
                /*  Scan thru comments in rebase file to separator ("//") again */
                ID[0] = 'k';   /* make sure the strcmp below fails the first time */
                while ((feof(fpinRE) == 0) && (strncmp(ID, "//",2)!=0))  {
                    e =fscanf (fpinRE, "%2s", ID);
                    /* if (e <= 0) {fprintf(stderr, "ReadMatrix.c:379 Failure to read\n"); exit(1); } */
                    ct = 'm'; /* and read to end of line - this needs to be able to deal */
                    line++; /* for debugging */
                    while (ct != '\n' && ct != '\r' && (feof(fpinRE) == 0)) {  /* with end of line conditions for Mac, Win, and unix - what have I missed? */
                        ct = fgetc(fpinRE);    /* I think this is where the the lock-up happed with bad REBASE files */
                    }
                } /* we're at "//" again, so go into file and start slurping directly data into struct */
            }
        } /* end of if statement that checks for ';' at beginning of RE name */
    } /* End of while loop that reads in  RE enzyme data into struct ... whew! */

    if (SelREs==1) { /* used in both -r and -P flags; adjust the checking above to extend past 1st check */
        for (i=0; i<MAX_SEL_ENZ; i++) {   /* check for RE names used in -r flag but not matched in REBASE used */
            if ((SelEnz[i].match != 1) && (SelEnz[i].PName[0] != '\0')) { /* match = 1 indicates a match */
                fprintf(stderr, "Missing or Misspelled enzyme ('-x' flag): %s\n", SelEnz[i].PName);
            }
        }
    }
    /*    fclose(tmpfp);  Still going to need this tmpfp for trnasfering RE names? */
    return *NumREs; /* Nprotos = *NumREs in this case */
}



