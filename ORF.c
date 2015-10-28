/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright © 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: ORF.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

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
#include "tacg.h" /* contains all  the defines, includes, function prototypes for both main() and functions */


/* call Translate () on the dna sequence BEFORE call to ORF_Analysis, so to separate
   functions better and to keep from having to pass extra stuff thru functions just to pass on
   to deeper functions (don't have to pass *sequence and mode to translate if call it before
   call to this */


/*  NB:  this will still not track ORFs that cross the origin. Neat is it could... */

void ORF_Analysis(char *prot, long protlen, short thisframe)
{
    int   NeedCalloc=1,     /* indicates whether you need to calloc a new buffer for a new orf (1) or not (0) */
          protCnt=0,        /* counter that moves along the submitted protein til it hits end */
          aaCnt=0,          /* aa counter for the discovered orf */
          MinORF=F.ORFs, /* Min ORF Siz - ORF has to be at least this # of aas long or we ignore it */
          ORFNum=0,         /*  # of the orf in sequence (thisframe only) */
          EOORF = 1,        /* 1 = start out NOT being in an ORF; 0 = start out being IN an ORF */
          Structsize = 100,
          ORFstep=100, /* how much mem is calloc'ed at each incr for the orf length */
          ORF_End=100;     /* pointer to the end of ORF[][].orf */
//         ccnt = 0;

    thisframe--; /* to account for frame N referring to ORFs[N-1] */
    /* BUT ALSO have to allocate mem for the whole ORF struct and keep testing if need more mem */
    /* Get initial mem for the orf */
    ORFs[thisframe] =  calloc(Structsize, sizeof(*ORFs[0]) );
    if (ORFs[thisframe] == NULL) BadMem("ORF.c: calloc fails initial mem for ORFs!!\n",1);

    /* and for the initial length of the orf */
    ORFs[thisframe][ORFNum].orf = (char *) calloc(ORFstep, sizeof(char));
    if (ORFs[thisframe][ORFNum].orf == NULL)  BadMem("ORF.c: calloc fails @ ORFs[thisframe][ORFNum].orf!!\n",1);

    while (protCnt < (int)protlen) {
        switch (prot[protCnt]) {
        case 'M':
            if (EOORF == 1) {
                /* 1st - do we need more structs? */
                if (Structsize - ORFNum < 3) {
                    if (F.Verbose > 1) fprintf(stderr, "Getting more ORF structs ...\n");
                    Structsize += 100;
                    ORFs[thisframe] =  realloc(ORFs[thisframe], sizeof(*ORFs[0])*Structsize);
                    if (ORFs[thisframe] == NULL) BadMem("ORFs.c: calloc failed to get more mem for ORFs!!\n",1);
                }
                /* now reset all the counters and do the bookkeeping for the next ORF */
                aaCnt=0;
                ORFs[thisframe][ORFNum].frame = thisframe;

                ORF_End = ORFstep;  /* init the endpointer to the orf */
                if (thisframe < 3) {
                    ORFs[thisframe][ORFNum].B_offsetAA = (long)protCnt;   /* and the beginning offset in AAs*/
                    ORFs[thisframe][ORFNum].B_offsetBP = (long)(protCnt*3+thisframe+1);   /* and in BPs*/
                } else {
                    ORFs[thisframe][ORFNum].B_offsetAA = (long)(protlen-protCnt);   /* and the beginning offset in AAs*/
                    ORFs[thisframe][ORFNum].B_offsetBP = ((long)(protlen-protCnt)*3);   /* and in BPs*/
                }
                if (NeedCalloc == 1) {

                    /* if NeedCalloc == 0, this orf has already been calloc'ed, so no need to do it again
                       EOORF and aaCnt have already been re-zeroed, so should be no need to zero the space */
                    ORFs[thisframe][ORFNum].orf = (char *) calloc(ORFstep, sizeof(char));
                    if (ORFs[thisframe][ORFNum].orf == NULL) BadMem("ORFs.c: calloc failed to get mem for NEXT ORFs!!\n",1);
                }
            }
            ORFs[thisframe][ORFNum].orf[aaCnt++] = prot[protCnt++];  /* copy and incre */
            EOORF = 0; /* and of course, indicate that we're not at the End Of Orf any more */
            break;

        case '*':      /* at end of the orf, just note it */
            if (EOORF != 1) {  /* since we started out supposedly IN an ORF...*/
                EOORF = 1;     /* set it so it's no longer the case */
                if (aaCnt >= MinORF ) {
                    ORFs[thisframe][ORFNum].orf[aaCnt] = '\0'; /* terminate the current orf string */
                    ORFs[thisframe][ORFNum].orflength = aaCnt; /* note the length */
                    if (thisframe < 3) {
                        ORFs[thisframe][ORFNum].E_offsetAA = (long)protCnt;   /* and the ending offset */
                        ORFs[thisframe][ORFNum].E_offsetBP = (long)(protCnt*3+thisframe);   /* and the ending offset */
                    } else {
                        ORFs[thisframe][ORFNum].E_offsetAA = (long)(protlen-protCnt);   /* and the ending offset */
                        ORFs[thisframe][ORFNum].E_offsetBP = (long)((protlen-protCnt)*3+thisframe);     /* -thisframe+3;   and the ending offset */
                    }


                    ORF_Calcs(ORFNum, thisframe);
                    ORF_End = ORFstep;
                    aaCnt = 0;
                    protCnt++; /* and increment over it */
                    NeedCalloc = 1; /* and set the indicator to request more mem */
                    ORFNum++; /*  aaCnt=0; ORFs[thisframe][ORFNum].frame = thisframe;  just moved */
                } else { /* it's too short - forget it and reset counters to the last ORF */
                    EOORF = 1;
                    aaCnt = 0;
                    if (ORFNum > 0) {
                        NeedCalloc = 0;
                    }
                }
            } else {  /* then we've hit an in-frame stop NOT in an ORF soooo.. */
                protCnt++;  /* just keep incr the pointer til we get to a 'M' */
            }
            break;

        default:  /* unless we're between orfs, just add the next aa to the current orf */
            if (EOORF != 1) {  /* if not 1, then we're in an ORF, so keep adding aas to the current ORF */
                ORFs[thisframe][ORFNum].orf[aaCnt++] = prot[protCnt];
            }  /* else we're not in an ORF, so ignore all aas until we hit another 'M' */
            protCnt++; /* but keep the protCnt incrementing regardless */
            break;
        }

        /* And then grab more mem if we get close to running out of it */

        if (ORF_End - aaCnt < 2) { /* check for end of alloced mem + get more if nec */
            ORF_End += ORFstep;
            ORFs[thisframe][ORFNum].orf = (char *)realloc(ORFs[thisframe][ORFNum].orf, sizeof(char)*ORF_End);
            if (ORFs[thisframe][ORFNum].orf == NULL) BadMem("ORFs.c:~114 calloc failed to get mem for NEXT ORFs!!\n",1);
        }
    }  /* while (protCnt < (int)protlen) { ... */
    /* if there's an ORF in progress, see if it's greater than ORF and if so print it out as well.
       otherwise decr ORFNum, etc so that it will print out OK */
    if (EOORF != 1) {
        if (aaCnt >= MinORF ) { /* || ORFNum == 0  : hjm - removing the 'count the 1st ORF regardless  */
            ORFs[thisframe][ORFNum].orf[aaCnt] = '\0'; /* terminate the current orf string */
            ORFs[thisframe][ORFNum].orflength = aaCnt; /* note the length */
            ORF_Calcs(ORFNum, thisframe);
            if (thisframe < 3) {
                ORFs[thisframe][ORFNum].E_offsetAA = (long)protCnt;   /* and the ending offset */
                ORFs[thisframe][ORFNum].E_offsetBP = (long)(protCnt*3+thisframe);   /* and the ending offset */
            } else {
                ORFs[thisframe][ORFNum].E_offsetAA = (long)(protlen-protCnt);   /* and the ending offset */
                ORFs[thisframe][ORFNum].E_offsetBP = (long)((protlen-protCnt)*3+thisframe);     /* -thisframe+3;   and the ending offset */
            }
        } else {
            free(ORFs[thisframe][ORFNum].orf); /* even tho it doesn't count, have to free it */
            if (ORFNum > 0) ORFNum--;
        }
    } else { /* if we're NOT in an ORF, no changes but MIGHT have to free the last buffer (ORFNum+1) */
        if (NeedCalloc == 0) {
            free(ORFs[thisframe][ORFNum].orf);  /* also need to free the one WAITING for the */
        }                                         /* the next ORF, if NeedCalloc == 0 */
    }

    /* and now print out all the stats for thisframe, if there are any ORFs of acceptable size */

    if (ORFs[thisframe][0].orflength >= MinORF) {
        PrintORFs(ORFNum, thisframe, MinORF); /* print out all the ORFs for the entire frame */
    }
    /* ORFs memory is freed at the end of tacg.c so that the ORF data is available to other subs
    	that need it such as the grafmap thingie */
}

void PrintORFs(int ORFNum, int f, int MinORF)
{

    int i, j, k, reps, NO, seqleft; /* f=frame, just for brevity */
    int A[21] = {0,2,3,4,5,6,7,8,10,11,12,13,15,16,17,18,19,21,22,23,24};

    NO = 0;
    for (i=0; i<=ORFNum; i++) {
        if (ORFs[f][i].orflength >= MinORF) {
//         fprintf(stderr, "MinORF: %d, ORFNum: %d,  fr: %d, i: %d, cur orflength: %d \n", MinORF, ORFNum, f, i, ORFs[f][i].orflength);
            NO++; // real number of ORFs (adjusted for minlength)
        }
    }

    fprintf(stdout, "\n== ORF Analysis for Frame %d:   %d ORF(s) > %d AAs\n", f+1, NO, MinORF);  // ORFNum+1,
    fprintf(stdout, " F#  ORF#   Begin(bp/AAs)     End(bp/AAs)   #AAs    MWt(KDa)   pI\n");

    /* THIS is the only bit that has to be changed to allow wrapping the protein output
       use the same checking routing as in the '-w1' bit */
    /* STILL need to print out the absolute and % composition of the peptides. */
    for (i=0; i<=ORFNum && ORFs[f][i].orflength >= MinORF; i++) {
        fprintf(stdout, ">%2d  %4d %7ld /%6ld %7ld /%6ld   %4d %11.2f  %6.3f\n", f+1, i+1,
                ORFs[f][i].B_offsetBP, ORFs[f][i].B_offsetAA, ORFs[f][i].E_offsetBP,
                ORFs[f][i].E_offsetAA, ORFs[f][i].orflength, ORFs[f][i].MolWt, ORFs[f][i].pI);
        if (F.OrfXtra == 1) {
            fprintf(stdout, "#L     A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     X     Y\n##");
            for (k=0; k<21; k++) {
                fprintf(stdout, "  %4d", ORFs[f][i].sum_AAs[A[k]]);
            }
            fprintf(stdout, "\n#%%");
            for (k=0; k<21; k++) {
                fprintf(stdout, "%6.2f", (float)ORFs[f][i].sum_AAs[A[k]]/(float)ORFs[f][i].orflength*100);
            }
            fprintf(stdout, "\n");
        }
        /* now check to see whether the ORF gets printed as a single or wrapped line */
        if (F.Width > MAX_BASES_PER_LINE) {  /* if '-w1' (F.Width = 10000000), it gets printed as a single line.. */
            fprintf(stdout, "%s\n", ORFs[f][i].orf);  /* ..just as before */
        } else {  /* else it has to be wrapped in the following disreputable stanza */
            /* need to reference F.Width for how many AAs to put in each line
               and need to calc how many lines to put out - string should die on terminating '\0',
               regardless of what F.Width says  */
            reps = ((int)(ORFs[f][i].orflength / F.Width)) + 1;
            for (j = 0; j < reps; j++) {
                seqleft = ORFs[f][i].orflength - (int)(j*F.Width);
                if (seqleft < F.Width) {
                    fprintf(stdout, "%s\n", ORFs[f][i].orf + (int)(j*F.Width));
                } else {
                    /* note the format of the format specification below - the '.' in the "%.*s\n"
                       string SHOULD NOT be required, but it apparently is at least in the linux
                       version - have to check to see what it does in other unices */
                    fprintf(stdout, "%.*s\n", (int)F.Width, ORFs[f][i].orf + (int)(j*F.Width));
                }
            }
        }
    }
    fprintf(stdout, "\n\n");
}


void ORF_Calcs(int ORFNum, int f)   /* f = frame - just need to make it very short for lines below */
{
    /* this also should calculate pI, since now we can, and maybe do the calcs for % wt based on Sums of AAs */
    int i, j, k, AAi, AAiend, reps;
    int acids[5] = {2, 3, 4, 24, 27};  /* acids = CDEY>  for pI */
    int bases[4] = {7, 10, 17, 26};    /* bases = HKR<	for pI */
    float pIMax, pImin, pI, old_pI, result;

    /*                0     1  2    3     4     5  6  7     8  9 10    11 12    13 14 15    16 17    18    19    20 21    22 23 24    25 26   27  */
    /*                A     B  C    D     E     F  G  H     I  J  K     L  M     N  O  P     Q  R     S     T     U  V     W  X  Y     Z  <    >  */
    float pka[28]  = {0,    0, 9.0, 4.05, 4.45, 0, 0, 5.98, 0, 0, 10.0, 0, 0,    0, 0, 0,    0, 12.0, 0,    0,    0, 0,    0, 0, 10.0, 0, 7.5, 3.55 };
    float Nterm[26]= {7.59, 0, 0,   0,    7.70, 0, 0, 0,    0, 0, 0,    0, 7.00, 0, 0, 8.36, 0, 0,    6.93, 6.82, 0, 7.44, 0, 0, 0,    0 };
    float Cterm[26]= {0,    0, 0,   4.55, 4.75, 0, 0, 0,    0, 0, 0,    0, 0,    0, 0, 0,    0, 0,    0,    0,    0, 0,    0, 0, 0,    0 };

    i = ORFNum;
    j = 0; /* j = pointer to the AA under consideration */

    ORFs[f][i].MolWt = 18.0152;  /* this may have to be incr by 18 for the associated H20 */
    while (ORFs[f][i].orf[j] != '\0') {
        /* 2 lines below can be combined after it verifies */
        AAi = ((int)ORFs[f][i].orf[j]) - 65; /* 65 = decimal 'A' */
        ORFs[f][i].sum_AAs[AAi]++;   /* incr the appro AA counter */
        if (j==0) {
            if (Nterm[AAi] != 0) pka[26] = Nterm[AAi];
            /* AAiend is the index of the last AA in the ORF under considerat'n */
            AAiend = (int) ORFs[f][i].orf[ORFs[f][i].orflength-1] - 65;
            if (Cterm[AAiend] != 0) {
                pka[27] = Cterm[AAiend];
            }
            /*          if (Cterm[ORFs[f][i].orf[ORFs[f][i].orflength-1]] != 0) {
                        pka[27] = Cterm[ORFs[f][i].orf[ORFs[f][i].orflength-1]];
                     }
            */
        }

        switch (ORFs[f][i].orf[j]) { /* following - template to add other calcs in afterwards */
        case 'A':
            ORFs[f][i].MolWt +=  71.08;
            break;
        case 'C':
            ORFs[f][i].MolWt += 103.14;
            break;
        case 'D':
            ORFs[f][i].MolWt += 115.09;
            break;
        case 'E':
            ORFs[f][i].MolWt += 129.12;
            break;
        case 'F':
            ORFs[f][i].MolWt += 147.18;
            break;
        case 'G':
            ORFs[f][i].MolWt +=  57.06;
            break;
        case 'H':
            ORFs[f][i].MolWt += 137.15;
            break;
        case 'I':
            ORFs[f][i].MolWt += 113.17;
            break;
        case 'K':
            ORFs[f][i].MolWt += 128.18;
            break;
        case 'L':
            ORFs[f][i].MolWt += 113.17;
            break;
        case 'M':
            ORFs[f][i].MolWt += 131.21;
            break;
        case 'N':
            ORFs[f][i].MolWt += 114.11;
            break;
        case 'P':
            ORFs[f][i].MolWt +=  97.12;
            break;
        case 'Q':
            ORFs[f][i].MolWt += 128.14;
            break;
        case 'R':
            ORFs[f][i].MolWt += 156.2 ;
            break;
        case 'S':
            ORFs[f][i].MolWt +=  87.08;
            break;
        case 'T':
            ORFs[f][i].MolWt += 101.11;
            break;
        case 'V':
            ORFs[f][i].MolWt +=  99.14;
            break;
        case 'W':
            ORFs[f][i].MolWt += 186.21;
            break;
        case 'Y':
            ORFs[f][i].MolWt += 163.18;
            break;
        case 'X':
            ORFs[f][i].MolWt += 118.89;
            if (F.Verbose >= 1) fprintf(stderr, "UNKNOWN AA (X) in protein @ pos'n %d - CAREFUL!! Used average wt (118.89 Da)\n", j);
            break; /* need this to deal with translation of degenerate sequences */


            /* won't have a '*' as it's been changed to a \0 above and \0 shouldn't be reached */
        default:
            if (F.Verbose >= 1) fprintf(stderr, "Error (ORF_Calcs) - bad character (%c) in protein string..\n",
                                            ORFs[f][i].orf[j]);
            break;
        }
        j++;
    }  /* while (ORFs[f][i].orf[j] != '\0')  */

    /* need to give the Nterm and Cterm their due count too */
    ORFs[f][i].sum_AAs[26] = ORFs[f][i].sum_AAs[27] = 1;

    /* now all the counting and MolWt calcs are done, so now calc the pI */
    pIMax = 14;				/*  Theoretical max */
    pImin =  0;				/*  Theoretical min */
    pI    =  7;
    old_pI = 0;
    reps = 0;
    while (fabs(pI-old_pI) > 0.001) {	/*  Two correct decimals */
        if (reps++ > 15) {	/*  14/0.001 > 2^14, so this shouldn't happen */
            fprintf(stderr, "Too many reps in calc'ing pI! Continuing...\n");
        }
        result=0;
        for (k=0; k< 5; k++) { /* for the acidic part  */
            result += ORFs[f][i].sum_AAs[acids[k]] / (1+pow(10,(double)(pka[acids[k]]-pI))) ;
        }
        for (k=0; k< 4; k++) { /* for the basic part */
            result -= ORFs[f][i].sum_AAs[bases[k]] / (1+pow(10,(double)(pI-pka[bases[k]])));
        }
        old_pI = pI;
        if (result > 0) {
            pI=(pI+pImin)/2;		/* Go lower since charge is neg */
            pIMax=old_pI;
        } else {
            pI=(pIMax+pI)/2;		/*  Go higher; charge is pos */
            pImin=old_pI;
        }
    }
    /* so at this point, we should have pI for the ORF under consideration so store it! */
    ORFs[f][i].pI = pI;
}
