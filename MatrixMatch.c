/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright © 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: MatrixMatch.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

/* The use of this software (except that by Harald T. Alvestrand, which is described in 'udping.c')
   and that by James Knight, which is described in 'seqio.c' is bound by the notice that appears
   in the file 'tacg.h' which should accompany this file.  In the event that 'tacg.h' is not bundled
   with this file, please contact the author.
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

/* contains all the defines, includes, function prototypes for both main() and functions */
#include "tacg.h"

extern struct RE_struct *RE;
extern struct Prox_struct *PP;
extern struct Digest_Data *DD;
extern struct ORF_struct *ORFs[];

/* Below sets up the margins for running thru the seq - if it's circular (topo=0), have to
   include the repeated seqs on the beginning and end.  if topo=1, then just have to run thru
   the real seq (minus the repeats at begin end).  In both cases, have to go thru DD->Dat and
   'normalize' the position info to the actual sequence coordinates.  by subtracting
   BASE_OVERLAP from the index as written - altho this could be done on the fly without too
   much effort */

void  MatrixMatch(char *sequence, int NumREs)
{
    int ML, N, ml, FR ;
    short Success, KeepGoing;
    long CSP=0;
    int ScoreToHere;

    N = NumREs; /* delete when start doing all REs at one shot */

    for (CSP = DD->BOS; CSP < DD->EOS-4; CSP++) {   /* do the whole seq incr 1 nuc at a time */
        for (N = F.RE_ST; N < NumREs; N++) {  /* again, start at F.RE_ST, not 0 or 1 */
            Success = 0;
            KeepGoing = 1;
            ML = RE[N].E_len; /* clarity = M */
            for (FR = 0; FR <2; FR++) {  /* do Fwd and Rev Matrix at the same sequence index (but one after other) */
                ml=0;
                ScoreToHere = 0;
                KeepGoing = 1;
                while (ml < ML && KeepGoing == 1 && Success != 1) {
                    switch (sequence[CSP+ml]) {    /* a=0, c=1, g=2, t=3,  */
                    case 'a':
                        ScoreToHere += RE[N].DNAMatrix[FR][0][ml];
                        break;
                    case 'c':
                        ScoreToHere += RE[N].DNAMatrix[FR][1][ml];
                        break;
                    case 'g':
                        ScoreToHere += RE[N].DNAMatrix[FR][2][ml];
                        break;
                    case 't':
                        ScoreToHere += RE[N].DNAMatrix[FR][3][ml];
                        break;

                    case 'y':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][1][ml],
                                            RE[N].DNAMatrix[FR][3][ml]);
                        break;
                    case 'r':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][0][ml],
                                            RE[N].DNAMatrix[FR][2][ml]);
                        break;
                    case 'w':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][0][ml],
                                            RE[N].DNAMatrix[FR][3][ml]);
                        break;
                    case 's':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][1][ml],
                                            RE[N].DNAMatrix[FR][2][ml]);
                        break;
                    case 'm':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][0][ml],
                                            RE[N].DNAMatrix[FR][1][ml]);
                        break;
                    case 'k':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][2][ml],
                                            RE[N].DNAMatrix[FR][3][ml]);
                        break;

                    case 'b':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][1][ml],
                                            imax(RE[N].DNAMatrix[FR][2][ml],
                                                 RE[N].DNAMatrix[FR][3][ml]));
                        break;
                    case 'd':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][0][ml],
                                            imax(RE[N].DNAMatrix[FR][2][ml],
                                                 RE[N].DNAMatrix[FR][3][ml]));
                        break;
                    case 'h':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][0][ml],
                                            imax(RE[N].DNAMatrix[FR][1][ml],
                                                 RE[N].DNAMatrix[FR][3][ml]));
                        break;
                    case 'v':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][0][ml],
                                            imax(RE[N].DNAMatrix[FR][1][ml],
                                                 RE[N].DNAMatrix[FR][2][ml]));
                        break;

                    case 'n':
                        ScoreToHere += imax(RE[N].DNAMatrix[FR][0][ml],
                                            imax(RE[N].DNAMatrix[FR][1][ml],
                                                 imax(RE[N].DNAMatrix[FR][2][ml],
                                                      RE[N].DNAMatrix[FR][3][ml])));
                        break;

                        /* default - should be no surprises here, but just in case.. */
                    default:
                        fprintf(stderr, "What the hell is this damn char (%c) doing here?!? Ignoring..",
                                sequence[CSP+ml]);
                        break;
                    }
                    if (ScoreToHere >= RE[N].CutOff[FR]) { /* if > score, don't try any more */
                        Success = 1;
                        KeepGoing = 0;
                        /* and since IT'S A HIT, have to enter the hit info into DD.  */
                        RE[N].E_Ncuts++;      /* incr the hit conter  */
                        DD->Dat[(DD->Dat[0])++] = CSP-BASE_OVERLAP + RE[N].E_tcut[FR];   /* this may have to be offset a bit to give the correct # */
                        DD->Dat[(DD->Dat[0])++] = N;  /* and log the RE entry that generated the hit */
                    }
                    if (ScoreToHere < RE[N].ScoreMustBe[FR][ml]) KeepGoing = 0; /* if NO chance to make the match, ditto */
                    ml++;
                }
            }
        }    /* for (N = 1; N < NumREs; N++) */
        /* check need for more DD->Dat mem at the beginning of each new matrix */
        if (DD->Dat[1] - DD->Dat[0] < 50) {  /* Now check to see if we need more more memory */
            DD->Dat[1] = DD->Dat[1] + 5000;
            DD->Dat = (long *) realloc (DD->Dat, sizeof(long *) * (DD->Dat[1])); /* realloc the size 5000 bigger */
            if (DD->Dat== NULL) BadMem("MatrixMatch.c: realloc failed (DD->Dat mem)\n",1);
        }
    }
    DD->Dat[DD->Dat[0]] = -22222; /* mark the end of DD->Dat so we know where to end later on */
}


/* just returns the larger of the 2 values fed to it.  Mostly used to reduce space */
int imax(int v1, int v2)
{
    return ((v1 > v2) ? v1 : v2);
}
