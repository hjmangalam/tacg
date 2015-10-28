/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
   (hjm@tacgi.com, 949 856 2847) */

/* $Id: Proximity.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

/* The use of this software (except that by Harald T. Alvestrand, which is described in 'udping.c')
   and that by James Knight, which is described in 'seqio.c' is bound by the notice that appears
   in the file 'tacg.h' which should accompany this file.  In the event that 'tacg.h' is not bundled
   with this file, please contact the author.
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
/* #include <math.h> */

/* contains all the defines, includes, function prototypes for both main() and functions */
#include "tacg.h"

void Proximity(void)
{
    /* need to go thru all the PP struct entries, doing all the matching between
       pairs, writing the matches of each Enz to the ~.matches in interleaved format
    */
    int PPCnt=0, PairsPerLine, SinglesPerLine, NumReps, i, j, dif, MCnt=0, EOM=100, p0,
        p1, i0, i1, HTML, SameName=0;

    HTML =(int)F.HTML;

    /*    fprintf(stdout, "\n\nType of search = %c\n", PP[PPCnt].st); */

    /* need to match the ordered pairs of names in SelEnz with the corresponding RE index of the factors */
    for (PPCnt=0; PPCnt<F.Prox; PPCnt++) {
        PP[PPCnt].matches = (long *) calloc(EOM, sizeof(long));  /*get initial mem */
        if (PP[PPCnt].matches == NULL) BadMem("calloc fails - PP[].matches\n", 1);

        p0 = PP[PPCnt].REi[0];  /* clarity - the RE index of the 1st name/pattern of the prox pair */
        p1 = PP[PPCnt].REi[1];  /* ditto for the 2nd name/pattern */

        if (PP[PPCnt].REi[0] == (PP[PPCnt].REi[1])) SameName = 1; /* if comparing self to self, set indicator */

        for (i0=1; i0 <= RE[p0].E_Ncuts; i0++) { /* starts at 1 cuz [0] is the index */
            for (i1=1; i1 <= RE[p1].E_Ncuts; i1++) { /* can't these iters be started at 1 -> <= RE[p0].E_Ncuts ?  */
                if (!((SameName == 1) && (RE[p0].Sites[i0] == RE[p1].Sites[i1]))) { /* dont' do any cmps against same self site */
                    if (EOM - MCnt < 2) { /* if we're getting near the end of alloc'ed mem */
                        EOM += 100; /* bump the counter */
                        PP[PPCnt].matches = (long *) realloc(PP[PPCnt].matches, sizeof(long)*EOM);  /* and ask for more */
                        if (PP[PPCnt].matches == NULL) BadMem("calloc fails - PP[].matches\n", 1);
                    }
                    switch (PP[PPCnt].st) {

                    case '0':   /* upstr:0 gt:-1 range:-1 - default*/
                        if (abs(RE[p0].Sites[i0] - RE[p1].Sites[i1]) < PP[PPCnt].Dlo) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    case '1':   /* upstr:0 gt: 1 range:-1 */
                        if (abs(RE[p0].Sites[i0] - RE[p1].Sites[i1]) > PP[PPCnt].Dlo) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    case '2':   /* upstr:0 gt: 0 range:calc */
                        dif = abs(RE[p0].Sites[i0] - RE[p1].Sites[i1]);
                        if (dif >= PP[PPCnt].Dlo && dif <= PP[PPCnt].Dhi) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    case '3':   /* upstr:1 gt:-1 range:-1 */
                        dif = RE[p1].Sites[i1] - RE[p0].Sites[i0];
                        if ((dif > 0) && (dif < PP[PPCnt].Dlo)) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    case '4':   /* upstr:1 gt: 1 range:-1 */
                        if (RE[p1].Sites[i1] - RE[p0].Sites[i0] > PP[PPCnt].Dlo) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    case '5':   /* upstr:1 gt: 0 range:calc */
                        dif = RE[p1].Sites[i1] - RE[p0].Sites[i0];
                        /*    dif = RE[p0].Sites[i0] - RE[p1].Sites[i1];  */
                        if (dif >= PP[PPCnt].Dlo && dif <= PP[PPCnt].Dhi) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    case '6':   /* upstr:-1 gt:-1 range:-1 */
                        dif = RE[p0].Sites[i0] - RE[p1].Sites[i1];
                        if ((dif > 0) && (dif < PP[PPCnt].Dlo)) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    case '7':   /* upstr:-1 gt: 1 range:-1 */
                        if (RE[p0].Sites[i0] - RE[p1].Sites[i1] > PP[PPCnt].Dlo) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    case '8':   /* upstr:-1 gt: 0 range:calc */
                        dif = RE[p0].Sites[i0] - RE[p1].Sites[i1];
                        if (dif >= PP[PPCnt].Dlo && dif <= PP[PPCnt].Dhi) {
                            PP[PPCnt].matches[MCnt++] = RE[p0].Sites[i0];
                            PP[PPCnt].matches[MCnt++] = RE[p1].Sites[i1];
                        }
                        break;

                    default:
                        if (F.Verbose > 1) fprintf(stderr, "Illegal searchtype in Proximity()\n");
                        break;
                    }
                }
            }
        }
        PP[PPCnt].matchlen = MCnt-1;
        MCnt = 0;
    }
    PPCnt--;

    /* So, that said, 1st, let's do a CRUDE Table that will at least see if the damn code worked at all.. */
    fprintf(stdout, "\n\n\n"); /* prefix to header */
    if (HTML==1) fprintf(stdout, "<H3><A NAME=\"proxtab\">");
    fprintf (stdout,"== Proximity Matches:\n");
    if (HTML==1) fprintf(stdout, "</A></H3>\n");

    for (PPCnt=0; PPCnt < F.Prox; PPCnt++) {
        /* 8 spaces per each match + 1 for '/', rounded up to account for the margin stuff */
        PairsPerLine = (int)((F.Width+O_LMARG+O_RMARG)/17);
        SinglesPerLine = PairsPerLine*2;
        NumReps = (int)((PP[PPCnt].matchlen + 1) / SinglesPerLine) + 1;
        fprintf(stdout, "\n%s(%dE)/%s(%dE) * %d matches * found from search parameters: %s\n",
                RE[PP[PPCnt].REi[0]].E_nam, RE[PP[PPCnt].REi[0]].Err,
                RE[PP[PPCnt].REi[1]].E_nam, RE[PP[PPCnt].REi[1]].Err,
                (int)(PP[PPCnt].matchlen+1)/2, PP[PPCnt].optline);
        for (j=0; j<NumReps; j++) {
            for (i=0; (i<(SinglesPerLine) && (j*SinglesPerLine + i < PP[PPCnt].matchlen)); i+=2) {
                fprintf(stdout, "%8ld/%-8ld  ", PP[PPCnt].matches[i+j*SinglesPerLine], PP[PPCnt].matches[i+(j*SinglesPerLine)+1]);
            }
            fprintf(stdout,"\n");
        }
        //free(PP[PPCnt].matches); //can't free it here as it's needed in PrintGelLadderMap()
    }

}

/* MatchSelREs() takes RE names and (typically, if there was a hash collision, determines if it was a real
        match or just a random collision; returns a -1 if no match or the index of the matching SelEnz  */
int MatchSelREs(char *suspect, struct SE_struct SelEnz[MAX_SEL_ENZ], int NumSelREs)
{
    int i=0, match=-1;
    char *sng;
    suspect = DownCase(suspect);
    if ((sng = (char *)calloc((MAX_PAT_NAME_LEN), sizeof(char))) == NULL) {
        BadMem("MatchSelREs: No mem for sng", 1);
    }
    for (i=0; i<NumSelREs; i++) { /* go thru the whole array, NOT stopping at 1st match (for -P flag) */
        strcpy(sng,SelEnz[i].PName);
        if (strcmp(DownCase(sng), suspect) == 0) {  /* if it matches.. */
            match = i;
            /*   SelEnz[i][10] = '*';  mark it for figuring out if the enzyme was matched or not*/
            SelEnz[i].match = 1; /* mark it for figuring out if the enzyme was matched or not*/
        }
    }
    if (sng != NULL) free(sng);
    return match;
}

