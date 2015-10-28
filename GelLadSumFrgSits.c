/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: GelLadSumFrgSits.c,v 1.4 2005/02/04 04:30:17 mangalam Exp $  */

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

/* tacg.h contains all the defines, includes, function prototypes for both main() and functions */
#include "tacg.h"

/* OrfMap() prints out the ORFs found from the specification of the -O flag in the same kind of format as
seen in the -l output.  May add a option to do either ORFs or just methioneins and stops, like Strider so
user has the option of seeing both the ORFs greater than a certain size or all the possible ORFs a la Strider.
If the latter, would need --orfmap o|m (ORFs or mets/stops). ANd if the latter - will have to modify the ORF-finder code to log mets and stops to an external structure (sort of like DD) so that this fn() can use them. */
void PrintORFMap(long seq_len, int *metstop[])
{
    int   n, i, ii, j, Cnt=0, printnow=0, iNumOfDivs, iSpacesPerDiv, omc_start=0, omc_end=0, front_offset=0, ch,
                       dashptr=0, endpoint=0, HTML, n_orfs;
    long  iBasesPerDiv;
    float rBasesPerSpace=0, rBasesPerDiv, rNumOfDivs, rSpacesPerDiv, rtmp;
    char  orfbuf[MAX_LADGEL_WIDTH+4], number_line[MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG+4],
          tic_line[MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG+4];
//   char nums[10] = {'0','1','2','3','4','5','6','7','8','9'};

    //following calculates vars for printing out the map identically to laddermap
    memset(orfbuf, ' ', MAX_LADGEL_WIDTH+4); /* and initialize to ' ' */
    memset(number_line,' ',MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG+2);   /* and initialize to spaces/blanks */
    memset(tic_line,' ',MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG+2);      /* and initialize to spaces/blanks */
    /* vars, inits for the pretty print summary map */
    /* calc all the vars necessary for the ladder */
    HTML = (int) F.HTML;
    rBasesPerDiv = (float) seq_len / PREFERED_NUM_DIVS;
    iBasesPerDiv = (int) rBasesPerDiv;

    rtmp = rBasesPerDiv;
    for (Cnt=1; rtmp >= 10.0; Cnt++) rtmp = rtmp /10; /* 1st part of finding a suitable div size */
    iBasesPerDiv = (int) (rtmp);
    for (n=1; n<Cnt; n++) iBasesPerDiv *= 10; /* and once found, multiply it back up to the right size */

    rNumOfDivs = (float)seq_len / (float)iBasesPerDiv; /* real number of divs to cover whole sequence */
    iNumOfDivs = (int) rNumOfDivs; /* truncate the number to get the whole seq covered in the available space */
    rSpacesPerDiv = F.Width / rNumOfDivs;
    iSpacesPerDiv = (int) rSpacesPerDiv;
    rBasesPerSpace = (float)iBasesPerDiv / (float)iSpacesPerDiv;

    for (i=1; i<=iNumOfDivs; i++) { /* fill out number_line with the number markers */
        sprintf(number_line+((i-1)*iSpacesPerDiv + F.MaxPatNameLen), "%*ld", iSpacesPerDiv, i*iBasesPerDiv);
    }
    strcat(number_line, "\000"); /* and mark it with an end of string */
    for (i=1; i<=iNumOfDivs; i++) memcpy (&tic_line[i * iSpacesPerDiv + F.MaxPatNameLen - 1], ":", 1);
    tic_line[F.Width + F.MaxPatNameLen + 1] = '\000';  /* just to make sure that it ends OK */


    //print the (HTML) header
    if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"orfmap\">");
    fprintf(stdout, "== ORF Map Summary\n\n");
    if (HTML > -1) fprintf(stdout, "</A></H3>");

    // calc the # of orfs to limit output
    n_orfs = 0;
    i = j = 0;
    while (i< 6 && F.ORFs2Do[i] != 0) {
//      while (ORFs[i][j].orflength > 0) {
        while (ORFs[F.ORFs2Do[i]-1][j].orflength > 0) {
            n_orfs++;
            j++;
        }
        i++;
        j=0;
    }

    front_offset = F.MaxPatNameLen + 1;

    if (n_orfs > 0) {
        //and print the numbering headers
        fprintf(stdout, "%s\n", number_line);
        fprintf(stdout, "%s\n", tic_line);


        /* Since this is not a repeating operation, it's basically a 1 liner repeated up to 3 times for the
           1st 3 ORFs if requested, then for the remaining ORFs if requested, bracket by the number lines and tics, '
           so it could be represented as:
           if any of ORFS 1-3 need to be printed, print out number line and tics
           calculate position of ORFs & then print ORFs in range of 1-3
           if any of ORFs 4-5 need to be printed, print out the tic line as a separator
           calculate position of ORFs & then print ORFs in range of 4-6
           finish by printing tics and numberline again (or just tics?).
           */

        //spit out all the ORFS that need to be spat

        ii = 0;
        i = 0;
        ch = 49; // = ascii 1

        while (ii < 6 && F.ORFs2Do[ii] != 0) {
            i = F.ORFs2Do[ii] - 1;
            j=0; // j = index of the orfs in frame that are larger than the min
            while (ORFs[i][j].orflength != 0) {  //while there's still a ORF to do
                printnow = 1;
                //calc where the '>' or '<' goes and whether it needs leading/trailing '---'s
                omc_start = front_offset + (int) ((0.5 + ((float)ORFs[i][j].B_offsetBP) / rBasesPerSpace) -1); //orfmapchar & round
                omc_end   = front_offset + (int) ((0.5 + ((float)ORFs[i][j].E_offsetBP) / rBasesPerSpace) -1);
                orfbuf[F.MaxPatNameLen-1] = ch + i;
                if (i < 3) { // for frames 1-3, needs leading ---, if larger than 1 char
                    orfbuf[omc_end] = '>';  // put a '>' in the right el of buffer line
                    if (omc_start != omc_end) { // have to calc the # of '-'s needed to lead it.
                        dashptr = omc_start;
                        endpoint = omc_end ;
                        while (dashptr < endpoint) {
                            orfbuf[dashptr++] = '-'; //leading '-'s
                        }
                    }
                } else { //frames 4-6, needs lagging ---, if larger than 1 char
                    orfbuf[omc_end] = '<';
                    if (omc_start != omc_end) { // have to calc the # of '-'s needed to lag it.
                        dashptr = omc_end +1;
                        endpoint = omc_start;
                        while (dashptr < endpoint) orfbuf[dashptr++] = '-'; //lagging '-'s
                    }
                }
                j++;
            }
            if (printnow) {
                orfbuf[F.Width + F.MaxPatNameLen + 1 + O_RMARG] = '\0';
                fprintf(stdout, "%s\n", orfbuf);
                printnow = 0;
                memset(orfbuf,' ',MAX_LADGEL_WIDTH+4); /* and re-initialize to blanks */
            }
            ii++;
        }
        fprintf(stdout, "%s\n", tic_line);
        fprintf(stdout, "%s\n\n\n", number_line);

    } else {
        fprintf(stdout, "\tThere were no ORFs greater than the minimum requested\n\n");
    }


// and now the MET / STOP map

    //print the (HTML) header
    if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"metstopmap\">");
    fprintf(stdout, "== MET / STOP Map Summary\n\n");
    if (HTML > -1) fprintf(stdout, "</A></H3>");

    //and print the numbering headers
    fprintf(stdout, "%s\n", number_line);
    fprintf(stdout, "%s\n\n", tic_line);
    memset(orfbuf,' ',MAX_LADGEL_WIDTH+4); /* and re-initialize to blanks */

    ii = i = 0;
    while (ii < 6 && F.ORFs2Do[ii] != 0) {
        ch = 49;
        i = F.ORFs2Do[ii] - 1;
//		fprintf(stderr, "i=%d\n", i);
        j=1; // j = index of the mets / stops
        while (metstop[i][j]!= 0) {  //while there's still more to do
            printnow = 1;
            //calc where the '*' or '|' goes
            orfbuf[F.MaxPatNameLen-1] = ch + i;  //frame #
            if (i < 3) { // for frames 1-3,
                omc_start = front_offset + (int) (0.5 + ((float) (abs(metstop[i][j]) * 3  / rBasesPerSpace) - 1));
                if (metstop[i][j] > 0) orfbuf[omc_start] = '>';  // MET
                if (metstop[i][j] < 0) orfbuf[omc_start] = '|';  // stop
            } else { //frames 4-6, need to reverse the direction of the MET marker to < and rev the postioning of the #s as well
                omc_start = front_offset + (int) (0.5 + ((float) (((float)seq_len) - (float)abs(metstop[i][j])*3 ) / rBasesPerSpace) - 1);
                if (metstop[i][j] > 0) orfbuf[omc_start] = '<';  // MET
                if (metstop[i][j] < 0) orfbuf[omc_start] = '|';  // stop
            }
            j++;
//			fprintf(stderr, "%d ", j);
        }
//      fprintf(stderr, "out of loop\n");
        if (printnow) {
            orfbuf[F.Width + F.MaxPatNameLen + 1 + O_RMARG] = '\0';
            fprintf(stdout, "%s\n\n", orfbuf);
            printnow = 0;
            memset(orfbuf,' ',MAX_LADGEL_WIDTH+4); /* and re-initialize to blanks */
        }
        ii++;
    }
    fprintf(stdout, "%s\n", tic_line);
    fprintf(stdout, "%s\n\n", number_line);

}



/* PrintGelLadderMap does what it implies - prints a ladder map sort of
	like the GCG program, but in textmode only for now.  It uses the same width
	specification (-w) that the rest of the output does, and writes output to
	stdout rather than to the files that it did previously.  It has been
	somewhat munged to provide both a Summary map function at the same time and
	a Pseudo-gel map, as the data that it calculates is common to all these
	functions.  It has since been modified to allow the viewing of a 'log stanza'
   in the output, to expand the view of a section of the gel  */

/* MAX_LADGEL_WIDTH used as MAX parameter;  actual width will be set from flags via F.Width */

/* some changes resulted from porting it to DEC Ultrix - the system on which it was ported had
	all kinds of broken and non-standard libs.  It's sprintf() seemed to be completely
	broken, so I had to rewrite some lines to avoid using the returned value at all */

/* for adding the Prox output handling, need to modify only the ladder part and only to put out
	pairs of matches, instead of a ladder line for every RE for Proximity Matching,
   NumProtos = # of PP entries (<=10)
   seq_len   = same as for others - the len ofthe DNA seq
   *Protos   = array that holds the indices of the 'Good REs', now archaic - will use RE indices directly
   gel       = what are we trying to output
               0 = ladder map
               1 = a gel map,
               2 = ladder map without the summary map (for the --clone, frinstance)
               4 = proximity map
               5 = orfmap (maybe..)
   prox      = repetition indicator for the # of times that loops have to be traveled -
               to differentiate betw PROX and others
*/

void PrintGelLadderMap(int NumProtos, long seq_len, int *Protos, int gel)
{

    int   NumSumREs, itmp, m, n, p, i, j, k, Cnt=0, ladder_cuts[4][MAX_LADGEL_WIDTH+4], iNumOfDivs, iSpacesPerDiv,
                                             Summary_Enz[MAX_NUM_RES], ok2wr[MAX_OUTPUT_LINES], okline, min_okline, max_okline, label_len, HTML,
                                             p_width=0, spl, xtra, I_nlogs, minlog, maxlog, prox=1, SCuts=1;
    long  iBasesPerDiv, li, maxloglen;
    float rBasesPerSpace=0, rBasesPerDiv, rNumOfDivs, rSpacesPerDiv, rtmp, truval=0, F_nlogs, log_seq_len=1,
          log_frag_siz=1, fminlog;
    char  ladder_char[4][MAX_LADGEL_WIDTH+4], number_line[MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG+4], mark,
          tic_line[MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG+4],  O_txt[MAX_OUTPUT_LINES][O_RMARG+MAX_BASES_PER_LINE+4];
    char nums[10] = {'0','1','2','3','4','5','6','7','8','9'};

    if (F.GelHI == 0) {
        F.GelHI = (long)(log10((double)seq_len) + 1);
    }

    minlog = (int)F.GelLO;    /* for readability */
    maxlog = (int)F.GelHI;   /* for readability - had to jump up to the next free array el.*/
    fminlog = (float) minlog;  /* for readability */
    HTML = (int)F.HTML;     /* for readability */
    if (F.Prox != 0) prox = 4;

    /* set/clear all the common arrays */
    for (i=0; i<4; i++) {
        memset(ladder_cuts[i], 0, sizeof(int)*MAX_LADGEL_WIDTH+4); /* and initialize to 0 */
        memset(ladder_char[i], ' ', MAX_LADGEL_WIDTH+4); /* and initialize to 0 */
    }
    memset(number_line,' ',MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG+2);   /* and initialize to spaces/blanks */
    memset(tic_line,' ',MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG+2);      /* and initialize to spaces/blanks */

    if ((gel == 1) && (maxlog != 0)) { /* if maxlog has been set to something via -g Lo,Hi... */
        maxloglen = (long)pow(10, (double)maxlog);
        /* below is the *simple* hack that allows you to slice out the required 'log stanza'
           and allows you to set the max LARGER than the sequence length */
        seq_len = maxloglen;
    }

    /* Separate the Gel stuff from the Ladder stuff */
    if (gel==1)  { /* for gel map only */
        if (seq_len < 100) {
            fprintf(stderr, "Sequence is too short (%ld bases) for a decent Gel map - bye..\n", seq_len);
            exit(1);
        }
        /* have to calc #s for labeling the axis */
        log_seq_len =  (float)log10((double)seq_len);
        F_nlogs = log_seq_len - fminlog;  /* the difference betw start and seq_len */
        I_nlogs =  (int) F_nlogs;
        spl = (int) (F.Width / F_nlogs); /* spl = spaces per log */

        /* setting up the tics and numbers */
        li = 1;
        /*    next line allows ANY lower cutoff, not just 10 or 1000  */
        for (i=0; i<minlog; i++) li = li * 10; /* generate the original number for -g flag */
        sprintf(number_line, "%*ld", (F.MaxPatNameLen + 1), li);
        for (i=0; i<=I_nlogs; i++) {
            li *= 10;
            if (li <= seq_len) sprintf(number_line+(i*spl+(F.MaxPatNameLen + 1)), "%*ld", spl, li);
            for (j=1; (j<10 && truval<log_seq_len); j++) {
                truval = (float)i+minlog + (float)log10((double)j);
                Cnt = (i*spl + F.MaxPatNameLen - 1)+(int)(log10((double)(j)) * (float)spl + 1.5);
                if (Cnt < ((int)F.Width)+F.MaxPatNameLen +2) tic_line[Cnt] = '.';   /* sanity check... */
            }
        }
        tic_line[Cnt+1] = '\0';
        strcat(number_line, "\000");/* to terminate correctly */
        /* And print a minimally descriptive header */
        fprintf(stdout, "\n\n");
        if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"gelmap\">");
        fprintf(stdout, "== Pseudo-Gel Map of Digestions:  ");
        if (HTML > -1) fprintf(stdout, "</A></H3>");
        /* also need HTML for the prox stuff */

    } else {  /* for ladder map only */
        if (SUMMARY_CUTS == 0) { /* sets SCuts to increase over 10,000 bases */
            if (seq_len >10000) {
                SCuts = (int)log10((double)seq_len) - 2;
            }
        }

        if (seq_len < 10) { /* seqs less than 10 won't work here because of the way the tics are calc'ed */
            fprintf(stderr, "Sequence is too short (%ld bases) for a decent Ladder map - bye..\n", seq_len);
            exit(1);
        }
        /* vars, inits for the pretty print summary map */
        p_width = F.MaxPatNameLen - O_LMARG + F.Width+O_RMARG;
        /* calc all the vars necessary for the ladder */
        rBasesPerDiv = (float) seq_len / PREFERED_NUM_DIVS;
        iBasesPerDiv = (int) rBasesPerDiv;

        rtmp = rBasesPerDiv;
        for (Cnt=1; rtmp >= 10.0; Cnt++) rtmp = rtmp /10; /* 1st part of finding a suitable div size */
        iBasesPerDiv = (int) (rtmp);
        for (n=1; n<Cnt; n++) iBasesPerDiv *= 10; /* and once found, multiply it back up to the right size */

        rNumOfDivs = (float)seq_len / (float)iBasesPerDiv; /* real number of divs to cover whole sequence */
        iNumOfDivs = (int) rNumOfDivs; /* truncate the number to get the whole seq covered in the available space */

        rSpacesPerDiv = F.Width / rNumOfDivs;
        iSpacesPerDiv = (int) rSpacesPerDiv;

        rBasesPerSpace = (float)iBasesPerDiv / (float)iSpacesPerDiv;
        for (j=0; j<prox; j++) {
            memset(ladder_cuts[j],0,sizeof (int)* MAX_LADGEL_WIDTH); /* and re-initialize to 0 */
            memset(ladder_char[j],'-',MAX_LADGEL_WIDTH+4); /* and re-initialize to dashes */
        }

        for (i=1; i<=iNumOfDivs; i++) { /* fill out number_line with the number markers */
            sprintf(number_line+((i-1)*iSpacesPerDiv + F.MaxPatNameLen), "%*ld", iSpacesPerDiv, i*iBasesPerDiv);
        }
        strcat(number_line, "\000"); /* and mark it with an end of string */
        for (i=1; i<=iNumOfDivs; i++) memcpy (&tic_line[i * iSpacesPerDiv + F.MaxPatNameLen - 1], ":", 1);
//      tic_line[F.Width + F.MaxPatNameLen + 1] = '\000';  /* just to make sure that it ends OK */
        tic_line[i * iSpacesPerDiv + F.MaxPatNameLen] = '\000';  /* just to make sure that it ends OK */
//		fprintf(stderr, "memcpy length = %d, end posn=%d\n", (i * iSpacesPerDiv + F.MaxPatNameLen - 1), (F.Width + F.MaxPatNameLen + 1));

        /* Print a minimally descriptive header */
        fprintf(stdout, "\n\n");
        if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"ladmap\">");
        fprintf(stdout, "== Ladder Map of Enzyme Sites:  ");
        if (HTML > -1) fprintf(stdout, "</A></H3>");
    }
    /* need to print the number_line once every 50 iterations and the tics every TIC_REPEAT
       iterations, so yet another level of iffiness is required */
    /* COMMON CODE */
    /* if minimum cuts (-m) flag has been set */
    if ( F.Min != 1) fprintf(stdout, "  *Minimum* Hits: %d",  F.Min);
    /* if maximum (-M) flag has been set */
    if ( F.Max != 32000) fprintf(stdout, "  *Maximum* Hits: %d", F.Max );
    if (gel==1) fprintf(stdout, "\n%*s\n", (int)(F.MaxPatNameLen + O_RMARG+F.Width),"# frags");
    else  fprintf(stdout, "\n%*s\n", (int)(F.MaxPatNameLen + O_RMARG+F.Width),"# hits");

    Cnt = 0;
    NumSumREs=0;
    for (m = 0; m < NumProtos; m++) { /* for Prox, NumProtos == PPCnt */
        if (RE[Protos[m]].E_Ncuts <= SCuts) {
            Summary_Enz[NumSumREs++]= Protos[m]; /* gathering REs for Summary map */
        }
        if ((RE[Protos[m]].E_mag  >= (int)F.Mag) &&    /* E_mag must be >= to -n flag */
                (RE[Protos[m]].E_Ncuts >= (int)F.Min) &&    /* Ncuts >= -m, */
                (RE[Protos[m]].E_Ncuts <= (int)F.Max)) {    /* Ncuts <= -M */

            mark = '|';   /* set mark here, not below; overwrite if needed (1 source of garbled format before) */
            if ((F.Overhang == 1)      || /* if want all, only 1st has to be true to short circuit */
                    ((RE[Protos[m]].E_olap > 0) && (F.Overhang == 5)) || /* if olap < 0 && overhang = 5' */
                    ((RE[Protos[m]].E_olap < 0) && (F.Overhang == 3)) || /* if olap > 0 && overhang = 3' */
                    ((RE[Protos[m]].E_olap ==0) && (F.Overhang == 0))) { /* if olap = 0 && overhang = blunt */

                if (RE[Protos[m]].E_olap > 0) mark = '\\';  /* set up the markers for the output */
                else if (RE[Protos[m]].E_olap < 0) mark = '/';
                /* if need repeated headers or tics, print them.... */
                if ((Cnt % NUMBERLINE_REPEAT) == 0) fprintf(stdout, "\n%s", number_line);
                if ((Cnt++ % TIC_REPEAT) == 0)      fprintf(stdout, "\n%s", tic_line);
                /* and now the real data... */
                if (gel != 4) fprintf (stdout, "\n%*s ", F.MaxPatNameLen, RE[Protos[m]].E_nam); /* print the name of the filtered RE */

                /* calculate where the cut marks go... */
                if (gel == 1) xtra = RE[Protos[m]].E_Nfrags;
                else xtra = RE[Protos[m]].E_Ncuts;

                if (gel==1)  { /* calc fragment migration via simple log approximation */
                    memset(ladder_char[0],' ',MAX_LADGEL_WIDTH); /* and re-initialize to spaces */
                    mark = '|'; /* reset mark for gel map - no diffs for 5', 3', blunts... */
                    if (F.Topo == 1) xtra++; /* for 'xtra' fragment in linear digest */
                    for (j=0; j<xtra; j++) {

                        /* this is where the MAX cutoff calc has to go for doing a window into the Gel
                        	(commented for now..)  */
                        if (RE[Protos[m]].Frags[j]>0) log_frag_siz = (float)log10((double)RE[Protos[m]].Frags[j]);
                        /* if (log_frag_siz >= fmaxlog) ladder_cuts[0][MAX]++; */   /* count the # of frags MORE than the MAX shown */
                        if (log_frag_siz <= fminlog) {
                            ladder_cuts[0][0]++;  /* count the # of frags LESS than the MIN shown */
                        } else {
                            n = (int) (0.5 + (float)F.Width * (log_frag_siz - fminlog)/(log_seq_len - fminlog));  /* debug line */
                            ladder_cuts[0][n]++;
                        }
                    }
                } else { /* it's for the regular or Proxmity ladder  */
                    if (gel == 4) { /* gel = 4 indicates a Prox Ladder type */
                        /* diff is that Prox handling requires that 2 lines be composed at the same time, leading to these
                           silly double 'for' loops ... */
                        /* Next stanza loads the MATCHES ladder_cuts and then ladder_char with an inner for loop that uses some
                           unpleasant offsets to hit them correctly.  The output for both 'cuts' and 'char' is supposed to be:
                           [0] - regular cuts for NameA done 2nd
                           [1] - MATCHES only for NameA done 1st
                           [2] - MATCHES only for NameB done 1st
                           [3] - regular cuts for NameB done 2nd
                        reminder - m=NProtos (=PP counter in the Proximity case) */
                        for (j=0; j<=PP[m].matchlen; j+=2) { /* loads both [1] and [2] above */
                            for (i=0; i<2; i++) {
                                itmp = ((int) (0.5 + ((float) PP[m].matches[j+i]  )/rBasesPerSpace));
                                if (itmp < 0) itmp = 0;  /* correcting for loss of cuts too close to the beginning */
                                ladder_cuts[i+1][itmp]++; /* loads both [1] and [2] */
                            }
                        }
                        for (i=0; i<2; i++) { /* loads both [0] & [3] - this stanza ~ident to below cept for indexing  */
                            for (j=1; j<=RE[PP[m].REi[i]].E_Ncuts; j++) {
                                itmp = ((int) (0.5 + ((float) RE[PP[m].REi[i]].Sites[j])/rBasesPerSpace));
                                if (itmp < 0) itmp = 0;  /* correcting for loss of cuts too close to the beginning */
                                ladder_cuts[i*3][itmp]++; /* loads both [0] and [3] */
                            }
                        }
                    } else { /* for ONLY a regular ladder - NO proximity stuff involved in this stanza  */
                        for (j=1; j<=RE[Protos[m]].E_Ncuts; j++) {
                            itmp = (int) ((0.5 + ((float) RE[Protos[m]].Sites[j] )/rBasesPerSpace) - 1);
                            if (itmp < 0) itmp = 0;  /* correcting for loss of cuts too close to the beginning */
                            ladder_cuts[0][itmp]++; /* loads only the [0] element */
                        }
                    }
                }
                /* now all the els of ladder_cuts are incr'd for the whole SF_data data; transfer them to ladder_char */
                for (j=0; j<prox; j++) { /* want this to go 0 -> 3 if prox(=4), only 0 if not (=1) */
                    for (i=0; i<=(int)F.Width; i++) { /* from 0 to print width (-w), allowing for 1 char of jitter at end */
                        if (ladder_cuts[j][i] != 0) {
                            if (gel == 4) mark = '|';
                            if (ladder_cuts[j][i] == 1) ladder_char[j][i] = mark;
                            else {
                                if (ladder_cuts[j][i] <10) {
                                    memcpy(&ladder_char[j][i], &nums[ladder_cuts[j][i]],1);
                                } else {
                                    ladder_char[j][i] = '*';
                                }
                            }
                        }
                    }
                }

                /* before printing ladder_char, also have to clear the '-'s that the sequence didn't get to */
                if (gel != 1) {
                    for (j=0; j<prox; j++) {
                        for (i=(int)(0.5 + seq_len/rBasesPerSpace); i<=(int)F.Width; i++) ladder_char[j][i] = ' ';
                    }
                }
                if (gel == 4) {
                    for (j=0; j<prox; j++) {
                        ladder_char[j][F.Width+1] = '\000';
                    }
                    fprintf(stdout,"\n\n%9s %*s %4d (All)", RE[PP[m].REi[0]].E_nam, (int)F.Width, ladder_char[0], RE[PP[m].REi[0]].E_Ncuts); /* print the prox names preceding the ladder */
                    fprintf(stdout,"\n%9s %*s %4d (Matched)", RE[PP[m].REi[0]].E_nam, (int)F.Width, ladder_char[1], (int)((PP[m].matchlen+1)/2)); /* print the prox names preceding the ladder */
                    fprintf(stdout,"\n%9s %*s %4d (Matched)", RE[PP[m].REi[1]].E_nam, (int)F.Width, ladder_char[2], (int)((PP[m].matchlen+1)/2)); /* print the prox names preceding the ladder */
                    fprintf(stdout,"\n%9s %*s %4d (All)\n", RE[PP[m].REi[1]].E_nam, (int)F.Width, ladder_char[3], RE[PP[m].REi[1]].E_Ncuts); /* print the prox names preceding the ladder */
                    fprintf(stdout, "\n%s", tic_line);
                    for (j=0; j<prox; j++) {
                        memset(ladder_cuts[j],0,sizeof (int)* MAX_LADGEL_WIDTH);
                        memset(ladder_char[j],'-',MAX_LADGEL_WIDTH+4);
                    }
                    free(PP[m].matches); //hjm
                } else {
                    /* but if either gel or regular ladder, don't print the name cuz it was already printed above (~178) */
                    ladder_char[0][F.Width+1] = '\000';
                    fprintf(stdout, "%*s %4d", (int)F.Width, ladder_char[0], xtra);  /* and then print out the corresp ladder */
                    memset(ladder_cuts[0],0,sizeof (int)* MAX_LADGEL_WIDTH+4); /* and re-initialize to 0 */
                    memset(ladder_char[0],'-',MAX_LADGEL_WIDTH+4); /* and re-initialize to dashes */

                }
            }
        }
    }
    /* at this point, have indices for all REs that cut <= SUMMARY_CUTS in Summary_Enz[], so can now append the summary..
       ladder_cuts, ladder_char are both reset, so can re-use them without further cost */
    /* Set up the output..*/
    if (gel == 0) {  /* conditional for printing the summary map */
        if (NumSumREs > 0) {
            memset(ladder_char[0],'-',MAX_LADGEL_WIDTH); /* and re-initialize to dashes */
            okline = min_okline = O_SEQ_LINE-2; /* set the min (highest on page) line of buffer to print out for Summary map */
            max_okline = O_SEQ_LINE + 3; /* how many more lines beyond dashes are required for output */
            memset(O_txt,' ', (MAX_OUTPUT_LINES * (O_RMARG+MAX_BASES_PER_LINE))); /* and initialize to spaces/blanks */
            for (i=0; i<MAX_OUTPUT_LINES; i++) ok2wr[i] = 0; /* set/reset all of ok2wr to O_LMARG for the new block*/
            /* before printing ladder_char, also have to clear the '-'s that the sequence didn't get to */
            for (i=(int)(0.5 + seq_len/rBasesPerSpace); i < (int)F.Width; i++) ladder_char[0][i] = ' ';

            fprintf(stdout, "\n\n");
            if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"summarycuts\">");
            fprintf(stdout, "== Summary of Enzymes that hit ** %d ** times or less:\n\n", SCuts);
            if (HTML > -1) fprintf(stdout, "</A></H3>");

            memcpy(&O_txt[O_SEQ_LINE][F.MaxPatNameLen+1], &ladder_char[0], F.Width); /* copy the dashes into place */
//			fprintf(stderr, "F.MaxPatNameLen=%d, O_LMARG=%d, start=%d, [%s]\n", F.MaxPatNameLen, O_LMARG, (F.MaxPatNameLen - O_LMARG), tic_line);
            memcpy(&O_txt[O_SEQ_LINE+1][0], &tic_line[F.MaxPatNameLen - O_LMARG], F.Width+30); /* copy the marker line into place */
            memcpy(&O_txt[O_SEQ_LINE+2][0], &number_line[F.MaxPatNameLen - O_LMARG], F.Width+30); /* copy the numbers into place */

            /* calculate where the cut marks go... */
            for (i=0; i<NumSumREs; i++) {    /* while there's still more REs to go thru */
                strcpy(ladder_char[0], RE[Summary_Enz[i]].E_nam);
                strcat(ladder_char[0], "@"); /* common for all cuts with this enz - next place to copy is at E_nam_l+1 */
                for (j=1; j<=RE[Summary_Enz[i]].E_Ncuts; j++) {      /* while there's still more cuts to go thru */
                    /* finish composing the marker name, using number_line, since we don't need it anymore  */
                    memset(number_line,' ',MAX_LADGEL_WIDTH+F.MaxPatNameLen - O_LMARG);  /* initialize to spaces/blanks */

                    sprintf(number_line, "%ld", RE[Summary_Enz[i]].Sites[j]); /* finish composing name by adding cut position  */

                    Cnt = strlen(number_line); /* stupid Ultrix segfaults on itmp = sprintf(); this is the alternative*/
                    memcpy(&ladder_char[0][RE[Summary_Enz[i]].E_nam_l+1], &number_line, Cnt+1); /* has no str termination tho */
                    label_len = strlen(ladder_char[0]);

                    itmp = ((int) ( F.MaxPatNameLen + 0.5 + ((float)RE[Summary_Enz[i]].Sites[j])/rBasesPerSpace)) ; /* where the cut maps to */
                    if (itmp < 0) itmp = 0; /* corrects for enz that map vclose to begin that would otherwise map to -1 */
                    O_txt[O_SEQ_LINE-1][itmp] = '|';   /* write the 'cut' indicator on the line above dashes */
                    /* locate the position for writing the name */
                    while ((itmp < ok2wr[okline]) && (okline != 0)) {
                        okline--; /* then go up to the lowest 'free' line */
                    }

                    if (okline == 0) {
                        fprintf(stdout, "Enz Names stack too high for the 'Summary' buffer space - increase '-w' value\n"
                                "to spread them out more.\n");
                        exit(1);
                    }
                    if (okline != O_SEQ_LINE-2)  { /* but if can't place the name on the lowest line.. */
                        p = O_SEQ_LINE-1;
                        k = 0;   /* then check if there's any space lower between previous names */
                        /* following 'while' must be ordered like it is - p will be decr only if k == 0 */
                        while (k == 0 && --p != okline) { /* leaving 1 space before and after the name */
                            if (O_txt[p][itmp-1] == ' '         &&     /* this can be replaced with something */
                                    O_txt[p][itmp+label_len] == ' ' &&     /* more efficient later */
                                    O_txt[p][itmp+4] == ' '         &&     /* must check for a space here..- this one might not be nec. */
                                    O_txt[p][itmp+5] == ' ' )  k = 1;      /* as well as here */
                        }
                        okline = p; /* p = vertical counter ~= okline, k = 1 when we find enuf space */
                    }
                    /* and memcpy the RE label from the struct to the output array */
                    memcpy (&O_txt[okline][itmp],&ladder_char[0], label_len); /* don't want the NULL term if it gets to here!! */
                    /* and incr the ok2wr posn for that line by the length of the RE name + 1 */
                    /* 'ok2wr[okline]' below is modified only if we didn't find any lower spaces */
                    if (itmp + label_len + 1 > ok2wr[okline])  ok2wr[okline] = itmp + label_len + 1;
                    if (okline < min_okline) min_okline = okline;
                }     /* for (j=1;j<=RE[Summary_Enz[i]].E_Ncuts; j++) .. */

            }   /* for (i=0; i<NumSumREs; i++) ... */
            /* Now print out the block we've composed, from min_okline to max_okline */

            for (i=min_okline; i<max_okline; i++) {
                fprintf(stdout, "%.*s\n",p_width,O_txt[i]);
            }
            /* And clean out the lines written into O_txt as well */
            memset(O_txt,' ',((O_RMARG + MAX_BASES_PER_LINE) * MAX_OUTPUT_LINES)); /* set all of O_txt to blanks */
        } else { /* there's no need to print a summary of no output, so just print an error-like message */
            if (F.Verbose >= 1) {
                fprintf(stderr, "\nNo patterns matched %d times or less, so the summary map was skipped.\n", SCuts);
            }
            fprintf(stdout, "\n\n");
            if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"summarycuts\">");
            fprintf(stdout, "== No patterns matched ** %d ** times or less, so Summary Map skipped.\n\n", SCuts);
            if (HTML > -1) fprintf(stdout, "</A></H3>");
        }
    } /*if (!gel) ... */
    fprintf(stdout,"\n");     /* make sure the output ends with a newline */
} /*  end of  PrintLadderMap() */



/* PrintSitesFrags prints out the fragments either sorted or unsorted, depending on whether the array has been sorted -
   it doesn't care, it just prints whatever's in the array nicely */
/* sites       indicator for whether the fn() should do site data (1) or fragment data (0)
*/
void PrintSitesFrags(int NumProtos, int reps, int *Protos, int sites, long seq_len)
{
    /* sites signals whether we want to print sites (=1), unsorted fragments(=0) or sorted fragments(=-1) */
    int b, i, j, k, l, m, xtra, HTML, namei;
    char *names[] =  {"== Fragment Sizes", "== Match Sites", "Fragment", "Match"};
    reps = reps-2; /* needs to be reduced by 2 for this fn() */
    HTML = (int)F.HTML; /* for readability */

    //fprintf(stdout, "\n\n"); /* prefix to PrintFrag() title */
    if (HTML > -1) { /* need to print the name tags */
        if (sites==1) {
            fprintf(stdout, "<H3><A NAME=\"sitetab\">");
        } else if (sites == 0) {
            fprintf(stdout, "<H3><A NAME=\"fragunsort\">");    /* unsorted frags */
        } else {
            fprintf(stdout, "<H3><A NAME=\"fragsort\">");    /* sites = -1 = sorted frags */
        }
    }
    /* then print the common headers */
    if (sites >= 0) namei = sites; /* sites = 0 || sites = 1 OK, but ... */
    else namei = 0; /* sites = -1 -> 0  */
    //fprintf(stdout, "%s by Enzyme or Pattern", names[namei]);
    if (sites == 0) {
        fprintf(stdout, " (UNSORTED)\n");
    } else if (sites == -1) {
        fprintf(stdout, " (SORTED)\n");
    }
    if (HTML > -1) fprintf(stdout, "</A></H3>"); /*  and the end of the 'name + header' for HTML  */

    if (F.Min != 1) {  /* if minimum cuts (-m) flag has been set */
        fprintf(stdout, "   *Minimum* Hits: %d",  F.Min);
    }
    if (F.Max != 32000) {  /* if maximum (-M) flag has been set */
        fprintf(stdout, "   *Maximum* Hits: %d", F.Max );
    }
    //fprintf(stdout, "Pattern: ");
    for (m = 0; m < NumProtos; m++) {
//		fprintf(stderr, "for %s,  hits = %d \n", RE[Protos[m]].E_nam, RE[Protos[m]].E_Ncuts);
        if (RE_Filters_OK(Protos[m])) {
//      		fprintf(stderr, " %s passed\n", RE[Protos[m]].E_nam);
            /* this 'xtra' stuff handles the topo differences (linear has an extra fragment), */
            xtra = RE[Protos[m]].E_Ncuts;      /* assign xtra to # cuts */
            if (sites <= 0 && F.Topo == 1) xtra++; /* if frags && topo is linear, # frags will be Ncuts+1 */

            if (F.Matrix != 0) { /* its a matrix match, so change the header slightly */
                fprintf (stdout, "\n%-10s  TFC = %s;  Fwd> = %s; Rev< = %s;  %d %s(s) \n",
                         RE[Protos[m]].E_nam, RE[Protos[m]].E_raw_sit, RE[Protos[m]].E_wsit[0], RE[Protos[m]].E_wsit[1],
                         xtra, names[namei+2]);
            } else {
                fprintf (stdout, "\n%-10s  %s (%d Err) - %d %s(s) found\n",
                         RE[Protos[m]].E_nam, RE[Protos[m]].E_raw_sit, RE[Protos[m]].Err, xtra, names[namei+2]);
            }

            l = (int)((xtra/reps) + 1);
            i=0;
            if (sites==1) {
                i = 1;    /* for the difference between how the arrays are organized */
                xtra++;
            }
            for (k=0; k<l; k++) {
                if (sites == 1) { /* we want a table of Sites */
                    /* if want to add '-' to nonpals that are in the oppo orientation, have to check
                       RE[x].Orient[i] for whether each site is on + or - strand.  It's SET if it's on the
                       - strand  */
                    /* use the one below instead of this:  fprintf(stdout, "%8ld", RE[Protos[m]].Sites[i]); */
                    for (j=0; j<(reps) && i<xtra; j++,i++) {
                        b = 1;
                        if (F.Sites == 2) {
                            if (BitArray(RE[Protos[m]].Orient,i-1,2)) { /* if the corresp bit is set (=1)  */
                                b = -1;
                            }
                        }
                        fprintf(stdout, "  %ld", RE[Protos[m]].Sites[i] * b);
                    }
                } else { /* else we want a table of Frags */
                    for (j=0; j<(reps) && i<xtra; j++,i++) {
                        fprintf(stdout, "  %ld",RE[Protos[m]].Frags[i]);
                    }
                }
                if (j == reps && i == xtra) k++; /* stop the extra \n at the end of an exact rectangle */
                fprintf(stdout, "\n");
            }
        }
        /* } */
    }
    /* fprintf(stdout, "\n\n");  */
}
