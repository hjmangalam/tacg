/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: SeqFuncs.c,v 1.5 2005/01/26 20:16:25 mangalam Exp $  */

/* The use of this software (except that by Harald T. Alvestrand, which is described in 'udping.c')
   and that by James Knight, which is described in 'seqio.c' is bound by the notice that appears
   in the file 'tacg.h' which should accompany this file.  In the event that 'tacg.h' is not bundled
   with this file, please contact the author.
*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <termios.h>

#include "seqio.h"
#include "tacg.h"

/* GrafMap is the 1st take on the long-delayed graphical plasmid/linear maps function */
/* on calling, runs thru RE and makes a 1st pass guestimate of how to map all the hits. If there
	are too many hits, it should truncate at a level which does not leave a pattern partially
   mapped.   As can be determined by examining the Postscript generated, I'm a PS newbie and the PS
   could be made significantly more compact by encoding the tic routines in PS for example and letting
   the PS engine do all the calculations. It's also full of debug info, which will remain until the
   code solidifies more.  The PDF output tho is quite robust.
*/

void GrafMap(long *Degen_Log, int *Protos, char *Description, char datestamp[80], float FrACGT[4])
{

    int i, j, n, bi, ti, sli, HPc, REi, FCnt, ORF_line_width, f, ORFNum,
        need_prolog=0, NLabels, RWCincr, previ, start, end, run,center, abegin,
        HPEnd = NUM_PLASMID_LABELS+30, debug=0, dstart=0, dend=15;
    /* set the debug vars */

    long HPi=0, HPMax=0, RWC, ldiff, BasesPerLabel, seq_len = F.Seq_Len, temp,
         HP[3][HPEnd], mdbl; /*  min dist between labels */
    char  *MapFileName,
          *o, star, *arrow,
          Label[20];
    float ORF_density, x_left, y_top, y_bot, yoffset;
    /* some constants, measurements in inches, arcs in degrees, */
    double Radius 		= 3.0,
               CenterX 		= 4.25,
                     CenterY 		= 4.5,

                           AnglePerLabel,

                           arc_end, hairline_radius,
                           PtCenterX, PtCenterY, PtRadius;

    FILE *PS;     /* output postscript file */

    struct NextHit_struct NH;
    struct stat *Fbuf;

    /* mangle the file name to be compatible with a cgi call if it was so called */
    if ((MapFileName = calloc(128, sizeof(char))) == NULL) {
        BadMem("SeqFuncs.c:MapFileName", 1);
    }
    if (F.HTML != -1) { /* then we need the /tmp-named files */
        sprintf(MapFileName, "%s/tacg_Map.ps", F.tmppath);
    } else {
        MapFileName = "tacg_Map.ps\0"; /* name by itself to write into the $PWD */
    }
    PtCenterX = CenterX * 72;
    PtCenterY = CenterY * 72;
    PtRadius  = Radius  * 72;
    x_left = y_bot = 25; /* ~1/3" up from the margin to start */
    y_top = 756; /* ~10.5" top margin */

    NLabels = NUM_PLASMID_LABELS; /*  how many Label slots there are */
    mdbl = seq_len / (long)NLabels;
    AnglePerLabel = 360 / NUM_PLASMID_LABELS;
    BasesPerLabel = (long)((long)AnglePerLabel * seq_len / 360);

    for (i=0; i<HPEnd; i++) {
        for (j=0; j<3; j++) {
            HP[j][i] = 0;
        }
    }

    if ((Fbuf = calloc(1, sizeof(*Fbuf))) == NULL) {
        BadMem("RecentFuncs:Fbuf",1);
    }
    if (lstat(MapFileName, Fbuf) == -1) {
        need_prolog = 1;
    }

    if ((PS=fopen(MapFileName,"a")) == NULL) {  /* or the standard/optional one from SetFlags() */
        fprintf(stderr,"Cannot open the output postscript file \"%s\" for writing!!\n", MapFileName);
        exit(1); /* print an error and die gracefully */
    }

    if (need_prolog == 1) {
        /* use the Illustrator PS prolog as it has a number of useful defs, so copy PROLOG into PS */
        fprintf(PS, "%s", tacg_PS_Prolog());
    }

    /* ------------------------  write the title, date, etc --------------------------------------- */
    /* write a header at top left */
    yoffset = y_top;
    fprintf(PS, "\n\n/Helvetica-Bold findfont 12 scalefont setfont\n");
    fprintf(PS, "newpath %f %f moveto\n (Description: %s) show      \n",
            x_left, yoffset, Description);
    yoffset -= 12;
    fprintf(PS, "%f %f moveto\n (Date: %s) show\n", x_left, yoffset, datestamp );

    fprintf(PS, "\n\n/Helvetica findfont 10 scalefont setfont\n");
    yoffset -= 10;
    fprintf(PS, "newpath %f %f moveto\n (Map by tacg v. %s (http://tacg.sf.net)) show\n",
            x_left, yoffset, TACG_VERSION);
    fprintf(PS, "\n\n/Helvetica findfont 10 scalefont setfont\n");
    yoffset -= 10;
    fprintf(PS, "newpath %f %f moveto\n (Suggestions to Harry Mangalam (hjm@tacgi.com)) show\n",
            x_left, yoffset);

    /* print a little legend at bottom left */
    yoffset = 30;
    fprintf(PS, "\n/Helvetica findfont 8 scalefont setfont\n");
    fprintf(PS, "newpath %f %f moveto\n (Label = Name[base#]overhang (*), where: ) show\n", x_left, (y_bot+yoffset));
    yoffset -= 8;
    fprintf(PS, " %f %f moveto\n (Name = Name of the RE or pattern,  ) show\n", x_left, (y_bot+yoffset));
    yoffset -= 8;
    fprintf(PS, " %f %f moveto\n ([base#] = cut site of the RE or ~center of the pattern) show\n",  x_left, (y_bot+yoffset));
    yoffset -= 8;
    fprintf(PS, " %f %f moveto\n (overhang: \\\\ = 5', / = 3', | = blunt. '*' = 1 hit only.) show\n",  x_left, (y_bot+yoffset));

    /* assume an 8x10.5 US letter page; will probably have to correct to various media later */

    /* set the line width and then draw a circle at the center of the page */
    fprintf(PS, "1 setlinewidth\nnewpath %f %f %f 0 360 arc stroke\n", PtCenterX, PtCenterY, PtRadius);
    /* outline circle to set off degens */
    fprintf(PS, "0.3 setlinewidth\nnewpath %f %f %f 0 360 arc stroke\n", PtCenterX, PtCenterY, PtRadius+3);

    /* write a small + at the center for a registration mark */
    fprintf(PS, " .3 setlinewidth %f %f moveto\n 10 0 rlineto stroke\n %f %f moveto\n 0 10 rlineto stroke\n",
            (PtCenterX - 5), PtCenterY, PtCenterX, (PtCenterY - 5));


    /* --------------------------   and the sequence stats -----------------------------------  */

    fprintf(PS, "\n/Helvetica-Bold findfont 10 scalefont setfont\n");
    fprintf(PS, "newpath %f %f moveto\n", (PtCenterX - 30), PtCenterY+30);
    fprintf(PS, " (%ld bases)  show\n", seq_len) ;


    fprintf(PS, "\n/Helvetica findfont 6 scalefont setfont\n");
    fprintf(PS, "newpath %f %f moveto\n", (PtCenterX - 100), PtCenterY+10);
    fprintf(PS, " (%ld A(%.2f %%)  %ld C(%.2f %%)  %ld G(%.2f %%)  %ld T(%.2f %%))  show\n",
            F.NBases[0], FrACGT[0]*100, F.NBases[1], FrACGT[1]*100, F.NBases[2],
            FrACGT[2]*100, F.NBases[3], FrACGT[3]*100);


    if (seq_len > (temp = F.NBases[0] + F.NBases[1] +F.NBases[2] + F.NBases[3])) {
        fprintf(PS, "newpath %f %f moveto\n", (PtCenterX - 90), PtCenterY-15);
        fprintf (PS, "\n(NB: sequence length > A+C+G+T due to %ld IUPAC degeneracies.) show\n", seq_len-temp);
        fprintf(PS, "newpath %f %f moveto\n", (PtCenterX - 80), PtCenterY-22);
        fprintf (PS, "(# of:  N:%ld  Y:%ld  R:%ld  W:%ld  S:%ld  K:%ld  M:%ld  B:%ld  D:%ld  H:%ld  V:%ld) show \n",
                 F.NBases[14], F.NBases[5], F.NBases[6], F.NBases[7], F.NBases[8], F.NBases[9], F.NBases[10],
                 F.NBases[10], F.NBases[11], F.NBases[12], F.NBases[13]);
    }

    /* write out the alternating light/dark 1KB bands ?  not yet.. */

    /* start counting the number of sites and enter them into a local array */
    HPc = 0;
    j=2;
    RWCincr = seq_len / NLabels;
    while (HPc < NLabels && j < (*DD).Dat[0] && ((*DD).Dat[j] < seq_len)) {
        RWC = (*DD).Dat[j];
        REi = (int)labs((*DD).Dat[j+1]);
        if (RE_Filters_OK(REi)) { /* if the RE is allowed by the filter */
            HPi = (int)((float)(*DD).Dat[j] / (float)seq_len * NLabels); /* HPi points to the relative position of the hit */
            /* bc circular sequences are padded, possible to find hits
            	before and beyond actual sequence, so have to check for over/underruns */
            if (HP[0][HPi] == 0) {     /* if there isn't anything already there */
                HP[0][HPi] = RWC; /* the RWC of the hit */
                HP[1][HPi] = REi;          /* the RE index */
                ldiff = HP[0][HPi-1] - HP[0][HPi];
                if (HPc > 0 && ldiff > 0 && ldiff < mdbl) {
                    HP[2][HPi] += RWCincr;
                }
            } else { /* have to find a nearby open slot in the appro position */
                /* since this is taken from DD.Dat, the sites WILL BE in order, so should just
                   have to step to the NEXT open slot */
                /* this routine does not guarantee that the lines won't cross, altho it should work in most cases
                   in order to guarantee this, would have to make it into a linked list structure - possible
                   but lets wait to see how it behaves in real life. */
                /*                fprintf(stderr,"\nCOLLISION @ %d!\n", HPi); */
                while (HP[0][HPi] != 0 && HPi < HPEnd) { /* incr until it finds an empty slot */
                    HPi++;
                }

                if (HPi >= HPEnd) {
                    /* then we ran into the end, but there are probably empty slots left before this
                       so let's see if we can find them */
                    sli=HPi;
                    while (HP[0][HPi] != 0) {
                        sli--;
                    }
                    if (sli <= 0) {   /* we really did use them all up, so die */
                        fprintf(stderr, "\nLooks like we ran out of Label Slots in GrafMap:HP\n");
//                        skip_this_map = 1;
                        exit(1); /* should skip, not exit, to handle multiple maps, but not yet */
                    } else {
                        /* we found one before [0], but to use it, have to move all the
                                  trailing data UP 1 so move everything from here to end up one */

                        ti=sli;
                        bi=sli+1;
                        while (bi < HPEnd) {
                            for (i=0; i<3; i++) {
                                HP[i][ti++] = HP[i][bi++];
                            }
                        }
                        /* so we should have moved everything up by 1..?
                        and HPi should still be pointing at the same (previously occupied) el. */
                    }
                }
                HP[0][HPi] = RWC;  /* the RWC of the hit */
                HP[1][HPi] = REi;  /* the RE index */
                HP[2][HPi] += RWCincr; /* the increment caused by the too-close labels */
            }
            HP[2][HPi] += HP[0][HPi]; /*  this makes the [2] teh entire diff instead of the increment */
            HPc++;
        }
        j+=2;
    }
    /* this should end with the (size NLabels) Holding Pen (HP) populated with
       RWC's (HP[0][x]) and REi's (HP[1][x]) of REs / patterns */

    if (HPi == HPEnd ) {
        fprintf(stdout, "\n\tMore than %d patterns met crits for plasmid map.\n"
                "\tIncrease the stringency and try again.\n", NLabels);
        exit(1);
    } else {
        HPMax = HPi + 1;
    }

    if (debug > 0) {
        /*  --------------------  DEBUG----------------------   */
        fprintf(stderr, "\n\nPass %d\n", debug++);
        for (i=dstart; i<dend; i++) {
            fprintf(stderr, "\n%3d: ", i);
            for (j=0; j<3; j++) {
                fprintf(stderr, " %5ld ", HP[j][i]);
            }
        }
        /*       -----------------------DEBUG -------------------------*/
    }

    previ = 0;
    /* 2nd pass - increase spacing of labels that have mapped too closely  */
    HPi = 1;
    while (HPi < HPEnd) {
        /*    for (HPi=1; HPi<HPEnd; HPi++)  */   /* original */
        if ( HP[0][HPi] != 0 && HP[0][HPi+1] != 0) { /* skip if not 2 in a row */
            start=HPi;
            while ((HPi < HPEnd)  && ((HP[0][HPi] != 0) || ((HP[0][HPi] != 0) && (HP[0][HPi+1] - HP[0][HPi] > mdbl)))) {

                /* could also examine the run itself to detect those spots where instead of taking the
                	next position, they are remapped to a spot corresponding to their original offset -
                   this happens at the end of a long run where the last few els are actually separated
                   by a long distance  from the rest of the run, but get mapped into the run bc they
                   aren't separated by any empty els */

                previ++;
                HPi++;
            } /* how many in a row? */
            end = --HPi;
            run = end - start + 1;

            if (debug > 0) {
                fprintf(stderr, "\n at HPi= %ld, start=%d end=%d run = %d\n",   HPi, start, end, run );
                for (j=0; j<3; j++) {
                    fprintf(stderr, " %5ld ", HP[j][HPi]);
                }
            }

            /* calculate the separation betw labels to be equal */

            center = ((HP[0][end] - HP[0][start]) / 2) + HP[0][start];
            abegin = center - (BasesPerLabel * (run/2)) ; /* back up the right amount */
            if (run == 1) {
                abegin = abegin - (BasesPerLabel*0.5);
            }
            if (run > 2) {
                abegin -= BasesPerLabel;    /* if have a larger run, back up some more */
            }
            if (abegin < 0) {
                abegin = 0;
            }
            HP[2][start] = abegin;  /* + BasesPerLabel; mark the 1st at the backup position */
            for (i = start+1; i <= end; i++) {
                abegin += BasesPerLabel;
                HP[2][i] = abegin; /* they should now be spread evenly across the group spread */
            }  /* and the NEXT round should detect a collision with this one. */
            previ = HPi;
        }
        HPi++;
    }

    if (debug > 0) {
        /*  --------------------  DEBUG----------------------   */
        fprintf(stderr, "\n\nPass %d\n", debug++);
        for (i=dstart; i<dend; i++) {
            fprintf(stderr, "\n%3d: ", i);
            for (j=0; j<3; j++) {
                fprintf(stderr, " %5ld ", HP[j][i]);
            }
        }
        /*       -----------------------DEBUG -------------------------*/
    }

    /* 3rd pass  */

    i = 2;
    while (i < HPEnd) {
        if ((HP[2][i-1] == 0) && (HP[2][i] - HP[2][i-2] < BasesPerLabel)) { /*  then there was mistake */
            /* and the empty el should then carry the corrected value */
            HP[2][i-1] = HP[2][i-2] + BasesPerLabel;
            HP[0][i-1] = HP[0][i];
            HP[1][i-1] = HP[1][i];
            HP[2][i] = 0; /* and set the bad slot to 0s so we don't get dups */
            HP[0][i] = 0;
            HP[1][i] = 0;
        }
        i++;
    }


    if (debug > 0) {
        /*  --------------------  DEBUG----------------------   */
        fprintf(stderr, "\n\nPass %d\n", debug++);
        for (i=dstart; i<dend; i++) {
            fprintf(stderr, "\n%3d: ", i);
            for (j=0; j<3; j++) {
                fprintf(stderr, " %5ld ", HP[j][i]);
            }
        }
        /*       -----------------------DEBUG -------------------------*/
    }


    /* now go thru all the sites in HP and mark the labels where they're supposed to go.  for each, have to
       calculate the degree of offset, draw a tic, draw a line to the label, write the label in a small font. */
    fprintf(PS, "/Helvetica findfont 6 scalefont setfont\n");
    NH.radius = PtRadius;
    NH.xc = PtCenterX;
    NH.yc = PtCenterY;   /* init NH struct */
    NH.xet = NH.yet = 0;
    NH.alpha = NH.beta =  NH.delta = 0;


    /* ------------------------------ Major and minor tics ---------------------------------------- */
    /* make the major tics -  this should be mechanized, no */
    if (seq_len <= 10000000) {
        if (seq_len <= 10000000 && seq_len > 1000000) {
            tacg_Draw_Seq_Tics(10, 0.5, seq_len, &NH, 1, 1000000,  -1,   "/Helvetica findfont 6 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(5,  0.2, seq_len, &NH, 1, 100000,  1000000, "/Helvetica findfont 4 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(2,  0.2, seq_len, &NH, 0, 10000,  100000, "/Helvetica findfont 4 scalefont setfont", PS);
        } else if (seq_len <= 1000000 && seq_len > 100000) {
            tacg_Draw_Seq_Tics(10, 0.5, seq_len, &NH, 1, 100000,  -1,   "/Helvetica findfont 6 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(5,  0.2, seq_len, &NH, 1, 10000,  100000, "/Helvetica findfont 4 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(2,  0.2, seq_len, &NH, 0, 1000,  10000, "/Helvetica findfont 4 scalefont setfont", PS);
        } else if (seq_len <= 100000 && seq_len > 10000) {
            tacg_Draw_Seq_Tics(10, 0.5, seq_len, &NH, 1, 10000,  -1,   "/Helvetica findfont 6 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(5,  0.2, seq_len, &NH, 1, 1000,  10000, "/Helvetica findfont 4 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(2,  0.2, seq_len, &NH, 0, 100,  1000, "/Helvetica findfont 4 scalefont setfont", PS);
        } else if (seq_len <= 10000 && seq_len > 1000) {
            tacg_Draw_Seq_Tics(10, 0.5, seq_len, &NH, 1, 1000,  -1,   "/Helvetica findfont 6 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(5,  0.2, seq_len, &NH, 1, 100,  1000, "/Helvetica findfont 4 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(2,  0.2, seq_len, &NH, 0, 10,  100, "/Helvetica findfont 4 scalefont setfont", PS);

        } else {
            tacg_Draw_Seq_Tics(10, 0.5, seq_len, &NH, 1, 100,  -1,   "/Helvetica findfont 6 scalefont setfont", PS);
            tacg_Draw_Seq_Tics(5,  0.2, seq_len, &NH, 1, 10,  100, "/Helvetica findfont 4 scalefont setfont", PS);
        }
    } else {
        fprintf(stderr, "\nSequence is > 10,000,000; tics jammed too close for this now\n\n");
    }

    /* & minor tics */
    i = 1;

    NH.x1 = PtCenterX;
    NH.y1 = PtCenterY + PtRadius;



    /* ------------------------------     RE labeling     ---------------------------------------- */

    HPi = 0;
    while (HPi < HPMax) {
        NH.radius = PtRadius+3; /* '+3' to allow a band to show degeneracies */
        while (HP[0][HPi] == 0 && HPi < HPMax) {
            HPi++;    /* walk up to the next populated slot */
        }
        if (HPi < HPMax) { /*  we're not to the end, so finish ... */
            RWC = HP[0][HPi];  /* get RWC */
            n = HP[1][HPi];  /* get REi */
            /*       fprintf(stderr, "\nDEBUG:HPi=%ld RWC = %ld RE=%s RWC adj=%ld\n", HPi, RWC, RE[n].E_nam, HP[2][HPi]);  */
            /* calc the coordinates for the next hit position  */

            fprintf(PS, "\n\ngsave\n   newpath  ");
            NH.alpha =  ((double)RWC / (double)seq_len) * 360; /* calc the offset in degrees */

            Calc_Next_Hit_Coords (&NH); /* now the new coords are in NH.xn, NH.yn */
            /* place a ticmark there */
            fprintf(PS, "   %6.2f %6.2f moveto      %% RE#: %d hit=%ld, angle=%f\n",
                    NH.xn, NH.yn, n, RE[n].Sites[i], NH.alpha);   /* do a direct moveto */
            /* write the tic */
            NH.radius += 10;
            Calc_Next_Hit_Coords (&NH); /* now the new coords are in NH.xn, NH.yn */

            NH.xet = NH.xn;
            NH.yet = NH.yn;
            fprintf(PS, "   0.1 setlinewidth %6.2f %6.2f lineto stroke   %% draw tic\n", NH.xet, NH.yet);
            if (HP[0][HPi] != 0) {    /* HP[2][HPi]  or put in a very small nonzero */
                NH.alpha =  ((double)(HP[2][HPi]) / (double)seq_len) * 360; /* calc the Label offset in degrees */
                NH.radius += 15;
                Calc_Next_Hit_Coords (&NH); /* now the new coords are in NH.xn, NH.yn */
                /* offset them by the correct amount */
                /*            fprintf(PS, "   %6.2f %6.2f moveto 0.1 setlinewidth  %6.2f %6.2f lineto stroke  %% angle=%f\n",
                                           NH.xet, NH.yet, NH.xn, NH.yn, NH.alpha );    */
                fprintf(PS, "   %6.2f %6.2f moveto 1 setlinewidth  %6.2f %6.2f lineto stroke  %% angle=%f\n",  NH.xet, NH.yet, NH.xn, NH.yn, NH.alpha );
                fprintf(PS, "   %6.2f %6.2f  moveto   %6.2f rotate  ",
                        NH.xn, NH.yn, (90 - NH.alpha));
            }

            /* compose & print the label */
            if (RE[n].E_olap > 0)      {
                o = "\\\\";
            } else if (RE[n].E_olap < 0) {
                o = "/";
            } else                       {
                o = "|";
            }
            star = ' ';
            if (RE[n].E_Ncuts == 1)    {
                star = '*';
            }
            sprintf(Label, "%s [%ld] %s %c", RE[n].E_nam, RWC, o, star);

            if (NH.alpha > 180) {
                fprintf(PS, " gsave\n 180 rotate\n ");                                         /*  */
                fprintf(PS, " (%s) stringwidth neg exch neg exch rmoveto\n (%s) show\n grestore \n grestore\n \n" , Label, Label);
            } else {
                fprintf(PS, "(%s) show\ngrestore\n", Label); /*  show the label in correct position */
            }
            HPi++;
        }
    }


    /* ------------------------------       Plot ORFs       ---------------------------------------- */

    /* now if ORFs were requested, print out the ORFS greater than the requested minimum */
    fprintf(PS, "\n gsave\n");
    NH.radius = PtRadius - 30; /* start point for the ORF arcs */
    ORFNum = 0;
    ORF_density = 0.0;
    FCnt = 0;
    ORF_line_width = 3; /* temp est  */

    /* this has to be modified to do 1-3 1st and 3-6 after, compacting output so that if 1246, it will draw
       1,2 then go directly to 4 then to 6, eliminating circles for 3,5. Do this by sorting the ORFs and compact
       to be contiguous.  Also have to reverse the arrowhead symbol for frames 3-6  */

    fprintf(PS, "\n\n/Helvetica findfont 8 scalefont setfont\n");
    while (F.ORFs2Do[FCnt] != 0) { /* do all at first, then follow the F.ORFs2Do array as to which to print */
        f = F.ORFs2Do[FCnt] - 1;
        ORFNum = 0;
        NH.radius -= 10 ; /* make radius a bit smaller for each frame. */
        hairline_radius = NH.radius + 7;
        /* draw the hairline circle as a separator */
        NH.alpha = 0;
        Calc_Next_Hit_Coords (&NH);
        if (f<3) {
            arrow = ">";
        } else     {
            arrow = "<";
        }
        fprintf(PS, "\n newpath %6.2f %6.2f moveto  (%d%s) show \n", NH.xn, (NH.yn-2), (f+1), arrow);
        fprintf(PS, "\n newpath 0.1 setlinewidth 0.5 setgray\n %f %f %f  0 360  arc stroke\n",
                NH.xc, NH.yc, hairline_radius-2 );
        ORF_density += 0.07;
        /* ORF_struct ORFs contains 6 elements, populated in order..? that contain the ORFs found.
           at each ORF iter, have to move to the appro radius and draw the arcs with the appro shading
           or hatching (identifying table at the side*/
        while (ORFs[f][ORFNum].orflength > 0) { /* indicator that there are no more ORFs in this frame */
            /* calc ORF angular offsets in degrees for the beginning of the ORF */
            if (f > 2) {
                temp = ORFs[f][ORFNum].E_offsetBP;
                ORFs[f][ORFNum].E_offsetBP = ORFs[f][ORFNum].B_offsetBP;
                ORFs[f][ORFNum].B_offsetBP = temp;
            }
            NH.alpha =  ((double)ORFs[f][ORFNum].E_offsetBP / (double)seq_len) * 360; /* calc the offset in degrees */
            Calc_Next_Hit_Coords (&NH); /* now the new coords are in NH.xn, NH.yn */
            arc_end = NH.alpha;
            /* and move there */
            fprintf(PS, "\n %f %f moveto\n", NH.xn, NH.yn); /* unless we have to reset to the origin 1st  */
            /* possibly add an arrowhead here to show direction of transcription/translation */
            NH.alpha =  ((double)ORFs[f][ORFNum].B_offsetBP / (double)seq_len) * 360; /* same thing with the End */
            Calc_Next_Hit_Coords (&NH); /* now the new coords are in NH.xn, NH.yn */
            /* now write the arc from the ORF origin to the end at line width and density set already */
            fprintf(PS, "%d setlinewidth %f setgray\n newpath  %f %f %f  %f %f  arc "
                    "%% frame=%d ORF=%d begin=%ld end=%ld len=%d\n  stroke\n",
                    ORF_line_width, ORF_density,  NH.xc, NH.yc, NH.radius, 90-arc_end, 90-NH.alpha,
                    f, ORFNum, ORFs[f][ORFNum].B_offsetBP, ORFs[f][ORFNum].E_offsetBP, ORFs[f][ORFNum].orflength);

            ORFNum++;
        }
        FCnt++;
    }
    /* and plot a final hairline to finish off the ORF */
    NH.radius -= 4 ; /* make radius a bit smaller  */
    NH.alpha = 0;
    Calc_Next_Hit_Coords (&NH);
    fprintf(PS, " newpath 0 setgray 0.2 setlinewidth\n %f %f %f  0 360  arc stroke\n",
            NH.xc, NH.yc, NH.radius);

    /* ------------------------------     Plot degeneracies    ---------------------------------------- */

    fprintf(PS, "\n\n%% Marking degeneracies\n3 setlinewidth 0.0 setgray\n");
    if (F.LogDegens) {
        NH.radius = PtRadius +1.5;
        NH.x1 = NH.xc;
        NH.y1 = NH.yc + NH.radius;
        for (i=2; (i<Degen_Log[1] && Degen_Log[i] < seq_len); i+=2) {
            NH.alpha = ((double)Degen_Log[i] / (double)seq_len) * 360; /* calc the offset in degrees to the START of the degen run */

            Calc_Next_Hit_Coords (&NH);
            arc_end = NH.alpha;
            fprintf(PS, "\n newpath %6.2f %6.2f moveto newpath %% Base %ld \n",
                    NH.xn, NH.yn,                   Degen_Log[i]);
            NH.alpha = ((double)Degen_Log[i+1] / (double)seq_len) * 360; /* offset to END of the degen run */
            fprintf(PS, "  %f %f %f  %f %f  arc stroke\n\n", NH.xc, NH.yc, NH.radius, 90-NH.alpha, 90-arc_end ); /* 90-NH.alpha, 90-arc_end */
        }
    }

    fprintf(PS, "\nshowpage\n");
    fclose(PS);
}

/* int tic_size - the length of the tics in pts
   float line_width - the width of the line in pts
   long seq_len - the sequence length in base pairs
   struct NextHit_struct *NH - the Next Hit struct
   int drawlabels - whether to print labels or just write the tics
   long repeat_period - the tic period (1000 -> 1000, 2000, 3000, etc)
   int skip_period - what period should be skipped.  If already done repeat_period -> 1000 and
         doing minor tics at 100, don't want to overwrite the 1000 label
   int font_string - the font characteristics as Postscript likes it:

   FILE *PS - the file pointer to the file we're writing to
*/
void tacg_Draw_Seq_Tics(int tic_size, float line_width, long seq_len, struct NextHit_struct *NH, int drawlabels,
                        long repeat_period, int skip_period, char *font_string, FILE *PS )
{
    int i, ntics, nk;
    /*   label_offset = -35;  how far should we move towards the center to compensate for the label */
    char *Label;
    if ((Label = calloc(128,sizeof(char))) == NULL) {
        BadMem("tacg_Draw_Seq_Tics:Label", 1);
    }

    tic_size *= -1;  /* tics go towards center */
    ntics = (int)(seq_len / repeat_period);
    /*   if ((seq_len % 1000) < 0.00001) { ntics--; }  knock off the last tic if it overlaps with the origin. */
    nk=0; /* num of k's, 1000s were the orig tic counter */
    fprintf(PS, "gsave\n %f setlinewidth  %s\n", line_width, font_string);
    for (i=0; i<= ntics; i++) {  /*  && ((*NH).alpha < 359) */
        (*NH).alpha = (double)nk / (double)seq_len * 360;
        if ((*NH).alpha < 360 ) { /* && (*NH).alpha > 0 */
            fprintf(PS, "\n\nnewpath\n gsave\n ");
            Calc_Next_Hit_Coords (NH);

            fprintf(PS, " %f %f moveto \t%% tic mark: %d, angle=%f\n",
                    (*NH).xn, (*NH).yn, nk, (*NH).alpha);
            fprintf(PS, " %f rotate\n gsave\n", (90-(*NH).alpha));        /* rotate */

            if (nk > 0) sprintf(Label, "%d", nk);
            else if (skip_period < 0 ) sprintf(Label, "     0 ");
            else sprintf(Label, " ");

            fprintf(PS, " %d 0 rlineto \n stroke\n grestore\n", tic_size);

            if (drawlabels) {  /* for the TICS */
                if ((*NH).alpha > 180) {
                    fprintf(PS, " gsave\n 180 rotate\n ");
                    fprintf(PS, " (%s) stringwidth neg exch .5 mul exch  rmoveto\n (%s) show\n grestore\n grestore\n" , Label, Label);
                } else {
                    fprintf(PS, " (%s) stringwidth neg exch dup 2.5 mul sub exch  rmoveto\n",  Label);
                    fprintf(PS, "(%s) show\n grestore\n\n", Label); /*  show the label in correct position */
                }
            } else {
                fprintf(PS, " grestore\n");
            }
            nk += repeat_period;
            if (skip_period > 0 && nk % abs(skip_period) < 0.00001) {
                nk += repeat_period;    /* skip the skip period */
            }
        }
    }
    fprintf(PS, " grestore\n");    /* to restore the context for whatever was set before */
    free(Label);
}


/* passed the pointer of the populated struct, it calculates the x,y coords of the next hit and passes them
   back as struct elements xn,yn */
void Calc_Next_Hit_Coords (struct NextHit_struct *NH)
{
    double pi = 3.141592653589, reflect = 1;
    if ((*NH).alpha > 180) {
        reflect = -1;
    }
    (*NH).x1 = (*NH).x1 - (*NH).xc;
    if ( fabs((*NH).x1) > (*NH).radius) {
        /*       fprintf(stderr, "abs(x1) (%f) > radius (%f)?? \n", (*NH).x1, (*NH).radius);  */
        (*NH).x1 = (*NH).radius; /* make sure it gets reset before next iteration */
    }
    (*NH).y1 = (*NH).y1 - (*NH).yc;
    (*NH).beta =  (pi * asin((*NH).x1 / (*NH).radius)) / 180;
    (*NH).delta = 90 - (*NH).alpha - (*NH).beta;
    /*    fprintf(stderr, "\n\tbeta=%f delta=%f alpha=%f asin=%f\n",
          (*NH).beta, (*NH).delta, (90 - (*NH).beta - (*NH).delta), asin((*NH).x1 / (*NH).radius)); */
    (*NH).yn = (sin((*NH).delta / (360 / (2 * pi))) * (*NH).radius); /* don't re-add the center value yet */
    (*NH).xn = (*NH).xc + (sqrt((*NH).radius * (*NH).radius - (*NH).yn * (*NH).yn) * reflect);
    (*NH).yn = (*NH).yc + (*NH).yn; /* now add it back */

    (*NH).x1 = (*NH).xn;
    (*NH).y1 = (*NH).yn;    /* update the 'old' coords as well for the next round */
}
/* This pitiful fn() just jams the psprolog header into a string and returns it to whatever calls it
   so that we can skipp all the file manipulation and problems when the file isn't in the right place
   when we expect it */

char *tacg_PS_Prolog(void)
{
    char *psprolog = "%!PS-Adobe-2.0 EPSF-1.2\n%%Creator: \n%%CreationDate: \
    \n%%BoundingBox: 22 171 567 738\n%%EndComments\n%%BeginProcSet:Adobe_Illustrator_1.2d1 0 \
    0\n/tacg3defs dup 100 dict def load begin\n% definition operators\n/bdef {bind def} bind \
    def\n/ldef {load def} bdef\n/xdef {exch def} bdef\n% graphic state operators\n/_K { 3 \
    index add neg dup 0 lt {pop 0} if 3 1 roll } bdef\n/_k /setcmybcolor where {\n  \
    /setcmybcolor get\n} {\n    { 1 sub 4 1 roll _K _K _K setrgbcolor pop } bind\n} \
    ifelse def\n/g {/_b xdef /p {_b setgray} def} bdef\n/G {/_B xdef /P {_B setgray} def} \
    bdef\n/k {/_b xdef /_y xdef /_m xdef /_c xdef /p {_c _m _y _b _k} def} bdef\n/K {/_B \
    xdef /_Y xdef /_M xdef /_C xdef /P {_C _M _Y _B _k} def} bdef\n/d /setdash ldef\n/_i \
    currentflat def\n/i {dup 0 eq {pop _i} if setflat} bdef\n/j /setlinejoin ldef\n/J /setlinecap \
    ldef\n/M /setmiterlimit ldef\n/w /setlinewidth ldef\n% path construction operators\n/_R {.25 sub \
    round .25 add} bdef\n/_r {transform _R exch _R exch itransform} bdef\n/c {_r curveto} \
    bdef\n/C /c ldef\n/v {currentpoint 6 2 roll _r curveto} bdef\n/V /v ldef\n/y {_r 2 copy \
    curveto} bdef\n/Y /y ldef\n/l {_r lineto} bdef\n/L /l ldef\n/m {_r moveto} bdef\n% path \
    painting operators\n/n /newpath ldef\n/N /n ldef\n/F {p fill} bdef\n/f {closepath F} \
    bdef\n/S {P stroke} bdef\n/s {closepath S} bdef\n/B {gsave F grestore S} bdef\n/b \
    {closepath B} bdef\n/rightshow % stack = fontsize string\n    { dup stringwidth pop   \
    % get len of string\n      exchange sub   %calc whitespace needed\n      0 rmoveto   \
    %move over that much\n      show } def  % and show it\n\nend\n%%EndProcSet\n%%EndProlog\n%%Page: 1 \
    1\ntacg3defs begin\n";
    return psprolog;
}


/* tacg_Eval_Rules_from_File takes care of:
	- reading rules from the file (path validated via SearchPaths()),
   - parsing the file into the Rules struct
   - validating that the rules are formed correctly (parens match up, logicals are correct)
   - making sure that each of the subpattern names has been loaded in RE (have to hash the
   	names as in satisfying the -x flag in ReadEnz().
   - eventually, create --rules 'thisrule,thatrule,mars,promoter, etc'
     (can be specified w/o --rulefile in which case it will read the deafult 'rules.data'.)
     as opposed to current --rule for specifying a rule at a time from the commandline.
     (which should write out that rule to the tacg.patterns file for incoporation into the
     main rules file.  require that it be named as well.
     ie: --rule 'Name,((()()())(())),width'

     - open the file, getting path name from F.RuleFile
     - read in the strings, concatting as we go, then validate the logic
         and stuff all the bits into the Rule_struct
     - at each iter, make sure there's enough mem
     - when finished with the file read, and all the bits are in Rule_struct, go thru
     		Rules[].Name and check against the names in RE[].E_nam.
     - make
 */
int tacg_Eval_Rules_from_File(struct SE_struct SelEnz[MAX_SEL_ENZ],  char datestamp[80])
{
    /*   long ProtoNameHash[], int NumREs, */

    FILE *fpRF;
    char junk[10], *ctemp1, ct, *rulestring, *packed_RS, *rstr;
    short NameCheck[RE_NAME_TAB_SIZ]; /* indices from realhash(PatName) to check against multiple identicals  */
    int e=0,SEi=0;
    long j, EOpRS=0,
            window_size = 0,
            NumRules=0,
            rstr_len, linenum=0,
                      RS_Cnt = 0, /* RuleStruct Counter */
                      EOAM_RS = INIT_NUM_RULES; /* End Of Allocated Memory for the Rule_Struct; allocated in tacg */

    short cont = 0; /* continue indicator */

    /* alloc the Rule_struct Rules */
    if ((Rules = calloc(EOAM_RS, sizeof(*Rules))) == NULL) {
        BadMem("Can't init mem for Rules", 1);
    }

    if ((ctemp1 = calloc(512,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'ctemp1'", 1);
    }
    if ((rulestring = calloc(MAX_LEN_RULESTRING,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'rulestring'", 1);
    }
    if ((packed_RS = calloc(MAX_LEN_RULESTRING,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'packed_RS'", 1);
    }
    if ((rstr = calloc(MAX_LEN_RULESTRING,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'rstr'", 1);
    }
    memset(NameCheck,0, sizeof(short)*RE_NAME_TAB_SIZ); /* make sure everything is set to 0 */

    if ((fpRF = fopen(F.RuleFile,"r")) == NULL) {  /* or the standard/optional one from SetFlags() */
        fprintf(stderr,"Cannot open the Rules file \"%s\" for reading!!\n", F.RuleFile); /* DON'T HIDE behind -V */
        exit(1); /* print an error and die gracefully */
    } else if (F.Verbose > 0) {
        fprintf(stderr,"Rule File %s opened OK for reading\n", F.RuleFile);
    }

    /*  Scan thru comments in rebase file to separator ("..")   */
    junk[0] = 'k';   /* make sure the strcmp below fails the first time */
    while ((feof(fpRF) == 0) && (strncmp(junk, "..",2) != 0))  {
        e = fscanf (fpRF, "%2s", junk);
        if (e <= 0) {
            fprintf(stderr, "SeqFuncs.c: 744 Failure to read\n");
            exit(1);
        }
        ct = 'm'; /* and read to end of line - this needs to be able to deal */
        while (ct != '\n' && ct != '\r' && (feof(fpRF) == 0)) {  /* with end of line conditions for Mac, Win, and unix - what have I missed? */
            ct = fgetc(fpRF);
            linenum++;
        }
    }

    /* we're at "..", so go into file and start slurping directly data into struct
       except for RE_rawsite which has to be filtered a couple ways before being acceptable */

    ct = 'm';
    while (fgets(rstr, 512, fpRF) && (NumRules <= MAX_NUM_RULES)) {
        /* check the size of the rule struct to see if it's going to run into the end of the alloc'ed mem */
        if (RS_Cnt == EOAM_RS) { /* then we need to alloc more mem for the Rules struct */
            EOAM_RS +=  RULES_INCR;
            if ((Rules = realloc(Rules, sizeof(*Rules)*EOAM_RS)) == NULL) {
                BadMem("Failed to get more mem for Rules.\n", 1);
            }
        }

        /* now the whole line's in rstr; use strsep to parse it, make decisions re: parsing, skipping */
        j = stripwhitespace(rstr,0,0);
        if (rstr[0] != ';'  && (strlen(rstr) != 0)) {  /* if the line hasn't been commented out, get all the bits */
            /* REMEMBER!  strsep chews the string to pieces and changes the start of the string */
            if (j>0 && F.Verbose > 0) {
                fprintf(stderr,"DEBUG: stripped %ld whitespace characters from RuleName\n", j);
            }

            strcpy(packed_RS, rstr);  /* copy over to pass to tacg_SlidingWindow - may not be needed */


            /* now have to break on ',' to get all the rulestring in 1 shot
            	- also have to detect if the last char is a continuation char '\' */

            /* if last char is '\', then rulestring is incomplete and there won't be a window term  */
            cont = 0;
            EOpRS = strlen(packed_RS) - 1; /* */
            if (packed_RS[EOpRS] == '\\') { /* then there's a continuation char */
                cont = 1;
                /* keep on going until there's no more continuation lines */
                while ((cont == 1) && fgets(rstr, 512, fpRF)) {  /*     */
                    j = stripwhitespace(rstr,0,0); /* no whitespaces, eols */
                    rstr_len = strlen(rstr); /* length excluding terminating \n */
                    EOpRS = strlen(packed_RS);
                    if (EOpRS + rstr_len >= MAX_LEN_RULESTRING) {
                        fprintf(stderr, "Total length of the rulestring exceeded; change MAX_LEN_RULESTRING \n"
                                "or check your continuation characters for the rule %s\n", ctemp1);
                        exit (1);
                    }
                    /* rstr should be the right size to copy mod the possible final \ and \n */
                    strcpy(packed_RS+EOpRS-1,rstr); /* copy all but the last char over */

                    if (rstr[(strlen(rstr)-1)] == '\\' ) {
                        cont = 1;
                    } else                                {
                        cont = 0;
                    }
                }
                if ((rstr = calloc(MAX_LEN_RULESTRING,sizeof(char))) == NULL) {
                    BadMem("Can't init mem for 'rstr'", 1);
                }
                strcpy(rstr,packed_RS);

            }
            /* coming out of the continuation loop, packed_RS contains the entire Name,RuleString,Window concat
               that the whole thing should look like. rstr has only the last line of a continued rule, so it needs to
               be synced, so copy packed_RS to rstr again */

            ctemp1 = strsep(&rstr, ","); /* using ctemp1 (=Name) for congruence with previous code  */
            if (strlen(ctemp1) > MAX_PAT_NAME_LEN) {
                fprintf(stderr, "Name %s is longer than allowed (%d) - check it!\n", ctemp1, MAX_PAT_NAME_LEN);
                exit(1);
            }


            /* so now packed_RS has the whole packed rulestring in it at this point; rstr has everything
               BUT the RuleName in it. So now we finish off the line, so break on ',' to get the rulestring
               and ONLY THEN break on the ',' to separate the rule term from the sliding window term
               (since already have processed the Name) */
            if (strlen(rstr) >= MAX_LEN_RULESTRING) {
                fprintf(stderr, "Total length of the rulestring exceeded; change MAX_LEN_RULESTRING \n"
                        "or check your continuation characters for the rule %s\n", ctemp1);
                exit (1);
            }

            /* break on the comma so we can ignore the window term - this probably needs some more error-checking */
            memset(rulestring, 0, MAX_LEN_RULESTRING);
            rulestring = strsep(&rstr,",");  /* rulestring now contains ONLY the packed rulestring */
            window_size = (long)atoi((char *)strsep(&rstr,","));

            if (window_size < 1) {  /* now test whether windowsize & rulestring are valid  */
                fprintf (stderr, "Window size for %s is less than 1, which is silly.\n", ctemp1);
                exit (1);
            }
            /* Are the component pattern Names valid? ie is this name represented in RE by being loaded
               from the RE database?  need to hash the Names of the individual patterns and see if they
               hit any of the Names that were loaded from whatever rebase was being used.
               Soooooo, go thru the rulestring and hash all words breaking on ()|^& to cut out the 'Name:min:Max' stanzas */
            tacg_Fill_SelEnz_fr_RuleString(SelEnz, &SEi, rulestring, NameCheck, ctemp1);



            /* stuff ALL the bits (Name, rulestring, window_size) into the struct (the validity of the rule itself
               tested in the evaluation code) */
            strcpy(Rules[RS_Cnt].Name, ctemp1);
            if ((Rules[RS_Cnt].RuleString = calloc((strlen(rulestring)+1),sizeof(char))) == NULL) {
                BadMem("Can't init mem for 'Rules[RS_Cnt].RuleString'", 1);
            }
            strcpy(Rules[RS_Cnt].RuleString, rulestring);
            Rules[RS_Cnt].Window = window_size;

            if ((rstr = calloc(MAX_LEN_RULESTRING,sizeof(char))) == NULL) {
                BadMem("Can't init mem for 'rstr'", 1);
            }
            RS_Cnt++;
        } else {
            /*          fprintf(stderr, "\nDEBUG: line ignored: %s\n", rstr); */
        }
    }
    F.Xplicit = -1; /* notify ReadEnz */
    /*   F.Rules = RS_Cnt-1;  */  /* this should NOT have to be decr by one -
            for F.Rules:
            0 = no rules
            1 = rules entered by cmdline
            2 = rules entered by file */
    F.Rules = 2; /* does this have to be decr by one 1st?? */
    F.NumRules = RS_Cnt;
    return SEi;
}



void BadMem(char* message, int exitcode)
{
    fprintf(stderr, "Error on memory call: %s\n", message);
    exit(exitcode);
}



/* BestHexWWPal() does much of the calculation to set the RE values related to which way the extended check
   should look (WW), what the best hexamer is, what the center of the site is, etc.  */
void BestHexWWPal(char *RE_rawsite, int len, int Cur)
{
    int i, j, dpl=0, ori=1, BestHexIndex=0;
    float f=0, MaxSoFar=0;

    /* for sites > 6, calc the sliding window of magnitudes to pick the best one */
    MaxSoFar = 0;
    if (len > 6) {
        for (i=0; i<len-5; i++) {
            f = MagCalc(RE_rawsite+i, 6);
            if (f > MaxSoFar) {
                MaxSoFar = f;
                BestHexIndex = i;
            }
        } /* so now the best hex is chosen - now copy it in */
        RE[Cur].WW[ori] = 0 - BestHexIndex;
        for (j=0; j<6; j++) RE[Cur].E_hex[ori][j] = RE_rawsite[j+BestHexIndex]; /* copy the best hex to RE_hex - */
    } else { /* the site is < 6 bases, so it has to be filled with n's */
        i=j=0;   /* this is from previous code - not great but it works */
        while (j<6) {
            if (RE_rawsite[i] == '\000') { /* if the site is shorter than a hexamer */
                while (j<6) RE[Cur].E_hex[ori][j++] = 'n'; /* pad out the site with n's */
            } else RE[Cur].E_hex[ori][j++] = RE_rawsite[i++]; /* else copy one to the other */
        }
        RE[Cur].WW[ori] = 0;
    }
    if (RE[Cur].E_pal == 0)  {
        /* if the RE is not a pal, do the things common to the normal ops
            and the ERRORs ops; do the rest below  */
        Anti_Par(RE[Cur].E_hex[ori], RE[Cur].E_hex[dpl], 6); /* antipar and pop in the hex  */
        /* topcut, WW require a little jig */
//		  if (F.Verbose > 1) fprintf(stderr, "BestHexWWPal: before: for RE=%d, E_olap=%d E_tcut[0]=%d, [1]=%d\n", Cur, RE[Cur].E_olap, RE[Cur].E_tcut[0], RE[Cur].E_tcut[1]);
        RE[Cur].WW[dpl] = 0 - RE[Cur].E_len + 6 - RE[Cur].WW[ori]; /* WW is negative */
        RE[Cur].E_tcut[dpl] = RE[Cur].E_len - RE[Cur].E_tcut[ori] - RE[Cur].E_olap;    /*  tcut IS MIDPOINT OF  *besthex*  */
//		  if (F.Verbose > 1) fprintf(stderr, "BestHexWWPal: after: for RE=%d, E_olap=%d E_tcut[0]=%d, [1]=%d\n", Cur, RE[Cur].E_olap, RE[Cur].E_tcut[0], RE[Cur].E_tcut[1]);
    }  /* end partial nonpalindrome handling */
}


/* tacg_Fill_SelEnz_fr_RuleString() takes the SelEnz struct and the rulestring (and a few other housekeeping vars
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
       and the ending SEi.
*/

void tacg_Fill_SelEnz_fr_RuleString(struct SE_struct SelEnz[MAX_SEL_ENZ], int *SEi,
                                    char *rulestring, short NameCheck[RE_NAME_TAB_SIZ],  char *ctemp1)
{
    char *Name_mM, *rstr, *PatName = "            ",
                           *amin = "           ",
                            *aMax = "           ";
// *Name_mM_addr not needed until need perfect memory cleanup.

    int ind;
    if ((Name_mM = calloc(33,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'Name_mM'", 1);
    }
    if ((rstr = calloc(MAX_LEN_RULESTRING,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'rstr'", 1);
    }
//    rstr_addr = rstr;
    strcpy(rstr,rulestring); /* re-using rstr; mem should be fine after re-calloc above */
//   fprintf(stdout,"DEBUG: rulestring = %s, rstr = %s\n", rulestring, rstr);
    Name_mM[0] = '\0'; /*  to pass on 1st pass */
    while (Name_mM[0] == '\0') {
        Name_mM = strsep(&rstr, "()|^&");
    } /* this should break out a 'Name:#:#' set */

    while (Name_mM != 0L) {
//   	fprintf(stdout,"DEBUG: Name_mM = %s\n", Name_mM);
        if (Name_mM[0] != '\0') {
            PatName = DownCase(strsep(&Name_mM, ":"));
            ind = (int) realhash(PatName, RE_NAME_TAB_SIZ);
//          fprintf(stderr, "DEBUG: realhash '%s' = %d\n", PatName, ind);
            if (strlen(PatName) <= MAX_PAT_NAME_LEN) {
                if (NameCheck[ind] == 0) {
                    strcpy(SelEnz[*SEi].PName,PatName);
                    NameCheck[ind] = 1; /* now mark it to make sure that we don't take it again */
                    if ((amin = strsep(&Name_mM, ":")) == '\0') {
                        SelEnz[*SEi].min =   0; /* RE[ii].min =  */
                    } else {
                        SelEnz[*SEi].min =  atoi(amin); /* RE[ii].min =  */
                    }
//               fprintf(stdout,"DEBUG: min = %d\n", SelEnz[*SEi].min);
                    if ((aMax = strsep(&Name_mM, ":")) == '\0') {
                        SelEnz[*SEi].Max = 32000; /* RE[ii].Max = */
                    } else {
                        SelEnz[*SEi].Max =  atoi(aMax); /* RE[ii].Max =   */
                    }
//               fprintf(stdout,"DEBUG: Max = %d\n", SelEnz[*SEi].Max);
                    (*SEi)++;

                } else { /* else it's already been seen in another rule or sub-rule and we don't care */
                    /*                fprintf(stderr, "DEBUG: The pattern named %s has already been seen previously; continuing .. \n", PatName); */
                }
            } else {
                fprintf (stderr, "\nERR: In Rule named %s, Pattern Name %s is too long (must be <=10 chars)\n",
                         ctemp1, PatName);
                exit(1);
            }
        }
        Name_mM = strsep(&rstr, "()|^&");  /* get the next one */
    }/* @ end of loop, all constituent patterns are loaded into SelEnz or the app has exited on an err */
    // original rstr is chewed to nothing; have to use a copy of the address to free it
    //free(rstr_addr);
    //free(Name_mM_addr);
} // end tacg_Fill_SelEnz_fr_RuleString()


/* ***************************************************************************************************** */





/* RE_Filter returns 1 if the RE indexed by Proto[REi] matches the tests for
	Magnitude
   Number of Cuts
   Cost
   Overlap
   etc that are can be used to fine-tune the selection of REs */
int RE_Filters_OK(int REi)    /* REi should be the ProtoCutters index, as that is the only RE el that we have to check */
{
    int OK=0;
    /*    fprintf(stdout, "\n REi = %d, RE name = %s \n", REi, RE[REi].E_nam); */
    if (((RE[REi].E_mag   >= (int)F.Mag)       &&    /* E_mag must be >= to -n flag */
            (REi >= F.RE_ST)                              &&    /* Don't EVER want dam/dcm sites printed */
            (RE[REi].E_Ncuts >= (int)F.Min)       &&    /* Ncuts >= -m, */
            (RE[REi].E_Ncuts <= (int)F.Max)       &&    /* Ncuts <= -M */
//        (RE[REi].Cost    >= (int)F.Cost)      &&
            (RE[REi].Cost    >= (int)F.Cost))      &&
            /* the dam/dcm part of this COULD be rolled in if the RWC (aka base_real_pos) was included as an optional
            	argument (or if -1, ignore it), but right now it doesn't seem to be used in the same way, so I'll let
               the code stay but commented out for now */
            /*        (DamDcmOK(base_real_pos, REi) == 1)   &&     print only RE's !masked by dam/dcm methylases */

            (F.Overhang == 1      || /* if want all, only 1st has to be true to short circuit */
             ((RE[REi].E_olap > 0) && (F.Overhang == 5)) || /* if olap < 0 && overhang = 5' */
             ((RE[REi].E_olap < 0) && (F.Overhang == 3)) || /* if olap > 0 && overhang = 3' */
             ((RE[REi].E_olap ==0) && (F.Overhang == 0)))) { /* if olap = 0 && overhang = blunt */
        OK = 1;
    } else {
        OK = 0;
    }
    if (F.Regex != 0) { OK = 1;}
    if (F.Verbose >= 2 && OK == 0) fprintf (stderr, " %s failed RE_Filters_OK()\n", RE[REi].E_nam);
    return OK;
}

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
   Hmmm - this could be the prototype for an overall RE filter to be used in a number of cases..
*/

/* Clone() goes thru RE[] and figures out which REs match the conditions of the option flag, where
   --clone '#_#,#x#,#x#,#_#' (up to MAX_NUM_CLONE_RANGES).  '_' indicates that the range cannot be cut;
   'x' indicates that it can be cut.
   There will probably be some interactions between --clone and other options, but
   for now, finish the basic functionality and sort out the interactions later.  In fact,
   RE[].. should not be a problem as it will be populated only with sites that are requested,
   so that if only a subset or special file of REs are used, only they will be used for
   for the clone() analysis.
   This fn() is called AFTER the digest is done in the normal way and should only
   need the number of REs (NumREs - not modified by this so it doesn't have to be
   passed as a pointer) as a max to go thru to determine the ranges that agree with the
   ranges specified.
*/
/* should really be running thru Protos[] from 0->ProtoCutters */
void Clone(int ProtoCutters, int *Protos)
{
    long i, c, j, m, p, r,
         Best[F.Clone+2][ProtoCutters], 	/* Best = array that stores and sorts the best REs for a set of cloning criteria
   												last el is for storing summary truths to check if ALL rules have been met */
         BadProtos[ProtoCutters];  	/* BadProtos keep track of Protos that map inside a forbidden zone and therefore have to
   										be removed from contention after the fact.  If a Protos index shows up in BadProtos,
                                 it cannot be a possible cloning RE (unless cloner wants to try a partial digest.. */
    short N_Crits;
    int LadderProtos[ProtoCutters];

    N_Crits =  F.Clone; /*  simplicity; F.Clone stores the # of ranges in the --clone option string */
    for (i=0; i< F.Clone+2; i++) {
        for (j=0; j<ProtoCutters; j++) {
            Best[i][j] = 0;
        }
    }
    memset(BadProtos,0,sizeof(long)*ProtoCutters);
    memset(LadderProtos,0,sizeof(int)*ProtoCutters);

    /* Need additional rules to elim REs.
    	If an RE hits in the forbidden zone, that RE is then added to a FORBIDDEN array to be checked post-processing.
    	now #x# = MUST HIT range  Should there be an option to set prefs for this? Or just show how many
       	there are in that zone?
       How should REs that hit in zones NOT explicitly marked as allowed (x) or forbidden (_) be treated?  In typical
       	cloning, don't want lots of extra cuts as it increases ends to contend with.
    */

    for (r=0; r<=N_Crits; r++) { /* for each range parsed out in SetFlags() */
        for (i=0; i<ProtoCutters; i++) {
            p = Protos[i]; /* p is the alias for the RE prototype index, so RE[p].etc will hit only prototypes */
            if (Clone_S[r].to_cut == 1) { /* if this is a range where we REQUIRE hits... */
                m = Clone_S[r].matching_REs[1]; /* [1] is the 1st free position of the array */
                for (j=1; j<RE[p].Sites[0]; j++) { /* RE[p].Sites[0] holds the index of the 1st FREE position in the Sites array */
                    if (RE[p].Sites[j] >= Clone_S[r].begin &&     /* site must match BOTH conditions */
                            RE[p].Sites[j] <= Clone_S[r].end ) {
                        /* do we need a per-range log for positives & negatives? */
                        Clone_S[r].matching_REs[m++] = i;  		/* incr m, store the Proto index */
                        Clone_S[r].matching_REs[1] = m;      	/* keep track of the 1st free el as well */
                        /* and check for realloc conditions */
                        if ((m+1) == Clone_S[r].matching_REs[0]) {
                            realloc_Clone_matching_REs(m, r);
                        }
                    } /* else we silently ignore the bad hits */
                }
                Clone_S[r].matching_REs[1] = m; /* and then store the position back into the array */
            } else { /* this is a range we where we will NOT TOLERATE hits */
                m = Clone_S[r].matching_REs[1]; /* [1] is the 1st free position of the array  */
                for (j=1; j<RE[p].Sites[0]; j++) { /* RE[p].Sites[0] holds the index of the 1st FREE position in the Sites array */
                    if (RE[p].Sites[j] <= Clone_S[r].begin ||      /* site can match either condition */
                            RE[p].Sites[j] >= Clone_S[r].end ) {
                        Clone_S[r].matching_REs[m++] = i;  /* incr m, store the Proto index */
                        Clone_S[r].matching_REs[1] = m;            /* keep track of the 1st free el as well */
                        /* and check for realloc conditions  */
                        if ((m+1) == Clone_S[r].matching_REs[0]) {
                            realloc_Clone_matching_REs(m, r);
                        }
                    } else { /* else we have to track the REs to mask out the results by adding their index to BadProtos[] */
                        BadProtos[i]++; /* if an el is > 0, it's EVIL */
                    }
                }
                Clone_S[r].matching_REs[1] = m; /* and then store the position back into the array */
            }
        }
    }

    /* now we've run thru all the criteria for judging whether the REs map into the ranges
     now have to output the data in some usable form.  1st, are there any REs that match ALL requirements?
    2nd output all ranges with all matching REs, sorted alphabetically.  Should --clone automatically set
    --sites, so can output all the sites at the same time??  Hmmm...
    */
    /* 1st are there any that match ALL crits? for all crits, are there any indices that appear in all of them??
    	set up array of length N_Crits and increment each RE index as go thru all the ranges.  Also has the effect of
       being able to sort those REs that match the MOST crits as well. */
    for (i=0; i<=N_Crits; i++) {   /*  cycle thru all the crits */
        for (j=2; j<Clone_S[i].matching_REs[1]; j++) { /* cycle thru all the prototype REs that match for each crit */
            Best[i][Clone_S[i].matching_REs[j]]++; 	/* this will count the # of times a proto came up positive in the criteria  */
        }
    } /*  so now have all the REs counted. */

    /* and print the results in some kind of visually appealing way. */
    /* 1st print out any REs that met ALL crits.. */
    fprintf(stdout,"\n\n== Here's the Clone output for the %d Criteria:\n", (F.Clone+1));
    for (i=0; i<=F.Clone; i++) {
        fprintf(stdout, "    Rule %ld:  %s\n", i+1, Clone_S[i].range_rule);
    }

    for (j=0; j <= F.Clone; j++) {  /* sum the rules results to generate a summary */
        for (i=0; i<ProtoCutters; i++) {
            if (Best[j][i] > 0 && BadProtos[i] == 0) {
                Best[F.Clone+1][i]++; /* for each rule passed, incr the truth element once */
            }
        }
    }
    /* and finally print out the Protos that match ALL the rules BUT this should probably be mod to allow selection and
    	filtering by the standard rules -o -n -m -M --cost, etc but later... this will require a test fn() that tests all the
       Protos against all the various restrictions, but it will be applicable to at least 4 cases */
    fprintf(stdout,"\nREs that met ALL criteria:\n");


    /*	DEBUG STANZA *****
       for (i=0; i<ProtoCutters; i++) {
          fprintf(stdout, "F.Clone+1 = %d, i = %d, Best = %d\n", F.Clone+1, i, Best[F.Clone+1][i]);
          if (Best[F.Clone+1][i] == F.Clone+1) {
       	   fprintf(stderr,"\n i= %d RE = %s", i, RE[Protos[i]].E_nam);
          }
       }
    */
    i = c = m = 0;

    while (i<ProtoCutters) {
        if (m >= F.Width) {
            m = 0;
            fprintf(stdout, "\n");
        }
        if (Best[F.Clone+1][i] == F.Clone+1 ) {
            fprintf(stdout,"%*s", F.MaxPatNameLen, RE[Protos[i]].E_nam);
            LadderProtos[i] = Protos[i];
            c++;
            m += 10;
        }
        i++;
    }
    if (c == 0) {
        fprintf(stdout,"There were no REs that met ALL criteria\n");
    } else {
        PrintGelLadderMap(i, F.Seq_Len, LadderProtos, 2);
    }
    memset(LadderProtos,0,sizeof(int)*ProtoCutters);

    fprintf(stdout,"\n\nResults (listed by Rules) that met SOME criteria & did not violate a hits-forbidden range.");
    for (j=0; j<=F.Clone; j++) {
        fprintf(stdout, "\n\n--- For Rule %ld (%s), the following REs matched:\n\n", j+1, Clone_S[j].range_rule);
        i = c = m = 0;
        while (i<ProtoCutters) {
            if (m >= F.Width) {
                m = 0;
                fprintf(stdout, "\n");
            }
            if (Best[j][i] > 0 && BadProtos[i] == 0) {
                fprintf(stdout,"%*s", F.MaxPatNameLen, RE[Protos[i]].E_nam);
                LadderProtos[i] = Protos[i];
                c++;
                m += 10;
            }
            i++;
        }
        if (c == 0) {
            fprintf(stdout,"There were no REs that met EVEN SOME of the criteria\n");
        } else {
            PrintGelLadderMap(i, F.Seq_Len, LadderProtos, 2);
        }
        memset(LadderProtos,0,sizeof(int)*ProtoCutters);
    }

    fprintf(stdout,"\n\n--- Results (listed by RE) that met SOME criteria & did not violate a hits-forbidden range\n");
    fprintf(stdout,"\n   RE Name");
    for (j=0; j<=F.Clone; j++) {
        fprintf(stdout,"  Rule %ld", j+1);
    }
    fprintf(stdout,"\n");
    for (i=0; i<ProtoCutters; i++) {
        if (BadProtos[i] == 0) {
            fprintf(stdout,"\n%*s", F.MaxPatNameLen, RE[Protos[i]].E_nam);
            for (j=0; j<=F.Clone; j++) {
                if (Best[j][i] > 0) {
                    fprintf(stdout,"    X   ");
                } else {
                    fprintf(stdout,"        ");
                }
            }
        }
    }
    fprintf(stdout,"\n\n");
};


void realloc_Clone_matching_REs(long current, long index)
{
    long newsize;
    /* and check for realloc conditions */
    if ((current+1) == Clone_S[index].matching_REs[0]) { /* if we're nearly up to the alloc'ed mem, need to realloc */
        newsize = Clone_S[index].matching_REs[0] + CLONE_INCR_MRE; /* this is how big the new size should be  */
        if ((Clone_S[index].matching_REs = (long *) realloc(Clone_S[index].matching_REs, sizeof(long) * newsize)) == NULL) {
            BadMem("Can't realloc mem for Clone_S[].matching_REs", 1);
        } else {
            Clone_S[index].matching_REs[0] = newsize;    /* new alloc range */
        }
    }
};


/* the example function included below is a mock-up to show how to include your own function.
	Obviously, there are 10s of functions included with this source code and some of them may be
   immediately applicable to your work.  Feel free to bash them into a useful state for yourself.
   If you think that you've created a widely useful function that you would like to have included
   with tacg, please give me a shout, and I'll include a pointer to it or provide other information
   about it on my web site.

   void Example(void){

	}
*/


/* LinearMap() is simply a functionization of the Linear Mapping chunk that
	should have been functionized long ago */

void LinearMap(long seq_len, int basesPerLine, int NumREs, char *sequence, int max_okline,
               char Codons[MAX_ORGS][N_CODONS][7], int codon_Table, int *Protos)
{

    /* Declarations */
    int i, d, k, mm=0, rc, n_letters=0, Xn=0, okline, HTML, spacer, spacermod,
                 min_okline, tic_line, Cur_RE, block_Cnt, block_repeat, block_cut_pos,
                 lastline = O_SEQ_LINE+2, hits, ok2wr[MAX_OUTPUT_LINES], p_width, abcp,
                 in_OL=-1, seq_offset, tripclick=0;

    long  *BackCheck, /* BackCheck[] checks proto entries in the Linear map output to remove duplicates */
          DCnt, Cnt_reset=0, base_real_pos, numstart;

    char  ctemp1[30], *s,  Prot_out[MAX_BASES_PER_LINE],
          O_txt[MAX_OUTPUT_LINES][O_LMARG+O_RMARG+MAX_BASES_PER_LINE];

    /*    char Codons[MAX_ORGS][N_CODONS][7]; */
    char os = '~';  /* started out as '<' to mark other strand, but causes problems with HTML generation */

    BackCheck = (long *)calloc((size_t)NumREs+2, sizeof(long));  /* Needs to be '+2' as combination cut will add one and 1 extra */
    if (BackCheck == NULL) BadMem("Can't init mem for BackCheck", 1);

    s = (char *) calloc(256, sizeof(char));
    if (s == NULL) BadMem("failed to calloc mem for 's'", 0);
    memset(O_txt,' ',(size_t)((O_LMARG+O_RMARG+MAX_BASES_PER_LINE)*MAX_OUTPUT_LINES));  /* set all of O_txt to blanks - works!*/
    p_width = O_LMARG + O_RMARG + basesPerLine;

    /* Routine to print out the text of the digestion, with names going where they should go
    a la Strider - have to keep track of how long the names are and where an allowable place
    to write the name - */

    /* In brief, the following function, does this:
    It reads thru DD->Dat (the array that keeps track of the linear sequence of cuts in the sequence)
    and matches the position of the enzymes to the current block being printed, taking into account
    the cut offset of the RE from the 1st base of the recognition sequence.  The tricky part is that
    because of this latter functionality, it has to check the sequence that will be the *following*
    block to see if there are any RE recognition seqs within the 1st few bases that will cause the actual
    cut site to be in the *previous* block.  Doing this seems to take a great amount of the cpu time - can
    probably be optimized a good deal more, but it's still >10 faster than GCG ;). */

    HTML = F.HTML;         /* simplicity */
    n_letters = F.Lab1or3; /* simplicity */
    Xn = F.XLinFr;         /* simplicity */
    numstart = F.NumStart;  /* simplicity */

    fprintf(stdout, "\n\n");
    if (HTML > -1) fprintf(stdout, "<H3><A NAME=\"linmap\">");
    fprintf(stdout, "== Linear Map of Sequence:\n"); /* Write the header */
    if (HTML > -1) fprintf(stdout, "</A></H3>\n");

    /* inits, etc */
    if (seq_len % basesPerLine == 0)  block_repeat = (int) seq_len/basesPerLine; /*  # times thru the compose/print cycle */
    else  block_repeat = (int) seq_len/basesPerLine + 1;
    DCnt=2;  /* start the DD->Dat counter at the beginning of REAL DATA [0]=Cnt; [1]=Siz */
    base_real_pos = DD->Dat[DCnt];

    /* spacer dictates the line spacing between the 1st 3 and 2nd 3 translation stanzas; if 6 there
        will be no space between the 2 stanzas, if 7 there will be 1 space; could make this another
        option but sheesh... */
    spacer = 7;
    spacermod = 0;
    if (F.Strands == 1) {
        spacermod++;
    }
    if (F.Tics == 0) {
        spacermod++;
    }
    spacer = spacer - spacermod;
//    if (spacer == 7) stanzaspace = spacermod + 1;

    for (block_Cnt=1; block_Cnt <= block_repeat; block_Cnt++) {
        /* hits is set to 1 if there are hits on the line; if no hits, decr the top buffer line to output
        	so that there is a constant 1 blank line spacer between stanzas */
        hits = -1;
        /* if there was an OVERLAP failure previously, fix it */
        if (in_OL == 1) {
            memset(BackCheck,0,sizeof(long)*(NumREs+2)); /* clear BackCheck[] as it's filled with stale data */
            DCnt = Cnt_reset;  /* and reset DCnt to that pointer */
        }
        in_OL = -1;  /* and reset 'in OVERLAP' marker for new block */
        base_real_pos = DD->Dat[DCnt];  /* Make sure that ' base_real_pos' is set */
        okline = min_okline = O_SEQ_LINE-2; /* set the min (highest on page) line of buffer to print out */

        /* Set up the seq and its rev compl in the output buffer */
        memset(ok2wr,0,sizeof(int)*MAX_OUTPUT_LINES);  /* set/reset all of ok2wr to 0's for the new block*/
        /* next line was deleted in a brief fork of tacg.c - note this if get unexlained segfs */
        memset(BackCheck, 0, sizeof(long)*NumREs); /* reset the BackCheck[] */

        /*write the correct numbering in Left margin  */
        if (F.NumStart < 1 && (((block_Cnt-1)*basesPerLine) + 1 + F.NumStart) == 0) {
            numstart++;
            tripclick++;
        }
        sprintf((char *)s, "%7ld",(((block_Cnt-1)*basesPerLine) + 1 + numstart)); /* numstart is just a + offset */
        memcpy(&O_txt[O_SEQ_LINE][0], s, 7);
        /* and in the right margin - combine with the above bit? */
        if (F.NumStart < 0 && (tripclick || (((block_Cnt)*basesPerLine + numstart)) > 0)) {
            if (!tripclick) {
                tripclick++;
                numstart++;
            }
        }
        sprintf((char *)s, "%7ld",(((block_Cnt)*basesPerLine + numstart)));
        memcpy(&O_txt[O_SEQ_LINE][O_LMARG+basesPerLine], s, 7);
        /* and drop in the sequence directly after it. */

        if (block_Cnt != block_repeat) k = basesPerLine; /* test for end of sequence to end gracefully */
        else k = (int)seq_len - (block_Cnt-1)*basesPerLine;
        if ((F.SeqBeg != 1 || F.NumStart != 0)  &&
                (F.Strands != 1 || F.Tics != 0)) { /* if subsequence or dif #ing, print the #s of the original seq too */
            /* write the correct numbering in Left margin  */
            sprintf((char *)s, "%7d",((int)(((block_Cnt-1)*basesPerLine)+F.SeqBeg))); /* these 2 lines could be combined...? */
            memcpy(&O_txt[O_SEQ_LINE+1][0], s, 7);
            /* and in the right margin - combine with the above bit? */
            sprintf((char *)s, "%7d",(int)(((block_Cnt)*basesPerLine - 1 + F.SeqBeg)));
            memcpy(&O_txt[O_SEQ_LINE+1][O_LMARG+basesPerLine], s, 7);
        }

        memcpy(&O_txt[O_SEQ_LINE][O_LMARG], &sequence[(block_Cnt-1)*basesPerLine+BASE_OVERLAP+1], (size_t)k);

        if (F.Strands == 2) { /* if we want the bottom strand as well */
            Rev_Compl(sequence+BASE_OVERLAP+1+((block_Cnt-1) * basesPerLine), s, (long)k);
            /* and then plunk the converted sequence into O_txt with a similar statement */
            memcpy(&O_txt[O_SEQ_LINE+1][10], s, (size_t)k);
        }
        if (F.Tics == 1) {  /* if we want tics  */
            /* write out the minor and major tics below the reverse complement */
            tic_line = O_SEQ_LINE+2;
            /* if there's only 1 strand AND we want tics, have to move them up 1 line */
            if (F.Strands == 1) tic_line -= 1;
            for (i=0; i<(int)(basesPerLine/10); i++) memcpy (&O_txt[tic_line][i*10+O_LMARG],"    ^    *",10);
        }

        /* Following call to Translate has to calculate the # of bases at each call; otherwise end up with
        	overrun condition at end; "k" in code above
         This stanza doesn't have to change at all, even after the -T/-t change of 9.2.98 */

        rc = seq_len % 3;  /* rc is the correction factor for the trans, if not a perfect div by 3;
         if it is a perfect div by 3, then rc = 0, and there won't be any change */
        if (n_letters != 0) {
            for (mm=0; (mm<(int)Xn && mm<3); mm++) {
                if (k % 3 != 0) k = ((int) k /3)*3;   /* also added this code to Translate() for better modularity */
                /* this chunk doen't care about rc */
                seq_offset = ((block_Cnt-1)*basesPerLine+BASE_OVERLAP+1);
                Translate(sequence+seq_offset+mm, Prot_out, k, n_letters, Codons, codon_Table);
                /* need to place the translation frame labels here before memcpy() call */
                d = 0;
                if (F.Strands == 1) {
                    d += 1;
                }
                if (F.Tics == 0) {
                    d += 1;
                }
                sprintf(ctemp1, "%d", mm+1);
                O_txt[O_SEQ_LINE+3+mm-d][0] = ctemp1[0];
                memcpy(&O_txt[O_SEQ_LINE+3+mm-d][O_LMARG+mm], &Prot_out, (size_t)k); /* and copy it into place */
                lastline = O_SEQ_LINE + 3 + mm - d; /* keep track of last line of Xl */

                if (Xn == 6) { /* but if it's 6 frames of Xlation, need to do oppo at the same time */
                    /* but not EXACTLY the same - seq to be has to be shifted slightly to account for
                       those seqs not exactly / 3 - should have to only add a correction to the expr below
                       AND add the same correction to the offset for printing it out, so this chunk DOES care
                       about rc */

                    Anti_Par(sequence+seq_offset-mm+rc, s, k); /* get the *antipar* sequence */
                    /* Translate() below shouldn't have to change */
                    Translate (s, Prot_out, k, n_letters, Codons, codon_Table); /* Xlate it into N-letter code */
                    if (n_letters == 1) Reverse (Prot_out); /* to match the 'backwards' DNA strand */
                    else Triplet_Reverse (Prot_out);
                    /* need to place the translation frame labels here before memcpy() call */
                    sprintf(ctemp1, "%d", mm+4);
                    O_txt[O_SEQ_LINE+spacer+mm][0] = ctemp1[0];
                    memcpy(&O_txt[O_SEQ_LINE+spacer+mm][O_LMARG-mm+rc], &Prot_out, (size_t)k); /* and copy it into place */
                    lastline = O_SEQ_LINE + spacer + mm; /* this is the last line of the Xl, so last to print */
                }
            }
        } else {
            d = 0;
            if (F.Strands == 1) {
                d += 1;
            }
            if (F.Tics == 0) {
                d += 1;
            }
            lastline = O_SEQ_LINE + 2 - d;
        }

        /* This 'while' loop tests each text block for REs and also prints out those blocks that have no
        cuts in them at all but which have to be printed anyway....*/
        while ((DD->Dat[DCnt] != -22222) && (base_real_pos < ((block_Cnt * basesPerLine) + BASE_OVERLAP))) {
            DCnt++;
            /* if we're not at the end and the real pos'n <  current block + the overlap.. */
            /* !! Now DD->Dat[DCnt] should point to the corresponding  RE #  */

            /* if following test is true (they're !=) then it's safe to write it; */
            /* if it's false, then the same proto has already written to the same space */
            if (BackCheck[labs(DD->Dat[DCnt])] != labs(DD->Dat[DCnt-1])) { /* if !=, do this, but 1st... */
                BackCheck[labs(DD->Dat[DCnt])] = labs(DD->Dat[DCnt-1]);  /*  ... set it = so that it won't do this one again */
                okline = O_SEQ_LINE-2;

                /* This is where the RE filter could be inserted so that the DD data is ignored unless the RE that's
                  referenced matches the parameters set in the filters (can take this straight from the ladder-map
                  functions.  */
                /* direct from GelLad..
                   these are the conditions that have to be passed in order to have the REs placed in the
                   linear map. */
                /*  DD->Dat[DCnt] is still referencing a possibly neg #, so have to take the abs() */
                Cur_RE = labs((int)DD->Dat[DCnt]);
                /*  the previous horrendous multi-line 'if' filter has been functionized to  RE_Filters_OK()
                   for those missing it in this version.. */
                if ((DamDcmOK(base_real_pos, Cur_RE) == 1)   &&  RE_Filters_OK(Cur_RE)) { /* if olap = 0 && overhang = blunt */
                    /* this is the overall test for whether the RE will be printed out in the Linear Map or not
                     		if the RE DOESN'T meet the crit, it just gets ignored here, but will be printed in other outputs */

                    block_cut_pos =  (int)base_real_pos - 1 - ((block_Cnt-1) * basesPerLine); /* see if it's this easy?! */
                    /* Now locate it in the block, checking for under- and over-runs */
                    if ((block_cut_pos >= 0)  && (base_real_pos <= seq_len)) { /* if the position is greater than the beginning of the block */
                        if (block_cut_pos <= basesPerLine) { /* and if it's less than the end of the block... */
                            abcp = O_LMARG + block_cut_pos - 1; /* abcp = adjusted block_cut_pos - used x times below */
                            O_txt[O_SEQ_LINE-1][abcp] = '\\';   /* write the 'cut' indicator */
                            hits = 1;
                            /* locate the position for writing the name */
                            while (block_cut_pos < ok2wr[okline]) okline--; /* then go up to the lowest 'free' line */
                            if (okline != O_SEQ_LINE-2)  { /* but if can't place the name on the lowest line.. */
                                i = O_SEQ_LINE-1;
                                k = 0;   /* check if there's any space lower between previous names */
                                /* following 'while' must be ordered like it is - i will be decr only if k == 0 */
                                while (k == 0 && --i != okline) { /* leaving 1 space before and after the name */
                                    if (O_txt[i][abcp-1] == ' ' &&                /* this can be replaced with something */
                                            O_txt[i][abcp+RE[Cur_RE].E_nam_l] == ' ' &&  /* more efficient later */
                                            O_txt[i][abcp+3] == ' ' &&       /* must check for a space here... */
                                            O_txt[i][abcp+4] == ' ' )  k = 1;  /* as well as here */
                                    /* i = vertical counter ~= okline, k = 1 when we find enuf space */
                                }
                                okline = i;
                            }
                            if (okline < 2) {
                                fprintf(stderr, "!! Output Buffer Blown around base %ld - degeneracy caused too many Enz's\n"
                                        "to be matched!! Increase the stringency of Enz selection w/ -n -o -M flags\n", base_real_pos);
                                exit(1);
                            }
                            /* and memcpy the RE name from the struct to the output array */
                            memcpy (&O_txt[okline][abcp],&RE[Cur_RE].E_nam, (size_t)RE[Cur_RE].E_nam_l);
                            if (DD->Dat[DCnt] < 0) { /* if the match is in the (-) orientation */
                                memcpy (&O_txt[okline][abcp+RE[Cur_RE].E_nam_l], &os, 1);
                            }

                            /* and incr the ok2wr posn for that line by the length of the RE name + 1 */
                            /* 'ok2wr[okline]' below is modified only if we didn't find any lower spaces */
                            if (block_cut_pos+RE[Cur_RE].E_nam_l+1 > ok2wr[okline]) {
                                ok2wr[okline] = block_cut_pos+RE[Cur_RE].E_nam_l + 1;
                            }
                        } /* if (block_cut_pos < basesPerLine).. */
                        else {  /* otherwise it's in the overlap area from the next block */
                            if (in_OL == -1) {   /* if this is the 1st time in OVERLAP, mark it */
                                in_OL = 1;        /* by setting in_OL to 1  */
                                Cnt_reset = DCnt - 1;  /* Have to back up to 1st RE in the OVERLAP that did not
			                                                resolve into the current block */
                            }   /* The above marker *should* back the DCnt to the position of the Cur_RE */
                        }
                    }     /* if (block_cut_pos >= 0) ... */
                }
                base_real_pos = DD->Dat[++DCnt]; /* make sure that this is set before the next loop */
                if (okline < min_okline) min_okline = okline;
            } else {
                base_real_pos = DD->Dat[++DCnt]; /* this gonna fix it? */
            }
        }  /*  while (DD->Dat..... */

        /* Now print out the block we've composed, from min_okline to MAX_OUTPUT_LINES */
        /* adj to min_okline is nec to provide a 1 blank line between stanzas  */
        for (i=min_okline-hits; i<= lastline; i++) fprintf(stdout, "%.*s\n",p_width,O_txt[i]);

        /* And clean out the lines written into O_txt as well */
        memset(O_txt,' ',(size_t)((O_LMARG + O_RMARG + MAX_BASES_PER_LINE) * MAX_OUTPUT_LINES));  /* set all of O_txt to blanks */
    }   /* for (block_Cnt=1; block_Cnt<block_repeat; block_Cnt++) ... */
    free(BackCheck); /* this line was also deleted from the bike fork of tacg.c at one point
                       along with the corresponding memset at ~ line 36 */
    free(s);
}





/* GenerateTOC() simply checks the status of all the major output headers and
	generates the appro HTML to make up a Table of Contents for the results that
   follow.  It does not actually check that they have completed correctly, but
   it does keep everything in one place.  Just a series of if()'s  */

void GenerateTOC(char *datestamp)
{

    /* print a header */
    fprintf(stdout, "<h2>Table of Results</h2>\n");
    /* maybe print some date, version, file, machine run on, domain run on? etc info */
    fprintf(stdout, "<b>tacg version %s; Date: %s </b><br>\n"
            "<i>Contact Harry Mangalam (hjm@tacgi.com) with comments, critiques</i><br>\n",TACG_VERSION, datestamp);

    /* Start the UNordered list */
    fprintf(stdout, "<ul>\n");

    if (F.Summary != 0) { /* this generates no-cuts and all-cuts */
        fprintf(stdout, "<li><a href=\"#nocuts\">Patterns that DO NOT MAP to this sequence</a>\n");
        fprintf(stdout, "<li><a href=\"#allcuts\">Total Number of Hits per Pattern</a>\n");
    }
    if (F.Ladder != 0) { /* generates ladder map and infrequent matches */
        fprintf(stdout, "<li><a href=\"#ladmap\">Ladder Map of Pattern Matches</a>\n");
        fprintf(stdout, "<li><a href=\"#summarycuts\">Summary Map of Infrequent Matches</a>\n");
    }
    if (F.GelLO != 0 || F.GelHI != 0) {
        fprintf(stdout, "<li><a href=\"#gelmap\">Pseudo Gel Map</a>\n");
    }
    if (F.Sites != 0) {
        fprintf(stdout, "<li><a href=\"#sitetab\">Sites by Pattern</a>\n");
    }
    if (F.Frags == 1) {
        fprintf(stdout, "<li><a href=\"#fragunsort\">Fragment Sizes (UNSorted)</a>\n");
    }
    if (F.Frags == 2) {
        fprintf(stdout, "<li><a href=\"#fragsort\">Fragment Sizes (Sorted)</a>\n");
    }
    if (F.ORFs != -1) {
        fprintf(stdout, "<li><a href=\"#orfs\">Open Reading Frame Analysis</a>\n");
    }
    if (F.Orfmap == 1) {
        fprintf(stdout, "<li><a href=\"#orfmap\">ORF Map & MET/STOP Map</a>\n");
    }
    if (F.LinMap != 0) {
        fprintf(stdout, "<li><a href=\"#linmap\">Linear Map</a>\n");
    }
    if (F.GrfBBin != 0) {
        fprintf(stdout, "<li><a href=\"#graphics\">Graphics Plotting Data</a>\n");
    }
    if (F.Prox != 0) {
        fprintf(stdout, "<li><a href=\"#pair_proximity\">Pairwise Proximity</a>\n");
    }
    if (F.Rules != 0) {
        fprintf(stdout, "<li><a href=\"#rules_proximity\">Complex Rules-based Proximity</a>\n");
    }
    if (F.Hookey != -1) {
        fprintf(stdout, "<li><a href=\"#hookeyunsort\">AFLP Analysis - UNSORTED</a>\n");
        fprintf(stdout, "<li><a href=\"#hookeysort\">AFLP Analysis - SORTED</a>\n");
    }

    if (F.GrafMap != -1) {
        if (F.PDF == 1) {
            fprintf(stdout, "<li><a href=\"%s/tacg_Map.pdf \">PDF Plasmid Map</a>\n", F.tmppath);
        } else {
            fprintf(stdout, "<li><a href=\"%s/tacg_Map.ps \">Postscript Plasmid Map</a>\n", F.tmppath);
        }
    }
    /* spares below */
    /*
       if (F. != ) {
       	fprintf(stdout, "<li><a href=\"#LABEL\">DESCRIPTION</a>\n");
       }
     */
    /* below marks the end of the UNordered list and the beginning of the PRE-formatted output */
    fprintf(stdout, "</ul><HR noshade size=5><pre>\n");
    /*    fprintf(stdout, "<>  </>\n"); */

}



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
   Hmmm - this could be the prototype for an overall RE filter to be used in a number of cases..
*/
int DamDcmOK(long RWC, int REi)
{
    int i, OK=0;
    /* 1st, check to see if this RE is sensitive at all and if not, skip it */
    /* also skip if it doesn't have ANY cuts at all */
    if ((RE[REi].E_Ncuts > 0) && ((F.dam == 1 && RE[REi].dam != 100) || (F.dcm == 1 && RE[REi].dcm != 100))) {
        if (RE[REi].Sites[RE[REi].Sites[0]] == 0) {
            RE[REi].Sites[RE[REi].Sites[0]] = 1; /* 1st time thru, set this to 1 (point it at the 1st real value */
        }
        OK = 0;
        i = RE[REi].Sites[RE[REi].Sites[0]]; /* shorten: i = the last el to have been checked */
        while (RWC >= RE[REi].Sites[i] && OK == 0) {
            if (RWC == RE[REi].Sites[i]) {
                OK = 1;
            }
            i++;
        }
        RE[REi].Sites[RE[REi].Sites[0]] = i;
    } else OK = 1;
    if (RE[REi].E_Ncuts == 0) OK = 0; /* just in case */
    return OK;
}

/* Add_dam and Add_dcm just fill the RE[1] and RE[2] entries for these special
   methylation enzymes. */

void Add_dam(int NumREs)
{
    strcpy(RE[NumREs].E_nam,"dam\0");
    if ((RE[NumREs].E_raw_sit = calloc(5, sizeof(char))) == NULL) BadMem("Add_dam-1", 1);
    RE[NumREs].E_raw_sit = "GATC\0";
    if ((RE[NumREs].E_wsit[1] = calloc(5, sizeof(char))) == NULL) BadMem("Add_dam-2", 1);
    RE[NumREs].E_wsit[1] = "gatc\0";
    strcpy(RE[NumREs].E_hex[1],"gatcnn\0");
    RE[NumREs].E_tcut[1] = 0;
    RE[NumREs].E_olap = 0;
    RE[NumREs].E_len = 4;
    RE[NumREs].E_pal = 1;
    RE[NumREs].E_dgen = 16;
    RE[NumREs].E_mag = 4;
    RE[NumREs].proto = NumREs;
    RE[NumREs].Err = 0;
    RE[NumREs].E_nam_l = 3;
}

void Add_dcm(int NumREs)
{
    strcpy(RE[NumREs].E_nam,"dcm\0");
    if ((RE[NumREs].E_raw_sit = calloc(6, sizeof(char))) == NULL) BadMem("Add_dcm-1", 1);
    RE[NumREs].E_raw_sit = "CCWGG\0";
    if ((RE[NumREs].E_wsit[1] = calloc(6, sizeof(char))) == NULL) BadMem("Add_dcm-2", 1);
    RE[NumREs].E_wsit[1] = "ccwgg\0";
    strcpy(RE[NumREs].E_hex[1],"ccwggn\0");
    RE[NumREs].E_tcut[1] = 0;
    RE[NumREs].E_olap = 0;
    RE[NumREs].E_len = 5;
    RE[NumREs].E_pal = 1;
    RE[NumREs].E_dgen = 8;
    RE[NumREs].E_mag = 4;
    RE[NumREs].proto = NumREs;
    RE[NumREs].Err = 0;
    RE[NumREs].E_nam_l = 3;
}

/***************************************************************************************************/
/* The following code performs the "pad the beginning sequence with the end and the end
   with the beginning" regardless of whether it's the whole sequence or a subsequence,
   using the vars 'begin'  and 'end'.  *Bare_Seq_Len is the the exactly "Real Sequence
   Length" as determined from the SEQIO routines, no padding added or subtracted -
   it's incremented to TotPaddedSeqLen in here and passed back to calling fn() via
   *Pad_Seq_Cnt  */

/* TotPaddedSeqLen    length of sequence + bracketing overlaps
   Bare_Seq_Len       actual number of bases read
*/
/* will try to contruct this so that you give it the exact sequence you want padded, and it
   does all the calculations, work to alloc a new sequence space and do the padding, returning
   exactly what you want */

char *PadEndForEnd(char *SeqIn, long *Pad_Seq_Len, long *Bare_Seq_Len)
{
    long m, j, l, BareLen, TotPaddedSeqLen, FrPadSeqLen;
    /*   int  f1=0, iseq; */
    char *sequence, blank[BASE_OVERLAP+1];
    memset(blank, 'a', BASE_OVERLAP+1);

    BareLen = strlen(SeqIn);
    if (F.Verbose > 1) fprintf(stderr, "SEQIn = %s\n", SeqIn);
    /* know exactly how much mem to allocate for version 2 from SEQIO call */
    if ((sequence = (char *)calloc((BareLen+(2*BASE_OVERLAP)+3), sizeof(char))) == NULL) {
        BadMem("Error alloc'ing *sequence in PadEndForEnd", 1);
    }
    FrPadSeqLen = BareLen + BASE_OVERLAP +1;  /*  that '+1' may be off tho */
    sprintf(sequence+BASE_OVERLAP+1, "%s", SeqIn);
    if (BareLen < BASE_OVERLAP)  { /* if it's a v short seq, pad it with junk */
        memcpy(&sequence[0], &blank, BASE_OVERLAP+1);   /* pad the beginning with a's */
        memcpy(&sequence[FrPadSeqLen], &blank, BASE_OVERLAP+1);  /* pad the back with a's */
    } else { /* else do the std tango of copy the end to the begin; begin to the end */
        for (l=BASE_OVERLAP+1,m=0; l <= BASE_OVERLAP*2+2 && m <= BASE_OVERLAP; l++,m++)  {
            sequence[FrPadSeqLen+m] = sequence[l]; /* copies seq past the end from the the beginning */
        }  /* note that in the above 'for', 'FrPadSeqLen' itself doesn't change */

        for (j=BASE_OVERLAP+1,m=0; j >= 0 && m <= BASE_OVERLAP; j--,l++,m++)  {
            sequence[m] = sequence[FrPadSeqLen-j]; /* copies seq to before the beginning from the end */
        }  /* note that in the above 'for', 'FrPadSeqLen' itself doesn't change */
    }

    if ((F.SeqBeg != 1 || F.SeqEnd != 0) && F.Verbose >= 1) fprintf(stderr, "\nSubsequence is %ld bases.\n", BareLen);

    TotPaddedSeqLen = BASE_OVERLAP + FrPadSeqLen;
    sequence[TotPaddedSeqLen+1] = '\0';  /* and term the sequence string, then decr it to point to end of the seq  */
    /* and free the mem not used by the sequence string - shouldn't be needed anymore because it was
    	EXACTLY alloc'ed at the beginning, BUT, it SEQIO doesn't handle sequences the way I handle
    	sequences, it could still be somewhat different - ie SEQIO could pass in spaces, numbers,
    	newlines, etc and count them in the size of the array, so this should be left in until it can be
    	absolutely verified that it's correct.  */

    /* sequence is padded with extra bases to allow circular cutting, but anything that depends on DNA length
       will reference the specific variable 'BareLen' */

    /* now we have the sequence in one long string, with 'BASE_OVERLAP' bp buffers at each end, so cut it as
       overlapping hexamers, starting at sequence[0], but only record cuts in the real seq */
    *Pad_Seq_Len = TotPaddedSeqLen; /* send it back to main() ... a greenie weenie way to do it */
    *Bare_Seq_Len = BareLen;       /*     ditto     */
    if (F.Verbose > 1) fprintf(stderr, "sequence = %s\n", sequence);

    return sequence;  /*  return the address of sequence ) */
}   /* End of PadEndForEnd() */



/* ReverseTranslate() is another eponymous fn(), taking a string of single letter aminia acids, it
   munges them in to reverse translated DNA, initially towards contributing to '--silent' (the
   fn() that looks for SILENT SITES, RE sites that can be introduced by mutation that will not
   change the coding sequence.
   WARNING: In doing this, there are some translational ambiguities produced.  It depends on which
   Codon table is being used, but for teh Standard one, Arg, Leu, and Ser are reverse translated to
   mgn, ytn, and wsn respectively, which will not forward translate unambiguously, so it remains
   to the user to be careful about interpreting th eresults.  If I can come up with a plan to
   verify the sequence once this is done, I'll implement it. */

char *ReverseTranslate(char *Prot, char Codons[MAX_ORGS][N_CODONS][7], int OrgTable)
{

    short int Add8[4]= {0,8,0,8}, i, j, k, aa, sum, codon, tripbase, basei=0, basev=0, count,
                                                                     AAi; /* Amino Acid index */
    char TT[16]= {'a','c','g','v','t','h','d','b',' ','m','r','s','w','y','k','n'};
    /* ,'A','B','C','D','E','F','G' */
    char c, RevTransTab[26][4], *DegenDNA;
    short int AAB[26][3][4]; /* AAB = matrix of bases that make up an Amino Acid */
    int Plen = 0, Dlen = 0;

    /* set all AAB els to -1 */

    for (i=0; i<26; i++) {
        for (j=0; j<3; j++) {
            for (k=0; k<4; k++) {
                AAB[i][j][k] = -1;
            }
        }
    }

    Plen=strlen(Prot);
    Dlen=Plen*3;
    if ((DegenDNA = (char *)calloc(Dlen + 1, sizeof(char))) == NULL) {
        BadMem("no mem at ReverseTranslate:DegenDNA\n", 1);
    }

    for (codon=0; codon<64; codon++) {
        c = Codons[OrgTable][codon][0];
        AAi = ((int)c) - 65; /* 65 = decimal 'A' */
        if (AAi < 0) {  /* case for '*' */
            if (F.Verbose > 0) fprintf(stderr, "caught AAi < 0 (%c), correcting to 25.\n", Codons[OrgTable][codon][0]);
            AAi = 25;
        }
        for (tripbase=0; tripbase<3; tripbase++) {
            /* using the trailing part of the Codons label = "TThract" (thats why 'tripbase+4' */
            switch (c = Codons[OrgTable][codon][tripbase+4]) {
            case 'a':
                basev =    basei = 0;
                break;
            case 'c':
                basev =    basei = 1;
                break;
            case 'g':
                basev =    basei = 2;
                break;
            case 't':
                basev = 4;
                basei = 3;
                break;
            default:
                fprintf(stderr, "Caught a weird char (%c) in ReverseTranslate().\n", c);
                break;
            }
            AAB[AAi][tripbase][basei] = basev;
        }
    }

    /* now AAB is loaded with the matrix, so go thru AAB and calculate the RT IUPAC values (from TT) */
    for (aa=0; aa<26; aa++) {
        for (i=0; i<3; i++) {
            count = -1;
            sum = 0;
            for (j=0; j<4; j++) {
                if (AAB[aa][i][j] >= 0) {
                    sum += AAB[aa][i][j];
                    count++;
                }
                RevTransTab[aa][i] = TT[Add8[count] + sum];
            }
            RevTransTab[aa][3] = '\0';  /* terminate it so it can act like a string */
        }
        if (F.Verbose > 1) fprintf(stderr, "IUPAC code for %c = %s.\n", ((char)(aa+65)), RevTransTab[aa]);
    }
    /* at this point, the translation array is set up and ready to access, so Finally ready
       to start to reverse translate!!   */

    for (aa=0; aa<Plen; aa++) {
        AAi = ((int)Prot[aa]) - 65;
        if (AAi < 0) AAi = 25;
        for (i=0; i<3; i++) {
            DegenDNA[(aa*3)+i] = RevTransTab[AAi][i];
        }
    }
    DegenDNA[Dlen] = '\0'; /* NULL terminate like a good boy */
    if (F.Verbose > 1) fprintf(stderr, "Prot = %s\nDegen DNA = %s\n", Prot, DegenDNA);
    return DegenDNA; /* and return it */
}




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

int BitArray(char *bstr, int I, int mode)
{
    unsigned int BytOff, bitOff;
    unsigned char mask[8] = {  0200,   /* 128 in octal */
                               0100,   /*  64 in octal */
                               040,   /*  32 in octal */
                               020,   /*  16 in octal */
                               010,   /*   8 in octal */
                               04,   /*   4 in octal */
                               02,   /*   2 in octal */
                               01    /*   1 in octal */
                            };
    char ch;   /* cannot be ** unsigned ** as that allows an additional 128 dec chars to it */
    BytOff = (int)(I/8); /* calc the # of bytes/chars to offset */
    bitOff = I - (BytOff*8);  /* # of bits that have to be added to the Byte Offset */
    /* mask is the array that for each additional bit offset contains the mask that
       has to be applied to the byte to coerce it into the right combo of bits. */
    /* 1st, some crude error checking...BUT this doesn;t check the bounds of the bstr handed to it -
       assumes that the bounds are oK as is .. for now */
    if (mode < 0 || mode > 2) {
        fprintf(stderr, "BitArray(): mode is outside of range 0-2!\n");
        exit(1);
    }

    if (mode == 0) { /* if asking for a clear operation */
        ch = ~bstr[BytOff] | mask[bitOff];  /* this is extremely convoluted, but it works */
        bstr[BytOff] = ~ch;
        return 1;
    } else if (mode == 1) { /* if asking for a set operation */

        bstr[BytOff] = bstr[BytOff] | mask[bitOff]; /*  it's that ole '|' (inc OR) again.. */
        return 1;
    } else {
        ch = (bstr[BytOff] | mask[bitOff]); /* so after all that, it's '|' (incl OR) that masks it correctly */
        if (ch != bstr[BytOff]) { /* then must be asking for a query operation (mode == 2) */
            /* then the queried bit is NOT set (it doesn't match the bit set in mask[] */
            return 0;
        } else {  /* the queried bit IS the same as the mask ( the ^ op returns 0) */
            return 1;
        }
    }
}


void Read_NCBI_Codon_Data(char Codons[MAX_ORGS][N_CODONS][7], char *Codon_Labels[MAX_ORGS])
{
    /* Codons - array that holds all the codon preference data for the different
                organisms or mitos loaded from external file "codon.data"
       Codon_Labels - array that holds the organism labels (Universal, various Mitos, etc. */
    FILE *fpCodon;
    int e, i, j, k,  sum[256];
    char triplet[3], s[256], labels[4], *ptmp, ctmp[512], L[10];
    for (i=0; i<3; i++)triplet[i] = 'a';

    if ((ptmp = (SearchPaths("codon.data", "CODON PREFS"))) == NULL) {
        fprintf(stderr,"Can't open the 'codon.data' file!! - Check spelling, ownership, existance, etc - Bye!!\n");
        exit (1);
    }

    /*  open the Codon data file  */
    if ((fpCodon=fopen(ptmp,"r"))  ==  NULL)   {   /*  if  it's not there and readable */
        fprintf(stderr,"Cannot open the file \"codon.data\" for reading (in function  Read_NCBI_Codon_Data())!!\n");
        exit(1); /* print an error and die gracefully */
    }
    e=i=0;
    while (((fgets(s, 256, fpCodon)) != NULL) && (i < MAX_ORGS)) {
        if (s[0] != '#' && s[0] != '\n') {
            if (F.Verbose >2) {
                fprintf(stderr, "Accepted in 'codons.data': %s\n", s);
            }
            sscanf(s,"%s %s", L, ctmp); /* stuff the stanza header string into the label array */
            Codon_Labels[i] = (char *)calloc(strlen(ctmp)+1, sizeof(char));
//			fprintf(stderr, "%s\n", ctmp);
            if (Codon_Labels[i] == NULL) BadMem("No mem for Codon_Labels[i]\n", 0);
            strcpy (Codon_Labels[i],ctmp);
            for (j=0; j<64; j++) { /* for each doublet of values in the stanza (64 in all) */
                /* The next 2 lines !work with gcc on Linux if the -O2 flag is used.  if no -O2 flag, it's OK ?!?!? */

                e = fscanf(fpCodon,"%3s %c %3s",ctmp, labels, labels+1); /* grab the triplet and the corresp. label */
                if (e <= 0) {
                    fprintf(stderr, "SeqFuncs.c:1963 Failure to read\n");
                    exit(1);
                }
//					fprintf(stderr, "%d: %s %s  ", j, ctmp, labels);
                for (k=0; k<3; k++) triplet[k] = tolower(ctmp[k]); /* should change 'tolower' to a macro */
                hash(triplet,sum,3); /* hash the triplet to use as an index as to where to put the label */

                memcpy(&Codons[i][sum[0]][0],&labels,4); /* and copy the label to the right place in "Codons" */
                /* below should copy the actual triplet to the end of the 'TThr' label, so it can be picked up
                   in ReverseTranslate() */
                memcpy(&Codons[i][sum[0]][0]+4,&triplet,3); /* after the 'TThr' stuff, add the actual bases  */
            }
            i++;
//				fprintf(stderr, "\n");
        } else if (F.Verbose >2) {
            fprintf(stderr, "Ignored in 'codons.data': %s\n", s);
        }
    }
    free(ptmp);
}



/* ExtractSeqAtHit() takes the list of patterns passed to the program via the '-X' flag and extracts
   the sequences around the hits beginning 'b' bases 5' of the pattern midpoint and ending 'e' bases
   3' of the pattern midpoint, where 'b' and 'e' are the parameters of the --extract b,e,revcomp flag.
   ('revcomp' is the variable that indicates whether the sequence should be revcomplemented before
   being sent to stdout).  the only thing that it has to be passed is the number of REs to chew thru */
void ExtractSeqAtHit(SEQINFO *SEQinfo, int seq_len, char *sequence, int ProtoCutters, int *Protos)
{
    int eRE, PC, site, rc=0, clen, otherstrand;
    long rbos=0, reos=0, APO=0;
    char *xtract;
    char *RCLabel[] = {"->",   /*  [0]  */
                       "<-"
                      };  /*  [1]  */

    clen = (int)(1 + F.XtEnd + labs(F.XtBeg)); /* '1+' for the starting base */
    if ((xtract = (char *) calloc(clen+1, sizeof(char))) == NULL) BadMem("Can't init mem for clen", 1);

    if (F.XtBeg < 0) rc = 1; /* if [50]<0, want revcomp seq if pattern is on other strand */

    /* set up the output*/
    fprintf(stdout, "\n\n== Extracted Sequences Surrounding Hits\n");
    // MUST check for requests that try to extract sequences beyond actual sequence endpoints.
    // have to check at each request that borders of extract are within sequence
    // warn that the extract will not currently correctly handle circular extracts
    // (as tacg doesn't have a routine for backtracking seuence all the way -
    // this should be straightforward, but haven't done it)
    for (eRE = 0; eRE < ProtoCutters; eRE++) {
        PC = Protos[eRE];
        for (site = 1; site <= RE[PC].E_Ncuts; site++) {
            otherstrand = 0;
            if (BitArray(RE[PC].Orient, site-1, 2) == 1) otherstrand = 1;

            if (otherstrand == 0 ||rc == 0) { /* if F strand, (or don't want the RC, calc the seq to include */
                rbos = RE[PC].Sites[site] - labs(F.XtBeg) + 1; // works for -p patterns (w/wo errors)
                reos = RE[PC].Sites[site] + F.XtEnd ;
                if (otherstrand == 1) APO=1; /* and twiddle the offset to compensate for the other strand */
            } else if (rc == 1) { /* if R strand, swap the offsets if want the RC */
                rbos = RE[PC].Sites[site] - F.XtEnd + 1; // shortens the length of the xt seq
                reos = RE[PC].Sites[site] + labs(F.XtBeg + 2) ;
            }

            if (rbos < 1) {
                rbos = 1;
            }
            if (reos > seq_len) {
                reos = seq_len;
            }
            clen = reos - rbos; //  + 1?
            /* if that site was in the revcomp && we want to get the rc in that case... */
            if ((otherstrand == 1) && rc == 1) {
                /* cuz we want the RC but as the top strand, it's actually AntiParallel */
                /* BUT!! if it's AP, have to play with the coordinates so that we're extracting the same
                   coordinates relative to the F strand */
                // Here's where we have to check if the borders of the requested extract lie within the
                // actual sequence and truncate if they don't - cannot submit a wonky sequence to Anti_Par()


                Anti_Par(sequence+rbos+1+BASE_OVERLAP-2, xtract, clen);
            } else {
                strncpy(xtract, sequence+rbos-APO+BASE_OVERLAP, clen);
            }
            xtract[clen] = '\0'; /* make sure it's 0-terminated */
            if (SEQinfo != NULL) {

                if (SEQinfo->description == NULL) {
                    SEQinfo->description = "None";
                }
                //else { 	strncpy(tmp_seq_id, SEQinfo->description, 23 ); }
                if (SEQinfo->dbname == NULL)      SEQinfo->dbname = "None";
                //TODO break description into 1st word as ID, use that as naming component.
                fprintf(stdout, ">%s_%ld_%.20s:   Seq from %ld to %ld (%s), Descr: %s, DB: %s, \n%s\n",
                        RE[PC].E_nam, RE[PC].Sites[site], SEQinfo->description, rbos, reos, RCLabel[otherstrand], SEQinfo->description, SEQinfo->dbname, xtract);
            } else {
                fprintf(stdout, ">%s_%ld:   Seq from %ld to %ld (%s), Descr: Raw Sequence, DB: N/A\n%s\n",
                        RE[PC].E_nam, RE[PC].Sites[site], rbos, reos, RCLabel[otherstrand], xtract);
            }

        }
    }
    fprintf(stdout, "\n");
    free(xtract);
}


/* Hookey() allows the use of 'tagging' a single RE in a combination of REs (typically a tagged
   6 or 8 cutter, paired with an untagged 4 cutter) so that you can pull out those fragments that have
   only an end cut by the tagged RE.  John Hookey suggested this as an assist for some genome
   mapping work he was doing.  Hookey now prints out both the UNSORTED and the SORTED list of fragments */
void Hookey(int H, int AI, long seq_len)
{
    /* we got here if we've detected a '=' appended to one of the REs in the '-x' stream */
    /* where H = NumREs (post-incremented in the function call) and AI = 'All Included', the index
       that points to the RE entry that keeps track of ALL the RE hit */
    char t;
    int i, j=1, k=1, l=0, m=0, reps=0, HTML;
    long LastTagged=-1, EoAI, *ToSort;

    HTML = F.HTML; /* barely nec, but ... */
    if ((ToSort = (long *)calloc((3 * (RE[F.Hookey].E_Ncuts + 2)), sizeof(long))) == NULL) {
        BadMem("in Hookey(), calloc failed: ToSort\n", 1);
    }
    /* the RE[=] is ID'ed by F.Hookey and SelEnz array: SelEnz[F.Hookey] */
    RE[H].Sites = (long *) calloc((3 * RE[F.Hookey].E_Ncuts) + 2, sizeof(long)); /* '+ 2' for extra frag */
    if (RE[H].Sites == NULL) BadMem("in Hookey(), calloc failed: RE[H].Sites\n", 1);
    RE[H].Frags = (long *) calloc((2 * RE[F.Hookey].E_Ncuts) + 2, sizeof(long)); /* '+ 2' for extra frag */
    if (RE[H].Frags == NULL) BadMem("in Hookey(), calloc failed: RE[H].Frags\n", 1);

    if (RE[AI].Sites[1] == RE[F.Hookey].Sites[1]) {   /* check beginning condition */
        RE[H].Sites[k++] = RE[AI].Sites[1];    /* copy the site */
        RE[H].Sites[k++] = RE[AI].Sites[2];    /* and the site after to the RE[H] entry */
        LastTagged = RE[AI].Sites[2];
        /* and need to generate the frags associated with these sites as well */
        RE[H].Frags[l++] = RE[AI].Sites[1];     /* 1st frag = the 1st site ofcombined */
        RE[H].Frags[l++] = RE[AI].Sites[2] - RE[AI].Sites[1]; /*  2nd frag = 2nd site - 1st site */
        j++; /* and k = 3, j = 2, l = 2 coming out of this stanza;  */
    }
    /* print the header for the Hookey Output - done here exactly as wanted */
    /* obviously still needs the HTML header stuff for the web version */
    if (HTML > -1) {
        fprintf(stdout, "<A NAME=\"hookeyunsort\"><br><h3>");
    }
    fprintf(stdout, "\n== Hookey Output **UNSORTED** ('=' indicates tagged site):");
    if (HTML > -1) {
        fprintf(stdout, "</h3>");
    }
    fprintf(stdout, "\n  Fragment Size    Begin Site    End Site\n");

    EoAI = RE[AI].Sites[0];  /* End of AI, used below */
    for (i=2; i<EoAI; i++) { /* traverse all the combined cuts, starting from 2... */
        if (RE[AI].Sites[i] == RE[F.Hookey].Sites[j] && (j < RE[F.Hookey].Sites[0])) {
            /* if the combined sites == one of the '=' sites */
            /* but note that these have to be handled differently than a regular .Sites[] array because
               there are 3 sites per 2 frags instead of the usual 1:1 correspondence.  Actually don't have to
               deal with it at all for Hookey, but if you wanted to extract the sites, would have to present them
               a bit differently */
            RE[H].Sites[k++] = RE[AI].Sites[i-1];  /* copy the site before */
            RE[H].Sites[k++] = RE[AI].Sites[i];    /* the site itself */
            RE[H].Sites[k++] = RE[AI].Sites[i+1];  /* and the site after to the RE[H] entry */

            /* and generate the frags corresp'ing to these sites.. */
            /* next test prevents successive sites from the same RE from generating additional frags */
            /* this should be re-org'ed after it gets working to move the 'if's  so that they're
                checked only once per run */
            if (RE[AI].Sites[i-1] != LastTagged) {
                RE[H].Frags[l++] = RE[AI].Sites[i] - RE[AI].Sites[i-1];
                fprintf(stdout, "%10ld      %10ld  %10ld=\n", RE[H].Frags[l-1],
                        RE[AI].Sites[i-1], RE[AI].Sites[i]);
                ToSort[m++] = RE[H].Frags[l-1];
                ToSort[m++] = RE[AI].Sites[i-1];
                ToSort[m++] = RE[AI].Sites[i] * -1;
            }
            if ((i + 1) < RE[AI].Sites[0]) {
                RE[H].Frags[l++] = RE[AI].Sites[i+1] - RE[AI].Sites[i];
                if (RE[AI].Sites[i+1] == RE[F.Hookey].Sites[j+1]) {
                    t = '=';
                } else {
                    t = ' ';
                }
                fprintf(stdout, "%10ld      %10ld= %10ld%c \n", RE[H].Frags[l-1],
                        RE[AI].Sites[i], RE[AI].Sites[i+1], t);
                ToSort[m++] = RE[H].Frags[l-1];
                ToSort[m++] = RE[AI].Sites[i] * -1;
                ToSort[m++] = RE[AI].Sites[i+1];
            }
            LastTagged = RE[F.Hookey].Sites[j]; /* for clarity */
            j++;
        }
        /* test for ending condition - if the last of the RE[=] sites == the last of the RE[AI] sites .. */
        if (RE[F.Hookey].Sites[j] == RE[AI].Sites[RE[AI].Sites[0] - 1] ) {
            RE[H].Frags[l++] = seq_len - RE[F.Hookey].Sites[j]; /* the last frag = sequence length - the last site */
        }
    }
    /* now qsort ToSort by 3*long and re-print to get SORTED output */
    reps = (m)/3;
    qsort(ToSort, (size_t)reps, 3*sizeof(long), compare);
    if (HTML > -1) {
        fprintf(stdout, "<A NAME=\"hookeysort\"><br><h3>");
    }
    fprintf(stdout, "\n\n== Hookey Output **SORTED** ('=' indicates tagged site):");
    if (HTML > -1) {
        fprintf(stdout, "</h3>");
    }
    fprintf(stdout, "  Fragment Size    Begin Site    End Site\n");

    for (i=0; i<m; i += 3) {
        if (ToSort[i+2] < 0) {
            fprintf(stdout, "%10ld      %10ld  %10ld=\n", ToSort[i], ToSort[i+1], ToSort[i+2] * -1);
        } else {  /*                               ^                                       */
            fprintf(stdout, "%10ld      %10ld=  %10ld\n", ToSort[i], ToSort[i+1] * -1, ToSort[i+2]);
        }  /*                               ^                                              */
    }

    /* fill in the RE descriptive headers */
    sprintf(RE[H].E_nam, "Hookey");
    RE[H].Sites[0]    = k;     /* just to make it like all the others [0] points to the next free space */
    RE[H].E_Nfrags    = l-1;   /* Need this because Sites and Frags are asymmetric at [0]  */
    RE[H].E_Ncuts     = k-1;   /* and adjust the number */
    RE[H].proto       = H;     /* and make sure that it thinks it a proto */
    RE[H].E_mag       = 8;     /* and make sure the mag is enough to let it pass the default checks */
    RE[H].E_pal       = 1;     /* and make sure that it thinks it's a pal */
    RE[H].E_raw_sit   = "Tagged Site";   /* and give it a label so that it can output something */

    /* and spit out a coupla extra lines to separate it from the next stanza */
    /* fprintf(stdout, "\n\n"); */
} /* end Hookey() */



/* DumpDataForPlot() simply dumps data that is stored in internal arrays to stdout so
    that it can be manipulated by external plotting and analytical apps - takes care of all
    output as well as the array manipulations that are required to flip the internal arrays
    on the diagonal or make them long and skinny (as gnuplot prefers)   */

void DumpDataForPlot(long seq_len,  int ProtoCutters,  int *Protos)
{
    int Nbins, **Grfdata, l, j, W, i, HTML;
    HTML = F.HTML;
    W = (int)F.GrfFmt;  /* F.GrfFmt indicates the type of output wanted -
                                    1 = 'square' output (Bins X, REs Y)
                                    2 = 'square' output (Bins Y, REs X),
                                    3 = long output (2 cols, reiterated Bins (X), RE data(Y) */
    Nbins = (int) (seq_len/F.GrfBBin) + 1; /* Nbins = # bins required */
    Grfdata = (int **) calloc(Nbins, sizeof(int*)); /* init the rows of pointers for Grfdata */
    if (Grfdata == NULL) BadMem("Can't init mem for Grfdata", 1);
    for (i=0; i<Nbins; i++) { /* and all the columns */
        Grfdata[i] = (int *) calloc(ProtoCutters+1, sizeof(int));
        if (Grfdata[i] == NULL) BadMem("SeqFuncs.c: Can't init mem for Grfdata[i]", 1);
    }

    fprintf(stdout, "\n\n\n"); /* prefix to GraphData title */
    if (HTML==1) fprintf(stdout, "<H3><A NAME=\"graphdata\">");
    fprintf(stdout, " == Graphics Data\n\n"); /* Print the header title */

    l = F.GrfBBin;
    for (i=0; i<Nbins; i++) { /* loads the Bins #s for Grfdata */
        Grfdata[i][0] = l;
        l += F.GrfBBin;
    }
    /* now load the RE data part */
    for (i=0; i<ProtoCutters; i++) { /* for each of those REs that have more than 1 hit */
        for (j=1; j<RE[Protos[i]].Sites[0]; j++) { /* for each hit */
            Grfdata[(int)(RE[Protos[i]].Sites[j]/F.GrfBBin)][i+1]++; /* base hit / (#bases/bin) */
        }
    } /* on exit, data is loaded correctly */
    /* now the diff output forms, massaging the square array Grfdata into the right shape */
    if (W == 1) { /* square, Bins X, REs Y, (original format) */
        /* 1st print the Bins out */
        fprintf(stdout, "      Bins");
        for (i=0; i<Nbins; i++)  fprintf(stdout, "%8d ", Grfdata[i][0]);
        fprintf(stdout, "\n");
        for (i=1; i<ProtoCutters+1; i++) {
            fprintf(stdout, "%*s", F.MaxPatNameLen, RE[Protos[i-1]].E_nam); /*print RE label */
            for (j=0; j<Nbins; j++) fprintf(stdout, "%8d ", Grfdata[j][i]); /* and the # */
            fprintf(stdout, "\n");
        }
    } else if (W == 2) { /* square, Bins Y, REs X, (diagonal flip of original format) */
        printf("# Reminder: how to plot using gnuplot:\n#edit file to remove all headers down to the lines prefixed by '#'\n# start gnuplot\n# $ gnuplot\n# gnuplot> set key autotitle columnheader\n# (columnheader' uses the header column for key labels.)\n# gnuplot> plot 'filename' using 1:2\n# gnuplot> replot 'filename' using 1:3\n# gnuplot> replot 'filename' using 1:4\n# etc for each pattern...\n# for lines connecting the dots, use:\n# gnuplot> plot 'filename' using 1:2 with lines\n\n");

        fprintf(stdout, "      Bins"); /* 1st print out the header line with 'Bins' */
        for (i=0; i<ProtoCutters; i++) { /* and the RE labels */
            fprintf(stdout, "%*s", F.MaxPatNameLen, RE[Protos[i]].E_nam); /*print RE label */
        }
        fprintf(stdout, "\n");
        for (i=0; i<Nbins; i++) {
            for (j=0; j<ProtoCutters+1; j++) fprintf(stdout, "%10d", Grfdata[i][j]); /* and the # */
            fprintf(stdout, "\n");
        }
    } else { /* Long, skinny form, reiterated Bins(X) RE(Y) for each RE */
        for (i=1; i<ProtoCutters+1; i++) {
            fprintf(stdout, "\n#     Bins  %*s\n", F.MaxPatNameLen, RE[Protos[i-1]].E_nam);
            for (j=0; j<Nbins; j++) {
                fprintf(stdout, "%10d %10d\n", Grfdata[j][0], Grfdata[j][i]);
            }
        }
    }
    fprintf(stdout,"\n");
    free(Grfdata);
}   /* End of DumpDataForPlot() */



/* this horrible bit of rubbish is another of the ilk of HorribleAccounting() - nasty, but common code that
   is being sent off to the penal colony of a separate function so it won't contaminate the morals and
   elegance of the rest of ReadEnz() - returns the 'magnitude' of the string given it */
float MagCalc(char *s, int len)
{
    float fmag = 0; /* float counter for magnitude */
    int i;
    for (i=0; i<len; i++) { /*for each base in the site ...*/
        switch (s[i])  {
        default:
            fprintf (stderr, "Strange character (%c) in MagCalc() (in pattern %s).\n", s[i], s);
            break;
        case 'a':
        case 'c':
        case 'g':
        case 't':
            fmag += 1;
            break;
        case 'y':
        case 'r':
        case 'm':
        case 'k':
        case 'w':
        case 's':
            fmag += 0.5;
            break;
        case 'b':
        case 'd':
        case 'h':
        case 'v':
            fmag += 0.25;
            break;
        case 'n':
            break; /* 'n' doesn't cause the mag to be increased at all */
        }
    }
    return fmag;
}


/*
 * Re-installation of the original GetSequence() to allow reading of raw seq -
 * useful in a variety  of situations, and for compatibility with the version 2
 * input model.
 */

/* GetSequence reads and formats the sequence from stdin and returns the length of the sequence, as well as a pointer
   to the array that holds it (*sequence) along with the bracketing repeats - needs only a few extra variables to do it
   and it cleans up main considerably. begin and end are indices to the real world coordinates of the sequence */
/* tot_seq_Cnt    length of sequence + bracketing overlaps
   seq_length     actual number of bases read
*/
char *GetSequence(long *tot_seq_Cnt, long *seq_length, long *Degen_Log[])
{
    long eoseq = 30000, m, j, l, seq_len, seq_Cntr=1, begin = 1, end = 1000000000, totSeqCnt;
    int  indegen=0, DLi=2, /* DLi = 2, as [0] = alloc'ed size, [1] is next el to write  */
         c, realloc_iter = 0, f1=0;
    char *sequence, ct, blank[BASE_OVERLAP+1], *infile;
    memset(blank, 'a', BASE_OVERLAP+1);

    FILE *INFILE;
    if ((infile = (char *) calloc (256, sizeof(char))) == NULL) BadMem("initial calloc for infile",1)  ;
    /* Get mem for sequence */
    sequence = (char *) calloc (30000, sizeof(char));  /* init the pointer to 30,000 */
    if (sequence == NULL) BadMem("initial calloc for sequence", 1);

    /* load the sequence array, filtering as we go */
    totSeqCnt = BASE_OVERLAP + 1; /* to allow the wrapped seq to be inserted before the real seq
                                    (and to allow real seq to start at x1) */
    begin = F.SeqBeg;  /* set in SetFlags() */
    if (F.SeqEnd != 0) end = F.SeqEnd;    /* Ditto */

    if (F.infile == NULL) {
        INFILE = stdin;
    } else {
        if ((infile = SearchPaths(F.infile,"Seq Input File")) == NULL) {
            fprintf(stderr, "Can't find this file: %s - please check it and try again.\n", F.infile);
            exit(0);
        }
        INFILE = fopen(infile,"r");
    }

    while ((c = getc(INFILE)) != EOF) { /* new - stdin-oriented;  */
        if (eoseq - totSeqCnt <10) {
            sequence = (char *)realloc(sequence, sizeof(char)*(eoseq+30000));
            if (sequence == NULL) BadMem("realloc sequence", 1); /* this does not pass realloc_iter... */
            eoseq = eoseq + 30000;
            realloc_iter++;
        }
        /* Now real seq will start at 'BASE_OVERLAP'; before that is wrapped seq from the end and at the end is
           wrapped seq from the beginning - required to handle circular seqs */
        ct = tolower(c);
        if (ct == 'u') ct = 't'; /* converts u's to t's */
        /* no longer need if/else decision for different modes - slightly slower this way but more compact and better
           user interaction ...? */
        switch (ct) {
            /* this is where the originalbase-handling chunk went... */
        case 'a': /* favored 4 always get added */
            F.NBases[0]++;
            goto NowProcess;
        case 'c':
            F.NBases[1]++;
            goto NowProcess;
        case 'g':
            F.NBases[2]++;
            goto NowProcess;
        case 't':
            F.NBases[3]++;

NowProcess:
            /* if its acgt && there was a degen beforehand, have to mark the end of that */
            /* this shouldn't need to be checked for mem as it just adds 1 value to the end */
            if (indegen == 1) {
                (*Degen_Log)[DLi++] = seq_Cntr-2;
                (*Degen_Log)[1]++;
                indegen = 0;
            }

        case 'y':
            F.NBases[4]++;
            goto NowProcessDegens; /* but if it's a valid IUPAC code */
        case 'r':
            F.NBases[5]++;
            goto NowProcessDegens;
        case 'w':
            F.NBases[6]++;
            goto NowProcessDegens;
        case 's':
            F.NBases[7]++;
            goto NowProcessDegens;
        case 'k':
            F.NBases[8]++;
            goto NowProcessDegens;
        case 'm':
            F.NBases[9]++;
            goto NowProcessDegens;
        case 'b':
            F.NBases[10]++;
            goto NowProcessDegens;
        case 'd':
            F.NBases[11]++;
            goto NowProcessDegens;
        case 'h':
            F.NBases[12]++;
            goto NowProcessDegens;
        case 'v':
            F.NBases[13]++;
            goto NowProcessDegens;
        case 'n':
            F.NBases[14]++;
            goto NowProcessDegens;

NowProcessDegens:

            /* regardless of F.Degen, may want to log degens, so .. */
            if (F.LogDegens) {
                if (indegen == 0) {  /* if 0, 1st degen after a non-degen, so needs to be marked */
                    /* (*Degen_Log)[even] = start, [odd] = end for all degens; doesn't distinguish among dif degens */
                    /* always check for end of alloc'ed mem */
                    if ((*Degen_Log)[0] - (*Degen_Log)[1] < 10) { /* if we're close, realloc another INCR_DEGEN_LOG els */
                        *Degen_Log = (long *)realloc(*Degen_Log, sizeof(long)*((*Degen_Log)[0]+INCR_DEGEN_LOG));
                        if (*Degen_Log == NULL) BadMem("realloc Degen_Log", 1);
                        else (*Degen_Log)[0] += INCR_DEGEN_LOG;
                    }
                    (*Degen_Log)[DLi++] = seq_Cntr-1; /* set the start point */
                    (*Degen_Log)[1]++;  /* and keep track of where we are */
                    indegen = 1; /* don't forget to set this */
                }
            }
            if (F.Degen > 0) {
                if (F.Degen == 1) { /* 1 = default so maybe user forgot to set it */
                    F.Degen = 3; /* 3 = old -D flag; this should only visited once */
                    if (F.Degen == 1) fprintf(stderr, "\nIUPAC degeneracies in sequence - assuming '-D3'. \n"
                                                  "Type 'tacg -h' for brief description or 'man tacg' for more info.\n");
                }
                if (seq_Cntr++ >= begin) sequence [totSeqCnt++] = ct;   /* add it to the sequence */
                if (seq_Cntr > end) goto gotAllSeq; /* nasty but efficient way to break out of loops*/
            } else if (F.Verbose >= 1) { /* if verbose */
                if (f1++ == 0) {
                    fprintf(stderr, "\n!! NONdegenerate input mode '-D0', but the following valid IUPAC\n"
                            "characters are in the sequence (from headers, sequence, embedded comments, etc):\n");
                }
                if (ct != '\n' && ct != ' ') fprintf(stderr, "%c ", ct);
            }
            break;

        default:  /* bad char detect'n but if not acgt we don't count it and don't care - */
            if (F.Verbose >= 1) { /* unless we're being verbose, so spit bad chars to stderr */
                if (f1++ == 0) fprintf(stderr, "Bad characters: ");
                if (ct != '\n' && ct != ' ') fprintf(stderr, "%c ", ct);
            }
            break;
        }  /* end of switch/case statement */
    }  /* at this point 'totSeqCnt' should point just past the end of the sequence */

gotAllSeq: /* fast exit point */
    seq_len = totSeqCnt-BASE_OVERLAP-1; /* value to return to main() */
    /* the following commented stanza is useless until I rewrite the GetSequence routines to
       first grab ALL the sequence and THEN figure out where the -b and -e endpoints map in it.
       Since that isn't the way it's done right now and there are other more important things to do,
       I'm going to punt for now.

       if (F.SeqBeg != 1 && F.SeqBeg > seq_len-4){
          fprintf(stderr, "\nThe -b option specifies a start number either larger than the total sequence length or too close to the end.\n");
          exit(0);
       }
       if (F.SeqEnd != 0 && F.SeqEnd >= seq_len){
          fprintf(stderr, "\nThe -e option specifies an end number larger than or equal to the total sequence length.\n");
          exit(0);
       }
    */

    /* The following code performs the "pad the beginning sequence with the end and the end with the beginning"
       regardless of whether it's the whole sequence or a subsequence, using the vars 'begin'  and 'end' */

    /* and now pad out the beginning w/ the end and the end w/ the beginning of sequence[]  */
    /* I'm using 'BASE_OVERLAP' cuz it'll change if some RE is found with a longer recog seq offset*/
    /* l = index in  sequence' indicating absolute starting point for seq to be cut;
          starts at BASE_OVERLAP+1 (just past the buffer) if no subsequence chosen.  If a subseq was specified,
          it should be bumped up by the amount of 'begin', the user's idea of where the cut should start
          m = increment up to BASE_OVERLAP
          begin = beginning of the seq to be cut from the user's point of view.
          end = ending point for seq to be cut from the user's point of view.;
          was seq_len or totSeqCnt if no subsequence chosen */

    if (seq_len < BASE_OVERLAP)  { /* if it's a v short seq, pad it with junk */
        memcpy(&sequence[0], &blank, BASE_OVERLAP+1);   /* pad the beginning with a's */
        memcpy(&sequence[totSeqCnt], &blank, BASE_OVERLAP+1);  /* pad the back with a's */
    } else { /* else do the std tango of copy the end to the begin; begin to the end */
        for (l=BASE_OVERLAP+1,m=0; l <= BASE_OVERLAP*2+2 && m <= BASE_OVERLAP; l++,m++)  {
            sequence[totSeqCnt+m] = sequence[l]; /* copies seq past the end from the the beginning */
        }  /* note that in the above 'for', 'totSeqCnt' itself doesn't change */

        for (j=BASE_OVERLAP+1,m=0; j >= 0 && m <= BASE_OVERLAP; j--,l++,m++)  {
            sequence[m] = sequence[totSeqCnt-j]; /* copies seq to before the beginning from the end */
        }  /* note that in the above 'for', 'totSeqCnt' itself doesn't change */
    }

    if (F.SeqBeg != 1 || F.SeqEnd != 0) fprintf(stdout, "\nRaw Subsequence is %ld bases.\n\n", seq_len);

    totSeqCnt += BASE_OVERLAP;
    sequence[totSeqCnt+1] = '\0';  /* and term the sequence string, then decr it to point to end of the seq  */
    /* and free the mem not used by the sequence string */
    sequence = realloc(sequence, sizeof(char)*(totSeqCnt+1));
    /* sequence is padded with extra bases to allow circular cutting, but anything that depends on DNA length
       will reference the specific variable 'seq_len' */

    /* now we have the sequence in one long string, with 'BASE_OVERLAP' bp buffers at each end, so cut it as
       overlapping hexamers, starting at sequence[0], but only record cuts in the real seq */
    *tot_seq_Cnt = totSeqCnt; /* send it back to main() ... a greenie weenie way to do it */
    *seq_length = seq_len;       /*     ditto     */
    free(infile);
    return sequence;  /*  return the address of sequence ) */
}   /* End of GetSequence() */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   GetSequence() has been superceded by GetSequence2 which is does more and is better
   debugged.   To retrieve the code for GetSequence(), see above
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* *SEQ_len is the the exactly "Real Sequence Length" as determined from the SEQIO routines,
   no padding added or subtracted - it's incremented to totSeqCnt in here */

/* GetSequence2 reads and formats the sequence from stdin and returns the length of the
   sequence, as well as a pointer to the array that holds it (*sequence) along with the
   bracketing repeats - needs only a few extra variables to do it and it cleans up main
   considerably.  begin and end are indices to the real world coordinates of the sequence */

/* tot_seq_Cnt    length of sequence + bracketing overlaps
   seq_length     actual number of bases read
   Degen_Log		the array that logs the positions of the Degens
*/
char *GetSequence2(char *SEQseq, long *tot_seq_Cnt, long *SEQ_len, long *Degen_Log[])
{
    long m, j, l, seq_len, seq_Cntr=1, begin = 1, end = 1000000000, totSeqCnt;
    int  f1=0, iseq,
         DLi=2,        /* Degen_Log index starts @ 2 bc 1st 2 are mem counters */
         indegen=0;   /* 0 = not in degen, 1 IN a degen */
    char *sequence2, ct, blank[BASE_OVERLAP+1];
    memset(blank, 'a', BASE_OVERLAP+1);

    /* load the sequence2 array, filtering as we go */
    totSeqCnt = BASE_OVERLAP + 1;   /* to allow the wrapped seq to be inserted before the real seq
                                    (and to allow real seq to start at x1) */
    begin = F.SeqBeg;  /* set in SetFlags() */
    if (F.SeqEnd != 0) end = F.SeqEnd;    /* Ditto */

    /* know exactly how much mem to allocate for version 2 from SEQIO call */
    //if ((hjm = (char *)calloc((10000),sizeof(char))) == NULL) BadMem("kablooie!",1);
    if (F.Verbose > 0) {
        fprintf(stderr,"\nSeqFuncs.c:GetSequence2(): BEGIN sequence read..\n" );
    }
    sequence2 = 0X0;
    sequence2 = (char *)calloc((*SEQ_len + (2*BASE_OVERLAP) + 3), sizeof(char)); // this should allow a pad of 2 bytes at the end

    if (sequence2 == NULL) BadMem("Error alloc'ing *sequence2 in GetSequence2", 1);

    for  (iseq = 0; iseq < *SEQ_len; iseq++) { /*  while() replaced w/ for() cuz we know the limits */

        /* Now real seq will start at 'BASE_OVERLAP'; before that is wrapped seq from the
           end and at the end is wrapped seq from the beginning - required to handle circular seqs */
        ct = tolower(SEQseq[iseq]);
        if (ct == 'u') ct = 't'; /* converts u's to t's */
        /* this next construction assumes that all the sequence that's passed in from seqio is going to be
           GOOD sequence - probably a valid assumption, but should seed some bad seq into a file to
           check it thoroughly */
        if (seq_Cntr++ >= begin) { /* if we're not starting to slurp seq yet, don't do any of the rest */
            switch (ct) {
            case 'a': /* favored 4 always get added */
                F.NBases[0]++;
                goto NowProcess;
            case 'c':
                F.NBases[1]++;
                goto NowProcess;
            case 'g':
                F.NBases[2]++;
                goto NowProcess;
            case 't':
                F.NBases[3]++;

NowProcess:

                /* if its acgt && there was a degen beforehand, have to mark the end of that */
                /* this shouldn't need to be checked for mem as it just adds 1 value to the end */
                if (indegen == 1) {
                    /* Degen_Log[DLi++] = totSeqCnt-1;  */
                    (*Degen_Log)[DLi++] = seq_Cntr-2;
                    (*Degen_Log)[1]++;
                    indegen = 0;
                }

                sequence2[totSeqCnt++] = ct;   /* add it to the sequence */
                if (seq_Cntr > end) goto gotAllSeq; /* nasty but efficient way to break out of loops*/
                break;

            case 'y':
                F.NBases[4]++;
                goto NowProcessDegens; /* but if it's a valid IUPAC code */
            case 'r':
                F.NBases[5]++;
                goto NowProcessDegens;
            case 'w':
                F.NBases[6]++;
                goto NowProcessDegens;
            case 's':
                F.NBases[7]++;
                goto NowProcessDegens;
            case 'k':
                F.NBases[8]++;
                goto NowProcessDegens;
            case 'm':
                F.NBases[9]++;
                goto NowProcessDegens;
            case 'b':
                F.NBases[10]++;
                goto NowProcessDegens;
            case 'd':
                F.NBases[11]++;
                goto NowProcessDegens;
            case 'h':
                F.NBases[12]++;
                goto NowProcessDegens;
            case 'v':
                F.NBases[13]++;
                goto NowProcessDegens;
            case 'n':
                F.NBases[14]++;
                goto NowProcessDegens;

NowProcessDegens:

                /* regardless of F.Degen, may want to log degens, so .. */
                if (F.LogDegens) {
                    if (indegen == 0) {  /* if 0, 1st degen after a non-degen, so needs to be marked */
                        /* Degen_Log[even] = start, [odd] = end for all degens; doesn't distinguish among dif degens */
                        /* always check for end of alloc'ed mem */
                        if ((*Degen_Log)[0] - (*Degen_Log)[1] < 10) { /* if we're close, realloc another INCR_DEGEN_LOG els */
                            *Degen_Log = (long *)realloc(*Degen_Log, sizeof(long)*((*Degen_Log)[0]+INCR_DEGEN_LOG));
                            if (*Degen_Log == NULL) BadMem("realloc Degen_Log", 1);
                            else (*Degen_Log)[0] += INCR_DEGEN_LOG; /* keep track of the end of alloced mem */
                        }
                        (*Degen_Log)[DLi++] = seq_Cntr-1; /* set the start point */
                        (*Degen_Log)[1]++;  /* and keep track of where we are */
                        indegen = 1; /* don't forget to set this */
                    }
                }
                if (F.Degen > 0) {
                    if (F.Degen == 1) { /* 1 = default so maybe user forgot to set it */
                        F.Degen = 3; /* 3 = old -D flag; this should only visited once */
                        if (F.Verbose > 0) fprintf(stderr, "\nIUPAC degeneracies in sequence - assuming '-D3'. \n"
                                                       "Type 'tacg -h' for brief description or 'man tacg' for more info.\n");
                    }
                    /* should the degen logging go in here? */

                    sequence2[totSeqCnt++] = ct;   /* add it to the sequence */
                    if (seq_Cntr > end) goto gotAllSeq; /* nasty but efficient way to break out of loops*/
                } else if (F.Verbose > 0) { /* if verbose */
                    if (f1++ == 0) {
                        fprintf(stderr, "\n!! NONdegenerate input mode '-D0', but the following valid IUPAC\n"
                                "characters are in the sequence (from headers, sequence, embedded comments, etc):\n");
                    }
                    if (ct != '\n' && ct != ' ') fprintf(stderr, "%c ", ct);
                }
                break;

            default:  /* bad char detect'n but if not acgt (or a degen) we don't count it and don't care - */
                if (F.Verbose > 2) { /* unless we're being verbose, so spit bad chars to stderr */
                    if (f1++ == 0) fprintf(stderr, "Bad characters: ");
                    if (ct != '\n' && ct != ' ') fprintf(stderr, "%c ", ct);
                }
                break;
            }  /* end of switch/case statement */
        }
    }  /* end of the primary 'for' loop - at this point 'totSeqCnt' should point just past
         the end of the sequence */
    if (F.Verbose > 0) {
        fprintf(stderr,"\nSeqFuncs.c:GetSequence2(): END sequence read ..\nReal Sequence length = %ld\n\n", totSeqCnt );
    }

gotAllSeq: /* fast exit point */
    seq_len = totSeqCnt-BASE_OVERLAP-1; /* value to return to main() */

    /* The following code performs the "pad the beginning sequence with the end and the end with the beginning"
       regardless of whether it's the whole sequence or a subsequence, using the vars 'begin'  and 'end' */

    /* l = index in  sequence' indicating absolute starting point for seq to be cut;
          starts at BASE_OVERLAP+1 (just past the buffer) if no subsequence chosen.  If a subseq was specified,
          it should be bumped up by the amount of 'begin', the user's idea of where the cut should start
          m = increment up to BASE_OVERLAP
          begin = beginning of the seq to be cut from the user's point of view.
          end = ending point for seq to be cut from the user's point of view.;
          was seq_len or totSeqCnt if no subsequence chosen */

    if (seq_len < BASE_OVERLAP)  { /* if it's a v short seq, pad it with junk */
        memcpy(&sequence2[0], &blank, BASE_OVERLAP+1);   /* pad the beginning with a's */
        memcpy(&sequence2[totSeqCnt], &blank, BASE_OVERLAP+1);  /* pad the back with a's */
    } else { /* else do the std tango of copy the end to the begin; begin to the end */
        for (l=BASE_OVERLAP+1,m=0; l <= BASE_OVERLAP*2+2 && m <= BASE_OVERLAP; l++,m++)  {
            sequence2[totSeqCnt+m] = sequence2[l]; /* copies seq past the end from the the beginning */
        }  /* note that in the above 'for', 'totSeqCnt' itself doesn't change */

        for (j=BASE_OVERLAP+1,m=0; j >= 0 && m <= BASE_OVERLAP; j--,l++,m++)  {
            sequence2[m] = sequence2[totSeqCnt-j]; /* copies seq to before the beginning from the end */
        }  /* note that in the above 'for', 'totSeqCnt' itself doesn't change */
    }

    if (F.SeqBeg != 1 || F.SeqEnd != 0) fprintf(stdout, "\nSUBsequence is %ld bases.\n", seq_len);

    totSeqCnt += BASE_OVERLAP;
    /* and term the sequence string; not strictly nec as the mem was calloc'ed; therefore /everything/ is \0  */
    sequence2[totSeqCnt+1] = '\0'; //just to make sure....
    /* and free the mem not used by the sequence string - shouldn't be needed anymore because it was
       EXACTLY alloc'ed at the beginning, BUT, it SEQIO doesn't handle sequences the way I handle
       sequences, it could still be somewhat different - ie SEQIO could pass in spaces, numbers,
       newlines, etc and count them in the size of the array, so this should be left in until it can be
       absolutely verified that it's correct.  */

    /* sequence is padded with extra bases to allow circular cutting, but anything that depends on DNA length
       will reference the specific variable 'seq_len' */

    /* now we have the sequence in one long string, with 'BASE_OVERLAP' bp buffers at each end, so cut it as
       overlapping hexamers, starting at sequence[0], but only record cuts in the real seq */
    *tot_seq_Cnt = totSeqCnt; /* send padded seq cnt back to main() ... a greenie weenie way to do it */
    *SEQ_len = seq_len;       /* unpadded sequence count    ditto     */
    return sequence2;  /*  return the address of sequence ) */
}   /* End of GetSequence2() */


/* Degen_Calc ia a sleazy little fn() to calculate the degeneracy of the hexamer
   (possibly later, any string) without the added overhead of the other stuff in hash()
   that also does this but also makes uncomfortable calls to other more fragile functions  */

int Degen_Calc(char *s)
{
    int i, f1=0, TotDegen=1;
    for (i=0; i<6; i++) {
        switch (s[i]) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
            break; /* no increase in degeneracy */
        case 'y':
        case 'r':
        case 'm':
        case 'k':
        case 'w':
        case 's':
            TotDegen *= 2;
            break;
        case 'b':
        case 'd':
        case 'h':
        case 'v':
            TotDegen *= 3;
            break;
        case 'n':
            TotDegen *= 4;
            break;
        default:  /* bad char detect'n  */
            if (F.Verbose >= 1) { /* if we're being verbose, spit bad chars to stderr */
                if (f1++ == 0) fprintf(stderr, "Degen_Calc: Bad chars: ");
                fprintf(stderr, "%c ", s[i]);
            }
            break;
        }
    }
    return TotDegen;
}


/* SearchPaths() takes a filename and examines various environment variables to locate that
   filename. If it can find the file, it returns the full path name; if not, it returns NULL */

char  *SearchPaths (char *InputFileName, char *FileType4Err)
{

    struct stat *Fbuf;
    int i, File_Found=0, VERBOSE;
    char *FileToFind, *tacglib;
    static char *Env_Vars[] = {"PWD",     "$PWD (Current Directory)",
                               "HOME",    "$HOME (Home Directory)",
                               "TACGLIB", "$TACGLIB (TACG library)"
                              };
    VERBOSE = (int)F.Verbose; /* set VERBOSE */
    FileToFind = (char *) calloc (512, sizeof(char));
    if ((Fbuf = calloc(1, sizeof(*Fbuf))) == NULL) {
        BadMem("no mem for Fbuf",1);
    }
    /* specify an alternative input file - this is constructed so that if the name starts */
    /* with a '/', '~', or  '.' , assume it's a valid path name; otherwise 1st look in */
    /* $PWD, then in $HOME, then in $TACGLIB */

    /* 1st, see if it can be opened without any searching, regardless of what it starts with */
    strcpy(FileToFind, InputFileName);

    if ((stat(InputFileName, Fbuf)) == 0) { /* grab and check the next entry as the alt rebase file name */
        if (VERBOSE >= 1) fprintf(stderr,"Found %s file \"%s\"!!\n",FileType4Err, FileToFind);
        File_Found = 1;
    } else {  /* see if it's in either in the current directory $PWD (1st) $HOME (2nd) or $TACGLIB (3rd) */
        /* is it in $PWD ? */
        i=0;
        while (i < 5) {
            tacglib = getenv(Env_Vars[i++]);
            if (tacglib != NULL) {
                strcpy(FileToFind, tacglib); /* should overwrite FileToFind with tacglib */
                strcat(FileToFind, "/"); /* to allow strcpy'ing the path and the file together correctly */
                strcat(FileToFind, InputFileName); /* copy the file name in behind it */
                if ((stat(FileToFind, Fbuf)) == -1) { /*  if ! here, print err and check next env var */
                    if (VERBOSE >= 1) fprintf(stderr,"Hmmm.. %s file \'%s' not found in %s.\n", FileType4Err, FileToFind, Env_Vars[i]);
                    File_Found = 0;
                    i++;
                } else {
                    if (VERBOSE >= 1)  fprintf(stderr,"Found alternative %s file '%s' in %s.\n",FileType4Err, FileToFind, Env_Vars[i]);
                    i=100;
                    File_Found = 1;
                }
            } else if (VERBOSE >= 1)  fprintf(stderr, "Hmmm... %s isn't defined so I can't find it there!\n\n", Env_Vars[i++] );
        }
    }
    free(Fbuf);
    if (File_Found == 0) {
        free(FileToFind);
        FileToFind = NULL;
    }
    return FileToFind;
}   /* End of SearchPaths ()   */


/* compare is a dippy little function that qsort needs to perform its sordid little sort */
int compare(const void *n1, const void *n2 )
{
    return ( *((int *) n1) - *((int *) n2) );
}

/* abscompare is a dupe of compare above, except that it compares the absolute values
   of the variables instead of the given values */
int abscompare(const void *n1, const void *n2 )
{
    return ( labs(*((int *) n1)) - labs(*((int *) n2))  );
}



/********************************  Function Translate  *****************************
Translate() takes as input a pointer to a string of nondegenerate **OR degenerate**
DNA and returns the translated sequence as protein in either 1 or 3 letter
code in either 3 letter spacing or 1 letter spacing (for generating full length
protein strings.
For degenerate sequences, it's pretty dumbassed - if there's a degenerate base
in the triplet, it'll just ignore the sequence and move on.

Will translate according to a number of translation tables stored as
3D array of pointers in Codons [organism][codon][label].  The coordinate
points to a char string that holds the single letter xl at pos '0' and the
3 letter xl starting at pos 1.  Does not do 3 frame translation - have to
do 3 calls to it for that.  Also, does not do any output formatting - have
to handle that in the print section...  JUST translates.  However, will
truncate the translated string to the correct length if not a multiple of 3
**********************************************************************************/

void Translate (char *DNA_in, char *Prot_out, int len, int n_letters, char Codons[8][64][7], int organism)
{
    /* DNA_in   pointer to the DNA seq that needs to be translated, typically the entire DNA array,
                with the correct offset PREcalculated so that the pointer refers to the correct start
                site
       Prot_out pointer to the array that is filled by Translate() in either 1 or 3 letter code,
                typically passed to the print routines to be instantly printed out.
       len      the length of the DNA sequence to be translated
       n_letters   codes the number of letters of the aa label and the spacing:
                1 - A  L  H  V
                2 - ALHV
                3 - AlaLeuHisVal
       Codons   array that holds all the codon preference data for the different organisms or mitos
                loaded by subroutine Read_Codon_Prefs() in tacg.c from external file "codon.data"
                that has to be in the same dir as program currently.
       organism var that indicates which organism's codon prefs should be used for this translation
    */

    int sum[256], i, ii, match, j, degens, prot_len;
    char XXX[3] = {'X','X','X'};
    char  TestAA;
    if (len % 3 != 0) len = ((int) len /3)* 3;   /* round len to a multiple of 3 */
    prot_len = len;
    if (n_letters == 2) prot_len = len/3;  /* and div by 3 if we're creating the compact form of translation */
    memset (Prot_out,' ',prot_len); /* set Prot_out to blanks */
    for (i=0,j=0; i<len && j<len; i += 3,j++) {  /* for the length of the DNA sequence */
        if ((degens = hash(DNA_in+i, sum, 3)) == 1) {    /* hash the next triplet; if it's nondegenerate... */
            if (n_letters == 1)  { /* single letter, triple spacing for sub-DNA labeling */
                Prot_out[i] = Codons[organism][sum[0]][0]; /* look up the correct value and plunk it into the output */
            } else if (n_letters == 2)  {  /* single letters, with single spacing */
                Prot_out[j] = Codons[organism][sum[0]][0];
            } else { /* triple letters, triple spacing for sub DNA labeling */
                memcpy (&Prot_out[i], &Codons[organism][sum[0]][1],3); /* or memcpy the 3 letter code over */
            }
        } else { /* hash returns >1 */
            /* this needs to be modified to allow TRUE degenerate translation, like for validating '-silent'
               these should only  -> X || XXX if the hashes DO NOT point at the same AA.  ie, for reverse translation
               MOST of them WILL point to the same AA; only ARG, LEU, SER WILL NOT (using standard Codon Table)... */
            /* but only need to do extended checking for degenerates */
            match = 1;
            ii=0;
            TestAA = Codons[organism][sum[ii]][0];
            while (match == 1 && ii < degens ) {
                if (Codons[organism][sum[ii++]][0] != TestAA) {
                    match = 0;
                }
            }

            if (match == 0) {
                if (n_letters == 1)  { /* single letter, triple spacing for sub-DNA labeling */
                    Prot_out[i] = 'X'; /* look up the correct value and plunk it into the output */
                } else if (n_letters == 2)  {  /* single letters, with single spacing */
                    Prot_out[j] = 'X';
                } else { /* triple letters, triple spacing for sub DNA labeling */
                    memcpy (&Prot_out[i], &XXX,3); /* or memcpy the 3 letter code over */
                }
            } else {
                if (n_letters == 1)  { /* single letter, triple spacing for sub-DNA labeling */
                    Prot_out[i] = Codons[organism][sum[0]][0]; /* look up the correct value and plunk it into the output */
                } else if (n_letters == 2)  {  /* single letters, with single spacing */
                    Prot_out[j] = Codons[organism][sum[0]][0];
                } else { /* triple letters, triple spacing for sub DNA labeling */
                    memcpy (&Prot_out[i], &Codons[organism][sum[0]][1],3); /* or memcpy the 3 letter code over */
                }
            }
        }
    }
    /* and be sure to terminate it correctly */
    if (n_letters == 2) Prot_out[j] = '\0';
    else Prot_out[i] = '\0';
}  /* End of Translate()  */


/************************* NOT USED ANYMORE * Read_Codon_Prefs() *************************************************************************/
void Read_Codon_Prefs (char Codons[8][64][4], char Codon_Labels[8][20])
{
    /* Codons           array that holds all the codon preference data for the different organisms or mitos
                                    loaded from external file "codon.prefs"
       Codon_Labels   array that holds the organism labels (Universal, various Mitos, etc.
    */
    FILE *fpCodon;
    int e=0, i, j, k, sum[256];
    char triplet[3], labels[4], *ptmp, ctmp[512];
    for (i=0; i<3; i++)triplet[i] = 'a';

    /*   fprintf(stderr, "\n\n (BEFORE call to SearchPaths - 1st 200 bps of sequence:");
       for (i=0;i<200;i++) fprintf(stderr, "%c", sequence[i]);   */

    if ((ptmp = (SearchPaths("codon.prefs", "CODON PREFS"))) == NULL) {
        fprintf(stderr,"Can't open the CODON PREFS file!! - Check spelling, ownership, existance, etc - Bye!!\n");
        exit (1);
    } /* else  fprintf(stderr,"Going to try to read \"%s\"(in function Read_Codon_Prefs)!!\n", ctmp);     */

    /*   fprintf(stderr, "\n\n (AFTER call to SearchPaths - 1st 200 bps of sequence:");
       for (i=0;i<200;i++) fprintf(stderr, "%c", sequence[i]);   */


    /*  open the Codon Preference input file  */
    if ((fpCodon=fopen(ptmp,"r"))  ==  NULL)   {   /*  if  it's not there and readable */
        fprintf(stderr,"Cannot open the file \"codon.prefs\" for reading (in function Read_Codon_Prefs)!!\n");
        exit(1); /* print an error and die gracefully */
    }

    for (i=0; i<8; i++) { /* for each different Codon Pref stanza - also hard-coded to 8; will never change */
        e = fscanf(fpCodon,"%s", ctmp); /* stuff the stanza header string into the label array */
        if (e <= 0) {
            fprintf(stderr, "SeqFuncs.c:2962: Failure to read\n");
            exit(1);
        }

        strcpy (Codon_Labels[i],ctmp);
        for (j=0; j<64; j++) { /* for each doublet of values in the stanza (64 in all) */
            /* The next 2 lines ! work with gcc on Linux if the -O2 flag is used.  if no -O2 flag, it's OK ?!?!? */
            e = fscanf(fpCodon,"%3s %4s",ctmp, labels); /* grab the triplet and the corresp. label */
            if (e <= 0) {
                fprintf(stderr, "SeqFuncs.c:2968: Failure to read\n");
                exit(1);
            }
            /* fprintf(stderr,"fscanf returns : %d; The triplet is %s\n", l, ctmp);  */
            for (k=0; k<3; k++) triplet [k] = tolower(ctmp[k]); /* should change 'tolower' to a macro */
            hash(triplet,sum,3); /* hash the triplet to use as an index as to where to put the label */
            memcpy(&Codons[i][sum[0]][0],&labels,4); /* and copy the label to the right place in "Codons" */
        }
    }
    free(ptmp);
    fflush(stderr);
}


/*********************************  Function Anti_Par  ********************************
*   Anti_Par takes 2 char pointers, the 1st to the original seq, the 2nd to the       *
*   seq that it generates that is the reverse complement (formally, the Antiparallel  *
*   of the sequence), and the length of the seq to consider.   Compare with Rev_Compl *
*   below which would produce cctagtaaag from the sequence below                      *
*   if original = ggatcatttc, anti parallel = gaaatgatcc                              *
**************************************************************************************/
void Anti_Par(char *ori, char *anti, long len)
{
    /* ori   pointer to the beginning of the original seq
       anti  pointer to the beginning of the converted anti parallel sequence
       len   the length of the sequence, both original and converted
    */
    long chop = len-1, m;
    for (m=0; m<len; m++) {
        switch (ori[m]) {
        case 'a':
            anti[chop-m] = 't';
            break;
        case 'c':
            anti[chop-m] = 'g';
            break;
        case 'g':
            anti[chop-m] = 'c';
            break;
        case 't':
            anti[chop-m] = 'a';
            break;
        case 'r':
            anti[chop-m] = 'y';
            break;
        case 'y':
            anti[chop-m] = 'r';
            break;
        case 'w':
            anti[chop-m] = 'w';
            break;
        case 's':
            anti[chop-m] = 's';
            break;
        case 'm':
            anti[chop-m] = 'k';
            break;
        case 'k':
            anti[chop-m] = 'm';
            break;
        case 'b':
            anti[chop-m] = 'v';
            break;
        case 'd':
            anti[chop-m] = 'h';
            break;
        case 'h':
            anti[chop-m] = 'd';
            break;
        case 'v':
            anti[chop-m] = 'b';
            break;
        case 'n':
            anti[chop-m] = 'n';
            break;
        default:  /* bad character detection */
            fprintf(stderr,"Acck! In Anti_Par(), I don't like char# %ld= %c! \n",m, ori[m]);
            break;
        }  /* end of switch/case statement */
    }
    anti[len] = '\0';
}



/********************** Function Complement_In_Place  ********************************
*   Complement_In_Place takes a char pointer to the sequence to be converted and converts it
*   in-place, for the distance passed by 'len'.  Thus you can operate on sequences and
*   subsequences (by passing in the pointer to an offset and the length you want to
*   operate on) - if original = aggctgctggat,  Complement = tccgacgaccta
**************************************************************************************/
void Complement_In_Place(char *ori, long len)
{
    /* ori   pointer to the beginning of the original seq
       len   the length of the sequence to be operated on
    */
    long m;
    for (m=0; m<len; m++) {
        switch (ori[m]) {
        case 'a':
            ori[m] = 't';
            break;
        case 'c':
            ori[m] = 'g';
            break;
        case 'g':
            ori[m] = 'c';
            break;
        case 't':
            ori[m] = 'a';
            break;
        case 'r':
            ori[m] = 'y';
            break;
        case 'y':
            ori[m] = 'r';
            break;
        case 'w':
            ori[m] = 'w';
            break;
        case 's':
            ori[m] = 's';
            break;
        case 'm':
            ori[m] = 'k';
            break;
        case 'k':
            ori[m] = 'm';
            break;
        case 'b':
            ori[m] = 'v';
            break;
        case 'd':
            ori[m] = 'h';
            break;
        case 'h':
            ori[m] = 'd';
            break;
        case 'v':
            ori[m] = 'b';
            break;
        case 'n':
            ori[m] = 'n';
            break;
        default:  /* bad character detection */
            fprintf(stderr,"Acck! In Complement(), I don't like char# %ld= %c! \n",m, ori[m]);
            break;
        }  /* end of switch/case statement */
    }
    ori[len] = '\0';
}


/***********************  Function Rev_Compl  *****************************************
*   Rev_Compl takes 2 char pointers, the 1st to the original seq, the 2nd to the      *
*   seq that it generates that is the reverse complement, and the length of the seq   *
*   to consider.  Used to generate the matching bottom strand in the Linear Map       *
*   ie  ggatcatttc => cctagtaaag                                                      *
**************************************************************************************/
void Rev_Compl (char *orig, char *rev, long len)
{
    /* variables same as in Anti_Par, above */
    Anti_Par (orig,rev,len);
    Reverse (rev);
}


/********************  Function Reverse  ********************************
* straight from K+R (p62) - reverses a string s in place                *
************************************************************************/
void Reverse (char *s)
{
    /* s  pointer to the beginning of the array that holds the seq to be reversed */
    int i,j;
    char c;
    for (i=0,j=strlen(s)-1; i<j; i++,j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}


/* Function Triplet_Reverse reverses a string triplet by triplet ie:
   ArgTrpPheAsnCys ==> CysAsnPheTrpArg  so as to make the 6 frame translations
   readable in the oppo orientation */
void Triplet_Reverse (char *str)
{
    int i, j, mid, length, itmp;
    char tmp[3];
    length =  strlen(str);
    mid = (int)(length / 2);
    for (i=0; i<=mid; i+=3) {
        for (j=0; j<3; j++) tmp[j] = str[i+j];
        for (j=0; j<3; j++) {
            itmp = length - i - 3 + j;
            str[i+j] = str[itmp];
            str[itmp] = tmp[j];
        }
    }
}


/*************************  Function hash  *******************************************************
*   function hash takes a pointer to the n-mer string and generates the                          *
*   integer equivalent degeneracies that n-mer can expand to (max of 256 in                      *
*   this instance (with a hexamer, and returns a pointer to that array of numbers (sum, below),  *
*   along with the integer number of degeneracies generated                                      *
*************************************************************************************************/
int hash(char *nmer, int sum[256], int num)
{
    /* nmer  pointer to begin of array that holds the nmer to be 'hashed'
       sum   array (assigned in main() to be [256]) that holds the all the possible variants of the
             hashed sequence
       num   length of the nmer to be hashed; yes I know it's redundant...
    */
    int key[6] = { 1, 4, 16, 64, 256, 1024};
    int el, h, N_degen, t_degen, degen, i;
    degen = N_degen = 1;
    memset(sum,0,sizeof(int)*256);  /* set all of 'sum' to 0 - should work but any faster?*/

    /* Big 'for' loop that calculates the hexamer, gets executed at every RE site, */
    for (el = 0;  el < num;  el++) {       /* and at every overlapping n-mer */
        switch (nmer[el])  { /* pointer passed to function already offset to starting position */
            /* a=0, c=1, g=2, t=3, degenerates are handled below */
        case 'a':
            break;   /* not really needed - (a) = 0 so no change in sum */
        case 'c':
            for (i=0; i<degen; i++) sum[i] = sum[i] + key[el];
            break;
        case 'g':
            for (i=0; i<degen; i++) sum[i] = sum[i] + 2*key[el];
            break;
        case 't':
            for (i=0; i<degen; i++) sum[i] = sum[i] + 3*key[el];
            break;

            /*  !!!! Now the degeneracies !!!!  */
            /* Double degeneracies first */

        case 'y':   /* c or t  */
            N_degen = degen*2;
            t_degen = N_degen-1;
            fill_out_sum (degen, N_degen, sum);
            for (i=0; i<(degen); i++)  {
                sum[i] = sum[i] + key[el];   /* counting up sum, incr for 'c'*/
                sum[t_degen-i] = sum[t_degen-i] + 3*key[el]; /* counting down sum,*/
            }
            break;                                 /*    incr for 't' */
        case 'r':   /* g or a */
            N_degen = degen*2;
            fill_out_sum(degen, N_degen, sum);
            for (i=0; i<(degen); i++)  {
                /*  'a' will give 0 - no change */
                sum[i] = sum[i] + 2*key[el];    /* counting up sum, incr for 'g' */
            }
            break;
        case 'm':   /* a or c */
            N_degen = degen*2;
            fill_out_sum(degen, N_degen, sum);
            for (i=0; i<(degen); i++)  {
                /*  'a' will give 0 - no change */
                sum[i] = sum[i] + key[el];    /* counting up sum, incr for 'c' */
            }
            break;
        case 'k':   /* g or t */
            N_degen = degen*2;
            t_degen = N_degen-1;
            fill_out_sum(degen, N_degen, sum);
            for (i=0; i<(degen); i++)  {
                sum[i] = sum[i] + 2*key[el];                 /* counting up sum, incr for 'g' */
                sum[t_degen-i] = sum[t_degen-i] + 3*key[el];  /* counting down sum, incr for 't' */
            }
            break;
        case 's':   /* c or g */
            N_degen = degen*2;
            t_degen = N_degen-1;
            fill_out_sum(degen, N_degen, sum);
            for (i=0; i<(degen); i++)  {
                sum[i] = sum[i] + key[el];                /* counting up sum, incr for 'c' */
                sum[t_degen-i] = sum[t_degen-i] + 2*key[el];  /* counting down sum, incr for 'g' */
            }
            break;
        case 'w':   /* a or t */
            N_degen = degen*2;
            fill_out_sum(degen, N_degen, sum);
            for (i=0; i<(degen); i++)  {
                /*  'a' will give 0 - no change */
                sum[i] = sum[i] + 3*key[el];    /* counting up sum, increment for 'c' */
            }
            break;

            /* Triple degeneracies  */
        case 'b':   /* not a - so c, g, or t */
            h = 2*degen;
            N_degen = degen*3;
            fill_out_sum(degen, N_degen, sum);
            /* Increment all the array values in sum, based on the degeneracies */
            for (i=0; i<degen; i++)  {
                /*  'a' will give 0 - no change */
                sum[i] = sum[i] + key[el];    /* counting up sum, increment for 'c' */
                sum[i+degen] = sum[i+degen] + 2*key[el];    /* counting up sum, increm for 'g' */
                sum[i+h] = sum[i+h] + 3*key[el];    /* counting up sum, incr for 't' */
            }
            break;
        case 'd':   /* not c - so a, g, or t */
            h = 2*degen;
            N_degen = degen*3;
            fill_out_sum(degen, N_degen, sum);
            /* Increment all the array values in sum, based on the degeneracies */
            for (i=0; i<degen; i++)  {
                /*  'a' will give 0 - no change */
                sum[i+degen] = sum[i+degen] + 2*key[el];    /* counting up sum, incr for 'g' */
                sum[i+h] = sum[i+h] + 3*key[el];    /* counting up sum, incr for 't' */
            }
            break;
        case 'h':   /* not g - so a, c, or t */
            h = 2*degen;
            N_degen = degen*3;
            fill_out_sum(degen, N_degen, sum);
            /* Increment all the array values in sum, based on the degeneracies */
            for (i=0; i<degen; i++)  {
                /*  'a' will give 0 - no change */
                sum[i] = sum[i] + key[el];    /* counting up sum, increment for 'c' */
                sum[i+h] = sum[i+h] + 3*key[el];  /* counting up sum, incr for 't' */
            }
            break;
        case 'v':   /* not t - so a, c, or g */
            h = 2*degen;
            N_degen = degen*3;
            fill_out_sum(degen, N_degen, sum);
            /* Increment all the array values in sum, based on the degeneracies */
            for (i=0; i<degen; i++)  {
                /*  'a' will give 0 - no change */
                sum[i] = sum[i] + key[el];    /* counting up sum, increment for 'c' */
                sum[i+degen] = sum[i+degen] + 2*key[el]; /* counting up sum, incr for 'g' */
            }
            break;

            /* And the big old quadruple degeneracy  */
        case 'n':   /* a,c,g, or t */
            h = 2*degen;
            N_degen = degen*4;      /* t_degen = N_degen-1; */
            fill_out_sum(degen, N_degen, sum);
            /* Increment all the array values in sum, based on the degeneracies */
            for (i=0; i<degen; i++)  {
                /*  'a' will give 0 - no change */
                sum[i] = sum[i] + key[el];    /* counting up sum, increment for 'c' */
                sum[i+degen] = sum[i+degen] + 2*key[el];    /* counting up sum, increment for 'g' */
                sum[i+h] = sum[i+h] + 3*key[el];    /* counting up sum, increment for 't' */
            }
            break;
        default:  /* bad character detection */
//      fprintf (stderr, "Acck! in hash(), I don't like %c (at el = %d of %d)! \n", nmer[el], el, num);
            break;
        }  /* end of switch/case statement */
        degen=N_degen;
    }  /* end of big 'for' loop */
    return degen;
}


/***************************  Function Palindrome  ********************************
*   Function definition of palindrome - returns 1 if the sequence is a pal,       *
*   0 if it's not - fails on 1st nonpal character, so ought to fail quickly if    *
*   it is going to.                                                               *
**********************************************************************************/
int palindrome (char *site, int length)
{
    /* site     the site to be 'palindromed'; name should really be more generic
       length   the length of the sequence to be palindromed; yes, redundant, but convenient
    */
    int pal = 1, i = 0, halflength;
    if (length%2 == 0) halflength = length/2;
    else halflength=(length/2)+1;
    length--; /* to make the '*(site+(length-i)' expresssion work */
    while ((pal==1) && i<halflength) {
        switch (*(site+i)) {
        case 'a':
            if (*(site+(length-i)) != 't') pal=0;
            break;
        case 'c':
            if (*(site+(length-i)) != 'g') pal=0;
            break;
        case 'g':
            if (*(site+(length-i)) != 'c') pal=0;
            break;
        case 't':
            if (*(site+(length-i)) != 'a') pal=0;
            break;
        case 'y':
            if (*(site+(length-i)) != 'r') pal=0;
            break;
        case 'r':
            if (*(site+(length-i)) != 'y') pal=0;
            break;
        case 'm':
            if (*(site+(length-i)) != 'k') pal=0;
            break;
        case 'k':
            if (*(site+(length-i)) != 'm') pal=0;
            break;
        case 'w':
            if (*(site+(length-i)) != 'w') pal=0;
            break;
        case 's':
            if (*(site+(length-i)) != 's') pal=0;
            break;
        case 'b':
            if (*(site+(length-i)) != 'v') pal=0;
            break;
        case 'd':
            if (*(site+(length-i)) != 'h') pal=0;
            break;
        case 'h':
            if (*(site+(length-i)) != 'd') pal=0;
            break;
        case 'v':
            if (*(site+(length-i)) != 'b') pal=0;
            break;
        case 'n':
            if (*(site+(length-i)) != 'n') pal=0;
            break;
        default:
            fprintf(stderr,"palindrome doesn't like %c! \n", *(site+i));
            break;
        }  /* end of switch/case statement */
        i++;
    }
    return pal;
}


/***************************  Function fill_out_sum  *************************************
*   Function definition for fill_out_sum - duplicates the degeneracy the correct # of    *
*   times in the array 'sum[]'.                                                          *
*****************************************************************************************/
void fill_out_sum (int O_dgen, int N_dgen, int s[256])
{
    /* O_dgen   old degeneracy that needs to be updated to the new degeneracy defined in...
       N_dgen   new degeneracy that's calculated here
       s        local version of sum that holds all the degeneracy values calculated here  */
    int i = 0, j = 0;
    for (i=O_dgen; i<N_dgen; i=i+O_dgen) {
        do  {
            s[i+j] = s[j];
            j++;
        } while (j<O_dgen);
        j = 0;
    }
}



/* DownCase() downcases all the letters in the submitted string and returns the pointer
	to the same (now downcased) string -it CHANGES the original string, so be careful.*/
char *DownCase(char *string)
{
    int i=0, len=0;
    /* check that there's something to downcase (thanks to Nicolas Joly for this
       bit) */
    if (string == NULL) {
        fprintf(stderr, "Boom!! in DownCase(), NULL string detected!\n");
        exit(1);
    }
    len = strlen(string);
    for (i=0; i<len; i++)  string[i] = tolower(string[i]);
    string[i] = '\0'; /* to terminate it properly */
    return string;
}


/* realhash() is a real hashing function, straight from Sedgewick (p233) that returns
        a good hash value within the scope of the allocated space, because of the mod
        fn. tablesize MUST be prime for reliable results */
unsigned realhash(char *string2hash, int tablesize)
{
    int hash;
    for (hash=0; *string2hash != '\0'; string2hash++) {
        hash = (64 * hash + *string2hash) % tablesize;
    }
    return hash;
}


/* Usage() just spits out a few lines about how to use the program - shouldn't be more than
        the smallest terminal screen (at a time) that program is expected to see - say,
        80x25. Well, in the old days it was..     */
void Usage(void)
{
    // if move this chunk out, needs the #include <termios.h> at top.
    // all this stuff down to the 'int i=0;' sets up the termio to grab single chars at a time
    //without waiting for a <ENTER> to terminate it.
    struct termios orig,now;
    setvbuf(stdout, NULL, _IONBF ,0);
    tcgetattr(0, &orig);
    now=orig;
    now.c_lflag &= ~(ISIG|ICANON|ECHO);
    now.c_cc[VMIN]=1;
    now.c_cc[VTIME]=2;
    tcsetattr(0, TCSANOW, &now);

    int i=0;
    IndexScreen();
    while ((i = getchar())) {
        switch (i) {
        default:
            fprintf (stdout, "Not a choice; choose 1-7, '(A)ll', '(i)ndex', or '(q)uit'\n");
            break;
        case 'i':
        case 'I':
            IndexScreen();
            break;
        case '\n':
            break;
        case '1':
            HelpScreen1();
            break;
        case '2':
            HelpScreen2();
            break;
        case '3':
            HelpScreen3();
            break;
        case '4':
            HelpScreen4();
            break;
        case '5':
            HelpScreen5();
            break;
        case '6':
            HelpScreen6();
            break;
        case '7':
            HelpScreen7();
            break;
        case 'q':
        case 'Q':
            if(tcsetattr(0, TCSAFLUSH, &orig) < 0) { // reset the terminal to something sane
                fprintf (stderr, "Can't reset terminal?!; This would be a bug.  Please report it.\n");
                exit(-1);
            }
            exit(0);
            break;
        case 'A':
        case 'a':
            HelpScreen1();
            HelpScreen2();
            HelpScreen3();
            HelpScreen4();
            HelpScreen5();
            HelpScreen6();
            HelpScreen7();
            break;
        }
    }
}

void IndexScreen(void)
{
    fprintf (stderr, "Usage:   tacg -flag option -flag option ... <infile >outfile\n"
             "tacg uses stdin/stdout/stderr; uses redirection or pipes for input and output;\n"
             "needs input specifier (| or <); output to screen (default), >file, | nextcmd\n"
             "Uses Knight's SEQIO for auto reformat on input; most ASCII formats accepted.\n"
             "1 or more of: -F -g -G -l -L -O -p -P --rules --rulefile -s -S -X flags must be\n"
             "specified for output.\n"
             "+-------+---------+-----------+--------+----------+-----------+----------+-------+\n"
             "| 1     |    2    |    3      | 4      |    5     |     6     |      7   |   A   |\n"
             "+-------+---------+-----------+--------+----------+-----------+----------+-------+\n"
             "| b     | dam     | i         | m/M    | rule     | silent    | X        |  ALL  |\n"
             "| e     | dcm     | infile    | n      | rulefile | tmppath   | #        | Pages |\n"
             "| c     | example | logdegens | o      | r regex  | T         | examples |       |\n"
             "| C     | f       | l         | O      | R        | v         +----------+------ +\n"
             "| clone | F       | L         | orfmap | raw      | V         | rev      |   i   |\n"
             "| cost  | g       | strands   | p      | s        | w         | comp     +-------+\n"
             "| D     | G       | notics    | P      | S        | W         | revcomp  | index |\n"
             "|       | h help  | numstart  | ps     |          | x         |          |       |\n"
             "|       | H HTML  |           | pdf    |          |           |          |       |\n"
             "+-------+---------+-----------+--------+----------+-----------+----------+-------+\n"
             " (leading dashes for flags above have been removed for brevity)\n"
             "Type one of the numbers indicated to show the help page for that flag.\n"
             "1-7 = page for those flags, 'A' = all pages, 'i' = index page, 'q' = quit\n"
            );

}


void HelpScreen1(void)
{
    ScreenSep();
    fprintf (stdout,
             "-b   {#}       beginning of DNA subsequence; 1* for 1st base of sequence.\n"
             "-e   {#}       end of DNA subsequence; 0* for last base of sequence.\n"
             "-c             order output by # of hits per Enz, else by order in REBASE file.\n"
             "-C   {0*-16}   Codon Usage table to use for translation:\n"
             "               0 Standard*      6 Echino_Mito      12 Blepharisma_Nucl\n"
             "               1 Vert_Mito      7 Euplotid_Nucl    13 Chlorophycean_Mito\n"
             "               2 Yeast_Mito     8 Bacterial        14 Trematode_Mito\n"
             "               3 Mold_Mito      9 Alt_Yeast_Nucl   15 Scenedes_Mito\n"
             "               4 Invert_Mito   10 Ascidian_Mito    16 Thrausto_Mito\n"
             "               5 Ciliate_Nucl  11 Alt_Flatworm_mito\n"
             "               (read file 'codon.data' for complete information)\n"
             "--clone {#_#,#x#..} find REs that don't cut in the range #_#, but do cut in\n"
             "             the range #x#. Returns RE names & sites matching each criteria.\n"
             "--cost {#}   use only REs that cost >= # units/$. Larger #s are cheaper.\n"
             "             Requires optional data be added to rebase file; REs lacking this\n"
             "             info are excluded.  #s>100 are cheap; #s<10 are v. expensive.\n"
             "-D {0|1*-4}    controls input and analysis of degenerate sequence where:\n"
             "          0  FORCES excl'n IUPAC degen's in sequence; only 'acgtu' accepted.\n"
             "          1* cut as NONdegenerate unless IUPAC char's found; then cut as '-D3'.\n"
             "          2  allow IUPAC char's; ignore in KEY hex, but match outside of KEY.\n"
             "          3  allow IUPAC char's; find only EXACT matches.\n"
             "          4  allow IUPAC char's; find ALL POSSIBLE matches.\n"
             "1-7 = page for those flags, 'A' = all pages, 'i' = index page, 'q' = quit\n"
             "\n"
            );
}

void HelpScreen2(void)
{
    ScreenSep();
    fprintf (stdout,
             "--dam          simulate cutting in the presence of Dam methylase (GmATC)\n"
             "--dcm          simulate cutting in the presence of Dcm methylase (CmCWGG)\n"
             "               for both above, REs not cutting due to methylation are listed as\n"
             "               not cutting at all in summary.\n"
             "--example {1-10} example code to show how to add your own flags and functions.\n"
             "               Search for 'EXAMPLE' in 'SetFlags.c' and 'tacg.c' for the code.\n"
             "-f   {0|1*}    form (or topology) of DNA - 0 (zero) for circular; 1 for linear.\n"
             "-F   {0*-3}    print/sort Fragments; 0*-omit; 1-unsorted; 2-sorted; 3-both.\n"
             "-g   {Lo(,Hi)} prints a gel map w/ lo cutoff of Lo, (opt) hi cutoff of Hi bp.\n"
             "               Lo&Hi can be any value to slice out the range of interest.\n"
             "   NB: if Hi is omitted, it will be set to the the seq length, but it can be \n"
             "       set to ANY length to provide identical output conditions for comparisons\n"
             "-G {#,X|Y|L}   streams numeric data to stdout for external analysis/plotting.\n"
             "   # = bases/bin (the hits for this many bases should be pooled).\n"
             "   X = bins on X axis; Y = bins on Y axis; L = Long output as 'bins(X) data(Y).\n"
             "-h (--help)    asks for (this) brief help page.\n"
             "-H (--HTML) {0*|1}  complete (0) or partial(1) HTML tags generated on the fly \n"
             "               for WWW output. See man page for appro usage. \n"
             "               0 = makes standalone HTML page, with Table of Contents.\n"
             "               1 = no page headers, only TOC, to embed in other HTML pages.\n"
             "\n"
            );
}

void HelpScreen3(void)
{
    ScreenSep();
    fprintf (stdout,
             "-i (--idonly) {0|1*|2} controls output for seqs that have no hits.\n"
             "   0 - ID line and normal output printed regardless of hits.\n"
             "   1 (default) ID line and normal output are printed ONLY IF there are hits.\n"
             "   2 - ONLY ID line is printed if there are hits.\n"
             "--infile {input filename}  OR tacg < input_filename OR stdout | tacg \n"
             "               the above are alt methods to specify an input to tacg.\n"
             "-l             prints a GCG-style ladder map.\n"
             "-L             print a Linear map - produces LOTS of output (~10x input).\n"
             "--logdegens    all degens logged for graphic output (mem intensive).\n"
             "               can be made more space efficient by using next 2 options.\n"
             "--strands {1|2*} in Linear map, print 1 or 2 strands per line of DNA.\n"
             "--notics       in Linear map, omit the tics that indicate 5, 10 bases.\n"
             "--numstart {#} in Linear map, start numbering at this number (+ or -).\n"
             "\n"
            );
}

void HelpScreen4(void)
{
    ScreenSep();
    fprintf (stdout,
             "-o {0|1*|3|5,#} overhang - 5=5', 3=3', 0 for blunt, 1(d) for all. If you append\n"
             "                ',#' (1<#<7), filters on the overhang length as well: '-o5,4' \n"
             "-O {###(x),min}  ORF table for selected frames; 'x' -> extra info re AA compos'n\n"
             "               ie: -O 135x,25 = fr 1,3,5 w/ xtra info; min ORF len = 25aas.\n"
             "               'x' -> 3 extra lines x 128 chars wide, mungs formal FASTA format.\n"
             "--orfmap       prints a pseudo graphic map showing the locations of all ORFs \n"
             "               requested with the -O flag (see above) & a MET / STOP map\n"
             "               that shows only MET and STOP codons.  Expand scale with '-w'\n"
             "-p {Label,pattern[,Err]} cmd line entry of (degenerate) patterns to search for\n"
             "               if Err is missing, it is set = 0, also sets -S for output.\n"
             "  eg:  -pFindMe,gyrttnnnnnnngct,1  looks for indicated pattern with 1 error.\n"
             "-P {Lab1,[+-lg]DistLo[-DistHi],Lab2}  Proximity match'g for 2 named patterns.\n"
             "               Lab1/2 patterns must be in a REBASE-format file in form:\n"
             " 'Lab1 1 IUPAC_pattern 0 Err !Comments '  where Err = max # errors allowed eg:\n"
             " 'FindMe 1 gyrttnnnnnnngct 0 1 !The pattern that I'm trying to find \n"
             "               can repeat to specify up to 10 relationships at once\n"
             "               + (-) Label1 is downstream (upstream) of Label2; default either\n"
             "               l (g) Label1 is < (>) or = to 'DistLo' from Label2\n"
             " 'DistLo-DistHi'  indicates an explicit distance range (obviating l,g)\n"
             "--ps           writes Postscript plasmid map (tacg_Map.ps), \n"
             "               forces circular DNA, notes degens around rim, can be \n"
             "               combined with -O (above) to plot ORFs in any frames. Will do \n"
             "               multiple pages with a multi-sequence file.\n"
             "--pdf          converts above PS plasmid map to PDF (tacg_Map.pdf) via exec\n"
             "               of 'ghostscript', if it exists. If it doesn't tacg will exit. \n"
             "-m/M {#}       minimum (-m) and/or Maximum (-M) # cuts/RE; 0* for all.\n"
             "-n {3*-10}     magnitude of recognition site; 3 = all, 5 = 5,6,7....\n"
             "\n"
            );
}

void HelpScreen5(void)
{
    ScreenSep();
    fprintf (stdout,
             "--rule 'RuleName,((LabA:m:M&LabB:m:M)|(LabC:m:M&LabD:m:M))|(LabE:m:M),window'\n"
             "               explicitly selects named patterns (<16) from a REBASE file with\n"
             "               per pattern min/Max limits. 'RuleName' = name for the pattern\n"
             "               len<=11 (unless #define changed), 'window' = sliding window \n"
             "               within which the rule must be true.\n"
             "               Parens () enforce logic; otherwise expressions are eval'd\n"
             "               L->R. LabX = Pat name; 'm' = min; M = Max; '&' = logical AND;\n"
             "               '|' = OR; '^' = XOR.  Enclosing single 'quotes' are REQUIRED.\n"
             "               Valid patterns are logged to file 'tacg.patterns' in current dir\n"
             "               for re-use. See manpage for info, examples\n"
             "--rulefile '/path/to/rulefile'\n"
             "              loads a series of complex rules like those described in --rule\n"
             "              above.  Format is as in the single quotes above : \n"
             "                 RuleName,(rulestring, as above),window \n"
             "-r (--regex) {'Label:RegexPat'} search for RegexPat; use 'Label' for naming;\n"
             "    translates IUPAC characters into std regex notation: \n"
             "                gy(tt|gc)nc{2,3}m -> g[ct](tt|gc).c{2,3}[ca] \n"
             "tacg uses Phil Hazel's pcre library and therefore perl syntax\n"
             "NB: regex-matching is incompatible with regular IUPAC and matrix-matching.\n"
             "-r (--regex) {'FILE:FileOfRegexPats'} open the FILE 'FileOfRegexPats' and search\n"
             "    for all Regex pat's in it; 'FILE' must be in CAPS to trigger this behavior.\n"
             "NB: -r REQUIRES the ':' separator and single quotes ' to enclose options.\n"
             "-R {alt pattern file} specifies alternative REBASE file in GCG format\n"
             "               OR use it to specify the MATRIX data file in TRANSFAC format.\n"
             "--raw          ANYTHING on STDIN is raw sequence; same behavior as pre-SEQIO.\n"
             "-s             summary - print Table of Zero Cutters, # hits of each Enzyme.\n"
             "-S {1*|2}      Sites - prints the the actual cut Sites in tabular form.\n"
             "               1* = sites noted as + offsets; fine for restriction mapping.\n"
             "               2  = note nonpals on bottom strand with '-' offsets.\n"
             "\n"
            );
}
void HelpScreen6(void)
{
    ScreenSep();
    fprintf (stdout,
             "--silent       searches for possible SILENT RE sites (those that won't cause\n"
             "           Xl'n to change.  use -LT1,1 to see rev-trans sequence re-trans'd\n"
             "           to verify that the seq is OK.  Arg, Leu, Ser codons will not RT/FT.\n"
             "           NB: this rev Xlates sequence affecting ALL results using this flag!\n"
             "--tmppath  passes a temporary path for cooperation with CGIs that need output\n"
             "           in a particular place. Currently used in plasmid map generation.\n"
             "-T {[0*|1|3|6],[1|3]} Translates frames 1, 1-3, 1-6 w/ Linear Map using 1 or 3\n"
             "               letter labels. -T3,3 xlates frames 1,2,3 with 3 letter labels.\n"
             "-v             prints version of the program, then dies.\n"
             "-V   {1-3}     Verbose/debug mode; spews diags to stderr (1=lots, 3=tons).\n"
             "-w   {1|#}     output width in bp's (60 < # < 210), truncated to a # mod 15\n"
             "               '-w 1' = 1 line output for easier parsing by external apps.\n"
             "-W (--slidwin) {#} Defines the sliding window for searching for pattern groups\n"
             "        use w/ '-x' as a looser alternative to -P (Proximity matching in pairs)\n"
             "-x {Label(,=),Label..(,C)} selects SPECIFIC REs (<=15) from the REBASE file;\n"
             "               If ',=' is appended to 1 RE Label, it will tag that RE for the\n"
             "               AFLP analysis.  See man page for details.  If ',C' is\n"
             "               appended to a list of >1 Labels, requests a multiple digest.\n"
             "\n"
            );
}

void HelpScreen7(void)
{
    ScreenSep();
    fprintf (stdout,
             "-X (--extract) {b,e,[0|1]} eXtracts the sequence around the pattern matched, \n"
             "    from b bases preceding, to e bases following the MIDDLE of pattern (if a \n"
             "    normal pattern, the START of the pattern if a regular expression. If the\n"
             "    pattern is found in the bottom strand AND the last field = 1, sequence is \n"
             "    anti-par'ed before it's extracted so all patterns are in same orientation;\n"
             "    if last field = 0, it is NOT reverse compl'ed.\n"
             "-#   {#}       Matrix matching, using TRANSFAC-format matrix input. Required #\n"
             "               is the matrix cutoff as a %%.  Uses '-R' to specify alternative\n"
             "               matrix file, if not 'matrix.data'. Also works with '-x'. \n"
             "NB: matrix-matching is incompatible with regular pattern and regex matching.\n"
             "--rev     reverses sequence(s) before analysis:  tacg -> gcat \n"
             "--comp    complements             \"             tacg -> atgc \n"
             "NB: above 2 options are useful in resolving sequencing or orientation errors\n"
             "--revcomp reverse complements     \"             tacg -> cgta \n"
             /*       "--revnum reverses numbering for some outputs (-L) \n"  */
             "ex:       tacg -f0 -n6 -T123,3 -sl -F3  <degen.input.file >output.file (to file)\n"
             "          tacg -f 0 -l  -n 5  -F 2 <input.seq.file  (to stdout/screen)\n"
             "          tacg -m 3 -T1,1 -s  |grep HindIII < Ecoli_genome.genbank >out\n"
             "          tacg --regex 'Rxname1:g(cc|tac)nr{2,4}yt <somefile  >outfile\n"
             "          tacg -x HindIII,EcoRV,bamhi,C -O 134,30 -w 90 -LlS <input.seq.file\n"
             "for more help, try 'man tacg'. Type 'Ctrl+C' if the program seems locked.\n"
             "Latest docs & code at: http://tacg.sf.net. Contact author:  hjm@tacgi.com\n"
             "\n"
            );
}

void ScreenSep(void)
{
    fprintf (stdout,
             "flag | opts  |  explanation (* = default; # = an integer)\n"
             "-----+-------+--------------------------------------------------------------\n");
}

