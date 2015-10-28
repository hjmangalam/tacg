/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: ReadEnzFile.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

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

/* ReadREs does just that - given a filename, it tries to open it and if successful, reads
   all the uncommented REs into the global RE[] struct, and then returns the # of entries in RE
   (NumREs).  Matching against the SELECTED REs specified on the command line takes place after
   this functions returns.  */

/* In Version 2, it also does a whole lot more, mostly related to creating the RE entries
   nec for doing the ERROR matching.  This involves reading a file, either a named REBASE
   formatted file or a temp file created in SetFlags to carry the patterns entered via the -p
   flag.  These entries are parsed as usual, but if there's an error term at the end, all the
   possible error terms are generated on the fly and entered into separate RE entries.  Each of
   these is compared with all previous RE entries back to the prototype to see if it's
   identical (fast and sloppy vs slow and careful, but it will match the majority without
   allowing too many duplicates).  When it verifies that a generated pattern is really novel,
   it appends it.  For each error allowed, it generates an 'n' scan, much like a linker scan
   across the sequence, replacing each base with an n.  Obviously for Err=2, it will do this 2x
   etc.  This is where a lot of the duplication occurs.  It also catches palindromes and only
   checks them once.  It also trims leading and lagging n's to decrease the size of the
   sequence and further reduce the size of succeeding scans.  As with the regular handling, it
   also finds the most significant hexamer profile in the sequence and uses the ~midpoint of
   that hex as the indicator site.  This reduces the amount of offset calc'ns and also
   indicates the most significant match more easily.  There is an inefficiency in the
   backchecking as it checks ALL the previously stored entries instead of just the ones which
   share the same hexamer hash.  This inefficiency does not show up until you ask for lots of
   degeneracy in large sequences, so I'm not going to make it more efficient just yet.  */

/* RebaseFile can be either the default or an optional filename - decided in/immed after SetFlags().
   RE[] is now global, so no need to declare it.  Fn() returns a pointer to hashtable[4096] which
   is composed in this fn() via a series of strange manipulations.  NumREs is now 'returned'
   via a pointer reference in the calling arguments */

/* this can now read either the old style GCG files, or can now read a tacg-extended GCG format that
   includes Errors Allowed in the pattern, dam-sensitivity, dcm-sensitivity,  units, and $ for that
   amount of enyme (file supplied with the high unitage amount from NEB, but you can supply any
   values you want).  The last 2 values are used to calculate the 'cost' factor and this can be used
   as a filter for selecting REs, along with overall min/Max, per RE min/Max, overhang, magnitude of
   recognition site, etc. */

int ReadEnz(char *EnzFileName, FILE *tmpfp, int *NumREs, struct SE_struct *SelEnz,
            int *hashtable[4096], long SeqLen)
{

    int e, i, j, k, l, m, n, p, z, lnct, cost,/* some of these counters can prob be re-used */
        r[5], /* 0 = err; 1 = dam; 2 = dcm; 3 = cost; 4 = #units */
        dpl=0, ori=1,  /* these are the indices for the doppel and the original sequences -
                        the 0 and 1 are reversed from what you might ordinarily expect
                        (o for ori, 1 for doppel) so that a 'for' loop will count up
                        correctly to the right # - there will always be a 0 (ori),
                        sometimes a 1 (doppel) */
                   EORE,          /* End of RE - marks and tracks the mem used by RE */
                   GotH = 0,      /* Indicator to say that I got the hookey RE */
                   Nprotos=0,     /* # of protos to return at end */
                   len, Degen, OLapLen=0, rep=0, rsl, thrashtable[4096],
                               itemp1=0, itemp2=0, dgn_sits[256], BOStanza, EOStanza,
                               SelREs, NumSelREs=0, RENameHash[RE_NAME_TAB_SIZ],
                                       SEi, Err=0, proto, here, E, same, PPi;

    /* thrashtable is a temp holder for hashtable data to get the numbers; afterwards all the
       info will be used to transferred to hashtable[], the pointer to which will be returned to
       main() and used thereafter.  thrashtable[] storage will die with this function, so its mem
       will be returned to the pool.  */

    float fRE_mag=0;
    char *ctemp1, *ctemp2, *RE_rawsite, *rstr, *tmp, ct, *tok, junk[10];

    FILE *fpinRE;
    memset(RENameHash,0,sizeof(int)*RE_NAME_TAB_SIZ);
    memset(dgn_sits,0,sizeof(int)*256);
    memset(thrashtable,0,sizeof(int)*4096);

    /* mem for rstr, rstrC */
    if ((rstr  = (char*)calloc(256, sizeof(char))) == NULL) BadMem("rstr", 1);
    tmp = rstr; /* IMPORTANT - this saves the beginning of the alloc'ed space to allow rstr to be reset to it
                  after it's chewed to pieces by strsep */

    if ((ctemp1 = calloc(256,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'ctemp1'", 1);
    }
    if ((ctemp2 = calloc(256,sizeof(char))) == NULL) {
        BadMem("Can't init mem for 'ctemp2'", 1);    // ctemp2 is the disposable version of ctemp1
    }


    if (F.Xplicit == -1 || F.Rules == 1) {
        SelREs = 1;
    }
    SelREs = (int)abs(F.Xplicit); /* this way because both '-r' and '-P' (puts a -1 in F.Xplicit) use SelEnz[] */
    /*  *NumREs=1;   this has to start at 1 rather than 0, because negation is used later to mark the doppels */
    EORE = INIT_RES;

    /* If the -x, -r, --rules, and -P flag has been set, hash the names into the name space */
    if (SelREs == 1) { /* -x, -r, --rules, and -P use SelEnz[] - should act as the correct input filter for both */
        while ((SelEnz[NumSelREs].PName[0] != '\0') && (NumSelREs < MAX_SEL_ENZ)) { /*   */
            RENameHash[realhash(SelEnz[NumSelREs].PName, RE_NAME_TAB_SIZ)]++; /* incr the index */
            NumSelREs++;
        } /* and reset the '-n' and '-o' flags to default, regardless of what they were set to in SetFlags() */
        F.Mag = 3;
        F.Overhang = 1;
    }

    /* open the REbase input file, specified either by flag or by setting above */
    if (F.Pat == 1) { /* either the temp file opened and fed in SetFlags() */
        fpinRE = tmpfp; /* point one to the other */
        rewind(fpinRE); /* and rewind it */
    } else if ((fpinRE=fopen(EnzFileName,"r")) == NULL) {  /* or the standard/optional one from SetFlags() */
        if (EnzFileName == NULL) {
            EnzFileName = "rebase.data";
        }
        fprintf(stderr,"\nCannot find or open the REbase file \"%s\" for reading!!\n"
                "(There's no readable-by-you 'rebase.data' file in this directory,\n"
                "your home directory, or in TACGLIB.\n"
                "You DID set the TACGLIB environment variable, didn't you?\n\n", EnzFileName); /* DON'T HIDE behind -V */
        exit(1); /* print an error and die gracefully */
    }

    /*  Scan thru comments in rebase file to separator ("..")   */
    junk[0] = 'k';   /* make sure the strcmp below fails the first time */
    while ((feof(fpinRE) == 0) && (strncmp(junk, "..",2) != 0))  {
        e = fscanf (fpinRE, "%2s", junk);
        if (e <= 0) {
            fprintf(stderr, "ReadEnzFile.c:138 Failure to read\n");
            exit(1);
        }
        ct = 'm'; /* and read to end of line - this needs to be able to deal */
        while (ct != '\n' && ct != '\r' && (feof(fpinRE) == 0)) {  /* with end of line conditions for Mac, Win, and unix - what have I missed? */
            ct = fgetc(fpinRE);
        }
    }

    /* we're at "..", so go into file and start slurping directly data into struct
       except for RE_rawsite which has to be filtered a couple ways before being acceptable */

    ct = 'm';
    /* using fgets / strsep to read in variables, instead of the cruddy old way w/ fscanf */
    while (fgets(rstr, 256, fpinRE) && (*NumREs < MAX_NUM_RES)) {
        /* now the whole line's in rstr; use strsep to parse it, make decisions re: parsing, skipping */
        if (rstr[0] != ';' && rstr[0] != '\n') {  /* if the enz hasn't been commented out, get all the bits */
            /* REMEMBER!  strsep chews the string to pieces and changes the start of the string */
            ctemp1 = strsep(&rstr, " \t"); /* using ctemp1 (=Name) for congruence with previous code  */
            if (strlen(ctemp1) > (MAX_PAT_NAME_LEN)) {
                fprintf(stderr, "Name %s is longer than allowed (%d) - check it!\n", ctemp1, MAX_PAT_NAME_LEN);
                exit(1);
            }
            tok = strsep(&rstr, " \t");
            while (tok[0] == '\0') tok = strsep(&rstr, " \t");
//			fprintf(stderr, "F.Pat=%d\n", F.Pat);
            if (F.Pat == 1) {
                itemp1 = 0;
            } else {
                itemp1 = atoi(tok);
            }

//         itemp1 = atoi(tok); /* next one - top offset, prev -> itemp1 */
            tok = strsep(&rstr, " \t");
            while (tok[0] == '\0') tok = strsep(&rstr, " \t");
            if (strlen(tok) > BASE_OVERLAP+2) {
                fprintf(stderr, "Recognition pattern for %s is longer than allowed [%d] - check it!\n", ctemp1,BASE_OVERLAP+2);
                exit(1);
            } else {
                RE_rawsite = tok;    /* Recognition pattern, prev -> RE_rawsite */
            }
            tok = strsep(&rstr, " \t");
            while (tok[0] == '\0') tok = strsep(&rstr, " \t");
            itemp2 = atoi(tok); /* next one - bottom offset, prev -> itemp2 */

            /* these are the std ones; now need to test for the next tacg-specific ones */
            z = 0;
            while (z < 5 && rstr != 0L) {  /* this stanza grabs the remaining vars */
                tok = strsep(&rstr, " \t");
                while (tok[0] == '\0') tok = strsep(&rstr, " \t");
                if (tok[0] != '?' && tok[0] != '!') {
                    r[z] = atoi(tok);  /* the 0th el would be ERR */
                    z++;
                } else if (z == 0) {    /* we've hit the end - set to the defaults, break it out of it */
                    r[0] = r[3] = r[4] = 0;
                    r[1] = r[2] = 100;
                    z = 6;
                } else { /* anything over 1 indicates either a truncated or erroneous dam/dcm / cost data */
                    r[3] = r[4] = 0;
                    r[1] = r[2] = 100;
                    break;
                }
            }
            /* else everything's peachy */
            /* for r[], 0 = err; 1 = dam; 2 = dcm; 3 = price; 4 = #units */
            Err = r[0];
            if (r[4] != 0) {
                cost = r[3]/r[4];  /* added 9.10.99 for cost filtering */
            } else {
                cost = 0;
            }
            /* cost is in (integer) units/$, not $/unit, to allow more intuitive entering of the option
               '--cost 300' means you want to allow cutting only with REs that are >= 300 units/$
               bc it's integer, fractionals won't register, but if you care to go lower than 1 unit/$
               you probably want to go for all of them anyway */

            /* This filters the entire recognition site from rebase.data to the REstruct.  Difference
               from below is that it does it for the whole sequence, not just the 1st 6 bases for the hashtable */
            /* grab the raw site for labelling purposes */
            /* if no '-r' or ('-r' && current name hashes to a spot that one of the SelREs hashed to) */
            SEi = -1;    /* SE index */
            /* next line matches the name just read in (ctemp1) with the names entered in SelEnz[].PName.  Returns
               the index of the SelEnz match if one is found.  For -P, this match has to be further matched
               to the PP struct entry and FURTHER matched as to RE index (calc'ed here) so that we can fill
               out the PP[].REi[] value for Proximity() */
            strcpy(ctemp2,ctemp1);
            if (SelREs == 1) {
                SEi = MatchSelREs(ctemp2, SelEnz, NumSelREs);
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
                rsl = strlen(RE_rawsite)+1; /* rsl = Raw Site Length  */
                RE[*NumREs].E_raw_sit = (char *) calloc(rsl, sizeof(char)); /* 1st grab some space */
                if (RE[*NumREs].E_raw_sit == NULL) BadMem("calloc - RE[*NumREs].E_raw_sit", 1);

                strcpy(RE[*NumREs].E_raw_sit, RE_rawsite);  /* then copy it to RE_rawsite */
                /* Stanza to calc overhangs for matching to filter-by-overhang */
                /* calc overlap value for exclusion clause below */
                if (F.Overhang != 1) {    /* if we want to restrict selection by overhang */
                    OLapLen =  itemp2;
                }
                i=len=0;
                while (RE_rawsite[i] != '\000') {  /* while not at end of string  */
                    if (RE_rawsite[i] != '_' && RE_rawsite[i] != '\'') { /*As long as it's a valid char */
                        RE_rawsite[len++] = tolower(RE_rawsite[i++]); /* The struct value gets the whole site */
                    } else i++;  /* incr i only, skipping the _ and ' chars  */
                }
                /* Here's where the 'nnn' detection/handling goes to trim leading & trailing n's */
                //this should be made a fn() and applied to the degen seqs generated from the error term ..

                if (RE_rawsite[0] == 'n') {
                    lnct=0; //set leading n counter
                    while (RE_rawsite[lnct] == 'n') {
                        lnct++;
                    }
                    itemp1 -= lnct;
                    itemp2 -= lnct;
                    i = 0; //reset this counter
                    while (lnct < len) {
                        RE_rawsite[i++] = RE_rawsite[lnct++];    // shift the pattern over
                    }
                }
                while (RE_rawsite[--len]=='n'); /* back up until hit a non 'n' */
                RE_rawsite[++len] = '\0';     /* then restore the end of string mark */




                /* Now calculate E_mag exactly to see if it should be included in the digestion */
                fRE_mag = MagCalc(RE_rawsite, len); /* calculate the 'magnitude' of the pattern */
                if (fRE_mag<2) {
                    fprintf(stderr, "Pattern named %s is too degenerate: Magnitude is only %4f.\n",
                            ctemp1, fRE_mag);
                }
                if (fRE_mag < F.Mag && (F.Verbose >= 1 || F.Xplicit == 1 || F.Pat == 1)) { /* only if you really want the verbiage... */
                    fprintf(stderr, "\nWARNING: Pattern %s is too degenerate: Magnitude is only %4f.\n",ctemp1, fRE_mag);
                }

                /* don't bother loading any REs that don't match these criteria - may have to add a check for
                   dam/dcm sensitivity as well as COST here as well */

                if (((fRE_mag >= (float)F.Mag) &&
                        (cost >= F.Cost)) &&

                        ( (F.Overhang == 1)  ||  /* either we want all of them or we have to go thru the possibilities */
                          (((F.Overhang == 5)  && (F.Overhanglen == 0) &&  (OLapLen > 0)) || /* if olap < 0 && overhang = 5' w/ no overhang length constraints */
                           ((F.Overhang == 5)  && (F.Overhanglen  > 0) &&  (OLapLen == F.Overhanglen)) || /* ditto above but WITH ohlen constraints */
                           ((F.Overhang == 3)  && (F.Overhanglen == 0) &&  (OLapLen < 0)) || /* if olap > 0 && overhang = 3' */
                           ((F.Overhang == 3)  && (F.Overhanglen  > 0) &&  (OLapLen == (F.Overhanglen * -1))) || /* ditto above but WITH ohlen constraints */
                           ((F.Overhang == 0)  &&                           OLapLen ==0)))) {

                    RE[*NumREs].E_mag = (int) fRE_mag; /* this truncates fRE_mag, but since it's upwardly inclusive ... */
                    RE[*NumREs].E_wsit[ori] = (char *) calloc(len+1, sizeof(char));
                    if (RE[*NumREs].E_wsit[ori] == NULL) BadMem("calloc - RE[*NumREs].E_wsit[ori]", 1);
                    strncpy(RE[*NumREs].E_wsit[ori], RE_rawsite, len); /* E_wsit = the recog site, stripped of _ and ' and extra n's */
                    /* need a \0 to terminate ? - shouldn't - it's been calloc'ed to 0's*/
                    RE[*NumREs].E_len = len;  /* set the real recog length since we calculated it */
                    RE[*NumREs].E_pal = palindrome (RE[*NumREs].E_wsit[ori], len); /* and determine if the whole site is a pal */


                    /* for Proximity() - Do all the matching for relating RE names in SelEnz and PP here,
                       as long as we have all the necessary vars to do it. */
                    if (F.Prox > 0) { /* if we need Proximity matching */
                        PPi = (int)SEi/2; /* the index is half of SEi, as there are 2 names per PP entry */
                        if (((SEi%2)-0.5) < 0.1) PP[PPi].REi[0] = *NumREs; /* which one of the 2 depends on the mod */
                        else PP[PPi].REi[1] = *NumREs;
                    }

                    /* Filter the raw site from rebase.data to the hexamer val in struct
                          that's submitted for the numeric conversion */
                    /* And assign the rest of the temps to the struct vars */
                    RE[*NumREs].E_nam_l = strlen(ctemp1); /* get the length of the RE name and pop it in */
                    strcpy (RE[*NumREs].E_nam, ctemp1);    /* and pop the RE name in too */
                    RE[*NumREs].E_tcut[ori] = itemp1; /* this is a temp val - corrected in BestHexWWPal() */
                    if (F.Verbose > 1) fprintf(stderr, "\nINFO:ReadEnz:for RE=%s, ori=%d, itemp1=%d, RE[*NumREs].E_tcut[ori]=%d\n", RE[*NumREs].E_nam, ori, itemp1, RE[*NumREs].E_tcut[ori]);
                    RE[*NumREs].E_olap = itemp2;
                    RE[*NumREs].Err = Err;



                    /* these have to be set to defaults if SelEnz is not used to set min/Max */
                    if (F.SelEnzUsed == 1) {                  /* added 8.19.99 for SlidWin */
                        RE[*NumREs].min = SelEnz[SEi].min;
                        RE[*NumREs].Max = SelEnz[SEi].Max;
                    } else {
                        RE[*NumREs].min = 0;
                        RE[*NumREs].Max = 32000;
                    }

                    RE[*NumREs].proto = *NumREs; /* so that it can point back to itself in the Cutting code */
                    RE[*NumREs].dam = r[1];  /* added 9.10.99 for dam methylation analysis */
                    RE[*NumREs].dcm = r[2];  /* added 9.10.99 for dcm methylation analysis */
                    RE[*NumREs].Cost = cost; /* added 9.10.99 for cost filtering */

                    Nprotos++; /* count the protos so we can return the number at end */
                    i = *NumREs;
                    /* this fn() does a lot of the calc'n / assignments common to normal and ERRORs code */
                    BestHexWWPal(RE_rawsite, len, i);
                    /* Now do the ugly handling of nonpalindromes - this handles non-pal's as another entry in
                       the SAME RE entry and only adds those entries that are different - big savings in space  */
                    if (RE[*NumREs].E_pal == 0)  { /* if the RE is not a pal, do the rest  */
                        /* reverse complement and pop in the pattern sequence */
                        RE[*NumREs].E_wsit[dpl] = (char *) calloc(RE[*NumREs].E_len+1, sizeof(char));
                        if (RE[*NumREs].E_wsit[dpl] == NULL) BadMem("calloc - RE[*NumREs].E_wsit", 1);
                        Anti_Par(RE[*NumREs].E_wsit[ori], RE[*NumREs].E_wsit[dpl], RE[*NumREs].E_len);
                    }  /* end nonpalindrome handling */
                    if (F.Verbose > 1) {
                        fprintf(stderr, "%s: %s -> %s: pal=%d tcut=%d olap=%d \n", RE[*NumREs].E_nam, RE[*NumREs].E_raw_sit, RE[*NumREs].E_wsit[1], RE[*NumREs].E_pal, RE[*NumREs].E_tcut[ori], RE[*NumREs].E_olap);
                    }
                    /* Here's where we figure out if we need to do the error matching and if so, create
                       and check the sequences that have to be checked.  Obviously, this is meant to be
                       for a very few sequences at once; otherwise the memory requirements explode */
                    if (Err > 0 && Err <= MAX_ERR) { /* if there is an ERROR term  */
                        /* now need to generate the list of names and validate that they're novel */
                        BOStanza = proto = *NumREs;   /* points to the prototype RE entry */
                        EOStanza = here = *NumREs+1;  /* pointer for the current RE entry */

                        for (E=0; E<Err; E++) {    /* 1 iteration for each error level  */
                            if (F.Verbose > 1) fprintf(stderr, "\tPrototype Expansion, Cycle %d on %s\n", E, RE[proto].E_wsit[ori]);
                            for (z=BOStanza; z<EOStanza; z++) { /* for Err=1, goes from 'proto' to 'here' */
                                /* next 2 'for's are the core of the ERROR handling - build the ERROR templates  */
                                /* 'n's at both ends effectively decr length by one - addressed below */
                                len = RE[z].E_len ;  /* this may not be nec here, but if this code is fn()ized, it may be */
                                for (j=0; j<len; j++) { /* for the length of the current (NOT nec the proto) pattern  */
                                    RE[here].E_wsit[ori] = (char *) calloc(RE[proto].E_len+1, sizeof(char));
                                    if (RE[here].E_wsit[ori] == NULL) BadMem("calloc - RE[here].E_wsit[ori]", 1);

                                    /* make an 'n' scan of the sequence - but for beginning and end just knock off the char- we know it'll match anything there; why test for it?  */
                                    // so do these seq edits 1st
//									for (i=0; i<(len-1);i++) { //shift over sequence 1st
                                    // fprintf(stderr, "j=%d, for here=%d, char i= %d = %c\n", j, here, i, RE[z].E_wsit[ori][i]);
//										RE[here].E_wsit[ori][i] = RE[z].E_wsit[ori][i+1]; //move sequence to beginnning
//									}
//									fprintf(stderr, "trimmed, moved seq: %s\n", RE[here].E_wsit[ori]);

                                    for (i=0; i<(len); i++) { /* for the length of the pattern */
                                        if (i==j ) RE[here].E_wsit[ori][i] = 'n';
                                        else RE[here].E_wsit[ori][i] = RE[z].E_wsit[ori][i];
//                              if (i>0 && i < (len) && i==j ) RE[here].E_wsit[ori][i] = 'n';
//                              else RE[here].E_wsit[ori][i] = RE[z].E_wsit[ori][i];
                                    }
//									fprintf(stderr, "n-scanned seq: %s\n", RE[here].E_wsit[ori]);

                                    /* now have a new, modified pattern - 1st trim n's from both ends and then ask is it novel? */
                                    RE[here].E_len = strlen(RE[here].E_wsit[ori]);  /*  the default len, unless otherwise mod'ed */

                                    /*  1st decide if it's a pal, then do both at once, rather than do the ori
                                       then the doppel */
                                    if (palindrome(RE[here].E_wsit[ori], RE[here].E_len) == 1) { /* if it's a pal.. */
                                        RE[here].E_pal = 1; /* just mark it - it will match the previous ones either way */
                                    } else { /* but if it's not a pal */
                                        RE[here].E_pal = 0;  /* also mark it as such and flip it to the AP in the dpl place */
                                        RE[here].E_wsit[dpl] = (char *) calloc(RE[here].E_len+1, sizeof(char));
                                        if (RE[here].E_wsit[dpl] == NULL) BadMem("calloc - RE[here].E_wsit[dpl]", 1);
                                        Anti_Par(RE[here].E_wsit[ori], RE[here].E_wsit[dpl], RE[here].E_len);
                                    }
                                    same = 0;
                                    /* and now check from where we STARTED on this stanza to where we ARE unless the seq is identical to
                                       a previous one along the way - this way does both ori/dpl for both the current pattern as well
                                       as the previous ones at the same time */
                                    /* THIS LOOP IS THE HUGE TIME SINK ON MULTIPLE ERRORS; maybe try a hash calculation to speed things up? */
                                    /* or use Knight's grepseq code to speed things up  */
                                    for (m=proto; m<here && same == 0; m++) { /* match vs the ori AND doppel at each step with..*/
                                        if (RE[here].E_len == RE[m].E_len) {
                                            p=1;
                                            if (RE[here].E_pal == 0) p=0; /* if !pal, have to do it 2x */
                                            for (n = 1; n>=p  && same==0; n--) {
                                                rep=1;
                                                if (RE[m].E_pal == 0) {
                                                    rep = 0;
                                                }

                                                for (l=1; l>=rep  && same==0; l--) {  /*  ..the newly generated seq */
                                                    same = Degen_Cmp(RE[m].E_wsit[l], RE[here].E_wsit[n], RE[here].E_len, 2, 1);
                                                }
                                            }
                                        } else same = 0; /* if they're ! the same size they can't be identical */
                                    }
                                    if (same == 0) { /* then it's scanned all the relevant RE bits and it really is novel.. */
                                        /* so fill out the rest of the bits that need to be filled out */
                                        BestHexWWPal(RE[here].E_wsit[ori], RE[here].E_len, here);
                                        RE[here].proto = proto;  /* set to credit the proto for cutting/labelling */

                                        if (F.Verbose > 2) fprintf(stderr, "for proto %s, degen seq %d: %s\n", RE[proto].E_wsit[ori], here, RE[here].E_wsit[ori]);
                                        here++;
                                    }  /* the above BestHexWWPal() should take care of the rest of the required values */
                                    if (EORE - (here + *NumREs) < 12) {
                                        if (F.Verbose > 2) fprintf(stderr, "\nNeed mem for RE (prototype expansion)..");
                                        EORE += RE_INCR;
                                        RE = realloc(RE, sizeof(*RE)*EORE);
                                        if (RE == NULL) BadMem("Failed to get more mem for RE.\n", 1);
                                        if (F.Verbose > 1) fprintf(stderr, "realloc success - RE NOW contains %d elements\n", EORE);
                                    }
                                }  /* for (j=0;j<len; j++) ... */
                            }  /* for (z=BOStanza; z<EOStanza; z++) */
                            BOStanza = EOStanza; /* and reset the Stanza markers to new values */
                            EOStanza = here;
                        }  /* for (E=0; E<Err;E++) ... */
                        if (F.Verbose > 1) {
                            fprintf(stderr, "RE NOW contains %d elements\n", (here-1));
                            for (i=1; i<here; i++) {
                                fprintf(stderr, "seq %d: %s\n", i, RE[i].E_wsit[ori]);
                            }
                        }
                        *NumREs = here; /* and sync the two - this handles all the incr of *NumREs in ERROR mode */
                    }  /* if (Err > 0 && Err < 4) ... */

                    if (Err == 0) (*NumREs)++;   /* if we're not handling ERRORs, incr it */
                    /* BUT if the RE is sensitive to dam or dcm and there's a --dam or --dcm flag, then don't   */
                    /* this might have to be changed to include the RE, but not allow any cutting */
                    /* if ((F.dam == 1 && RE[*NumREs - 1].dam == 101) || (F.dcm == 1 && RE[*NumREs - 1].dcm == 101)) (*NumREs)--;  */

                    if (EORE - *NumREs < 12) {
                        if (F.Verbose > 2) fprintf(stderr, "\n\t\t#REs=%d - Need mem for RE..", *NumREs);
                        EORE += RE_INCR;
                        RE = realloc(RE, sizeof(*RE)*EORE);
                        if (RE == NULL) BadMem("Failed to get more mem for RE.\n", 1);
                        if (F.Verbose > 2) fprintf(stderr, "..Got it! RE NOW contains %d elements\n", EORE);
                    }
                }  /* end of magnitude check */
            }  /* if (SelREs != 1 || SEi != -1) ...*/
        } /* end of if statement that checks for ';' at beginning of RE name */
        rstr = tmp;  /* IMPORTANT - re-set the pointer to the ORIGINAL start of rstr */
    } /* End of while loop that reads in  RE enzyme data into struct ... whew! */

    if (SelREs ==1) { /*                 used in both -r and -P flags; adjust the checking above to extend past 1st check */

        j = MAX_SEL_ENZ;
        for (i=0; i< j; i++) {   /* check for RE names used in -r flag but not matched in REBASE used */
            if ((SelEnz[i].match != 1) && (SelEnz[i].PName[0] != '\0')) { /* '*' indicates a match */
                fprintf(stderr, "Missing or Misspelled Pattern name: %s\n", SelEnz[i].PName);
                exit(1);
            }
        }
    }

    /* Now need to cycle thru the struct, generating the hashtable and chk_sits arrays */

    /* 1st, calc how many (NOT which - that's later) RE entries hash to a particular value in thrash/hashtable */
    for (m=1; m<*NumREs; m++) { /* start at 1, NOT F.RE_ST, NOT 0, fool! */
        rep = 1;
        if (RE[m].E_pal == 0) rep = 0;
        for (k=1; k>=rep; k--) { /* ugh - count down to catch the doppel */
            RE[m].E_dgen = hash(RE[m].E_hex[k], dgn_sits, 6); /* call to hash !!! */
            for (i=0; i<RE[m].E_dgen; i++) { /*  increment the thrashtable to reflect this degeneracy  */
                thrashtable[dgn_sits[i]]++;   /* incr the corresponding element of the array */
            }
        }
    }

    /* Now calloc the space that chk_sits (now incorp into hashtable) needs using the data from above  */
    for (i=0; i<4096; i++) {
        /* hashtable[4096][how many REs hash to this value + 1] */
        hashtable[i] = (int *) calloc((thrashtable[i]+1), sizeof(int));
        if (hashtable[i] == NULL) BadMem("ext mem for hashtable", 1);
        hashtable[i][0] = thrashtable[i];  /* and copy over the 1st element */
    }

    /* And finally rehash the enz list, filling out the hashtable array with the index of the RE struct
       that has to be checked in the event that the sequence hashes to the chk_sits index - ugly to do
       something twice like this, but it works pretty damn fast.  Using trashtable[1] as the pointer to keep
       track of where the next available space in hashtable is, rather than another element in hashtable itself
       as we already have that info in hashtable[0] */

    /* For next line, DON'T be tempted to use memset() - memset sets it on a per-byte basis, so for an
       int, each of the 4 or 2  bytes will be set to 1, yielding not 1 but some much larger number */
    for (j=0; j<4096; j++) thrashtable[j] = 1; /* set it to 1 NOT 0; have to skip 0th el of hashtable */

    for (i=1; i<*NumREs; i++) { /* start at 1, not F.RE_ST, NOT 0, fool! - have to take dam and dcm cm if needed */
        rep = 1;
        if (RE[i].E_pal == 0) rep = 0;
        for (k=1; k>=rep; k--) { /* ugh - count down to catch the doppel */
            Degen = hash(RE[i].E_hex[k], dgn_sits, 6);
            m = 1;
            if (k == 0) m = -1; /* set the multiplier for i below so that the doppel will be negated in
                                       the hashtable, so we can tell later on what caused the hit */
            for (j=0; j<Degen; j++) hashtable[dgn_sits[j]] [thrashtable[dgn_sits[j]]++] = i*m;
        }
    }
    /* the above results in hashtable[4096][....]  where the [0] element is the # of possible hits at
    that hashvalue (the value in the original hashtable, when it was a 1D arrary) and the [1-#]
    values are the indices of the RE entries that hash to this value (that in later exploits have to
    be examined to discover whether they are true hits or not).  This new version of hashtable
    eliminates chk_sits as well as shortening the allocation by one number ..wowee.  Not so much a
    great space saving, but a simplification  */

    free(ctemp1);
    free(ctemp2);/*  */
    fclose(tmpfp);
    if (F.Verbose > 0) {
        fprintf(stderr,"\nReadEnz: Nprotos = %d\n",Nprotos );
        fprintf(stderr,"ReadEnz: sizeof(RE) = %ld; Total primary alloc'ed size of RE = %ld\n", sizeof(RE), (Nprotos * sizeof(*RE)));
    }
    return Nprotos;
}


/* CalcEstSites() calculates the estimated number of hists that there should be based on the
   distribution of ACGT in the input sequence; returns the estimated # of hits depending on
   the entire sequence.
   NOW Does it for the site in BOTH strands (fixed 4.18.99) by Anti_Par'ing the
   site and running both sites thru the switch().
*/

float CalcEstSites(char *BindSite, int BSLen, float FrACGT[4])
{
    int i,j;
    char *BS[2];
    float P_Sites[2]; /* probability of the site occuring by chance in that sequence */
    if ((BS[0] = calloc(BSLen+1, sizeof(char))) == NULL) BadMem("No mem for BS", 1);
    if ((BS[1] = calloc(BSLen+1, sizeof(char))) == NULL) BadMem("No mem for BS", 1);
    strncpy(BS[0], BindSite, BSLen); /* copy in the orig site */
    Anti_Par(BS[0], BS[1], (long)BSLen); /* now AP it */
    P_Sites[0] = P_Sites[1] = 1.0;
    for (j=0; j<2; j++) {
        for (i=0; i< BSLen; i++) {
            switch (BS[j][i]) {
            default:
                break; /* n or anything else obnoxious */
            case 'a':
                P_Sites[j] *= FrACGT[0];
                break;
            case 'c':
                P_Sites[j] *= FrACGT[1];
                break;
            case 'g':
                P_Sites[j] *= FrACGT[2];
                break;
            case 't':
                P_Sites[j] *= FrACGT[3];
                break;
            case 'w':
                P_Sites[j] *= (FrACGT[0]+FrACGT[3]);
                break;  /* a or t */
            case 's':
                P_Sites[j] *= (FrACGT[1]+FrACGT[2]);
                break;  /* c or g */
            case 'r':
                P_Sites[j] *= (FrACGT[0]+FrACGT[2]);
                break;  /* a or g */
            case 'y':
                P_Sites[j] *= (FrACGT[1]+FrACGT[3]);
                break;  /* c or t */
            case 'm':
                P_Sites[j] *= (FrACGT[0]+FrACGT[1]);
                break;  /* a or c */
            case 'k':
                P_Sites[j] *= (FrACGT[2]+FrACGT[3]);
                break;  /* g or t */
            case 'b':
                P_Sites[j] *= (FrACGT[1]+FrACGT[2]+FrACGT[3]);
                break;  /* not a */
            case 'd':
                P_Sites[j] *= (FrACGT[0]+FrACGT[2]+FrACGT[3]);
                break;  /* not c */
            case 'h':
                P_Sites[j] *= (FrACGT[0]+FrACGT[1]+FrACGT[3]);
                break;  /* not g */
            case 'v':
                P_Sites[j] *= (FrACGT[0]+FrACGT[1]+FrACGT[2]);
                break;  /* not t */
            }
        }
    }
    free(BS[0]);
    free(BS[1]);
    return (P_Sites[0] + P_Sites[1]);
}


