/* tacg - a command line tool for the analysis of nucleic acids and protein  */
/* Copyright � 1994-2005 Harry J Mangalam, tacg Informatics
(hjm@tacgi.com, 949 856 2847) */

/* $Id: Cutting.c,v 1.3 2005/01/26 22:44:08 mangalam Exp $  */

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

/* Cutting() is a functionization of the matching routines needed to do the cutting analysis on
   degenerate (gyrwsctkmyycg) sequence - also uses HorribleAccounting() as a function.  This fn()
   includes routines for both 'lite' and 'heavy' degenerate cutting, where lite is defined as ignoring
   any key hexamer which has ANY degeneracy and 'heavy' will force degenerate matching of any degeneracies
   up to a compiled-in limit (currently 256 = 4 ns (or equivs) out of the keying hexamer)
*/

/* Remember:  it is NOW this way....
if F.Degen = 0  FORCEs exclusion of IUPAC and other characters in sequence)
                  (forced default = Degen_Cmp mode 0)
if F.Degen = 1  default - cut as nondegenerate unless IUPAC characters found; then cut as '-D3'
                  (default = Degen_Cmp mode 0)
if F.Degen = 2  allow IUPAC characters; ignore them in the KEY hexamer, but match outside of KEY
                  (Degen_Cmp mode 2) - same as before
if F.Degen = 3  allow IUPAC characters; find only exact matches
                  (Degen_Cmp mode 2) - same as before
if F.Degen = 4  allow IUPAC characters; find ALL POSSIBLE matches
                  (Degen_Cmp mode 1)
*/

/* Cutting the Sequence, using SHKey/shortcut approach.....  The shortcut is to hash the
   sequence hexamer only once and then use the calculated values of the hexamer to assist in
   calculating the next hash key.  Since 5 of the 6 characters have already been 'hashed' for
   the next key, why recalculate them and suffer the overhead of a function call, when you can
   (with just a little extra code) get the same value on the fly, inline....  It appears to be
   about 5-7x times faster for long sequences (where the cost is in the rehashing, not the base
   cost of hashing the RE data twice.  Surprisingly enough, this is the same approach (albeit
   implemented less elegantly) used in Udi Manber's latest approximate string matching, as well
   as BLAST.   */

/* Cutting() has now been improved in to adaptively switch between the shortcut approach
   and the hash() approach to yield speeds for degenerate sequence that are almost as fast as
   those for nondegenerate sequence.  */


/* Because of the way SHKey is calculated, sequence indices should start on a '1' boundary, rather than a '0'
   boundary so that a 1 can be subtracted without catastrophe - also makes it easier to understand where the
   sequence actually is */

void Cutting(char *sequence, int *hashtable[4096])
{
    int i,  degen, mode=2;
    int Key=0, lead=0, dgn_sits[256], lag=0;
    char RE_hex[6];
    long CSP=0, lj, progress=0, onep=0,
         DCnt=7; /* DCnt has to be 'long' to be able to be decr to at least -1,000,000,000
                     also started at 7 to force careful checking of 1st hexamer JIC */

    for (i=0; i<6; i++) RE_hex[i] = sequence[DD->BOS-1+i]; /* then 5 more of the real seq, since it's discarded anyway */
    degen = Degen_Calc(RE_hex); /* being careful cuz it might start with degenerates */
    if (degen <= MAX_DEGEN) {
        degen = hash(RE_hex, dgn_sits, 6);  /* calculate the initial hash value */
    }
    Key = dgn_sits[0];

    /* Below sets up the margins for running thru the seq - if it's circular (topo=0), have to include the
       repeated seqs on the beginning and end.  if topo=1, then just have to run thru the real seq (minus the
       repeats at begin end).  In both cases, have to go thru DD->Dat and 'normalize' the position info to the
       actual sequence coordinates.  by subtracting BASE_OVERLAP from the index as written - altho this could be
       done on the fly without too much effort */

    if (F.Degen < 2)  mode = 0;  /* set the mode for all cases here at start */
    else if (F.Degen > 3) mode = 1; /* mode = 2 is set at beginning */

    if (F.Verbose > 0) {
        fprintf(stderr,"\nCutting.c:Cutting(): Begin 'cut'..\n" );
        onep = (long)(DD->EOS / 100);
    }

    for (CSP=DD->BOS; CSP<DD->EOS-4; CSP++) {   /* do the whole seq incr 1 nuc at a time */
        if (F.Verbose > 0 && CSP > progress) {
            fprintf(stderr,"\nCutting.c:Cutting(): progress: %ld\n", CSP );
            progress += onep; // incr onep to another %
            fprintf(stderr, "Matches so far = %d\n", (int)(DD->Dat[0]-2)/2);
        }

        lj = CSP+5;  /* CSP+5  lj points to the 'lagging' base */
        /* calculate 'lag' value 1st to check if there's incoming degenerates */
        switch (sequence[lj]) {    /* a=0, c=1, g=2, t=3,  */
        case 'a':
            lag = 0;
            DCnt--;
            break;   /* lag values are numeric value * 1024 */
        case 'c':
            lag = 1024;
            DCnt--;
            break;   /* ..and decr the degen Cnt - when it goes below 0 we can use the FAST method */
        case 'g':
            lag = 2048;
            DCnt--;
            break;
        case 't':
            lag = 3072;
            DCnt--;
            break;
        default:
            DCnt = 7; /* degen char detection; bump up the Degen Cnter to deny fast cutting for the window  */
            break;
        }
        if (DCnt < 1) {
            /* calculate 'lead' value */
            switch (sequence[CSP-1]) {     /* [CSP-1]    a=0, c=1, g=2, t=3,  */
            case 'a':
                lead = 0;
                break;   /* lead values are numeric value * 1 */
            case 'c':
                lead = 1;
                break;   /* doesn't matter if they're degen - already handled by the 1st stanza */
            case 'g':
                lead = 2;
                break;
            case 't':
                lead = 3;
                break;
            default:
                break;  /* bad character detection - but it's already been done */
            }
            /* instead of calling 'hash' on each new hexamer, can calc it incrementally as below */
            Key = ((Key-lead) >> 2) + lag;  /*  bit shift version - will it work?  yes! */
            /* fprintf(stderr," r %c %d  ", sequence[CSP-1], Key);   */
            /* now check the generated key against the hashtable to see if there's an entry */
            if (hashtable[Key][0] != 0) { /* no degeneracy so only have to check Key position (dgn_sits[0] if using hash() */
                HorribleAccounting(sequence, Key, mode, hashtable, CSP /*, already */);
            }
        } else { /* there's still a degen in the window, so have to calc degen via hash() */
            degen = Degen_Calc(sequence+CSP);    /*  && F.Degen > 2  */
            if (degen <= MAX_DEGEN) { /* if not, just ignore it */
                degen=hash(sequence+CSP, dgn_sits, 6);
                /* below stanza is for ALL cases of degeneracy in the window - mode sets the different treatments */
                for (i=0; i<degen; i++)  {
                    Key = dgn_sits[i]; /* clarity and coherence with other Fn()s */
                    /* fprintf(stderr," x %c %d  ", sequence[CSP-1], Key);   */
                    /*               if (F.Degen > 2 || CSP < (DD->BOS)+6) { */
                    if (hashtable[Key][0] != 0) { /* could be lots of degeneracy !!!  This is gonna suck cycles */
                        HorribleAccounting(sequence, Key, mode, hashtable, CSP);
                    }
                    /*                } */
                }
            } /* if (degen <= MAX_DEGEN) .... */
        } /* else   there's still a degen in the window... */
    } /* for (CSP=DD->BOS; CSP<DD->EOS-3; CSP++)... */

    DD->Dat[DD->Dat[0]] = -22222; /* mark the end of DD->Dat so we know where to end later on */
    if (F.Verbose >= 1) fprintf(stderr, "\nTotal # Matches = %d\n", (int)(DD->Dat[0]-2)/2);  /* if verbose, how many matches? */
}


/* HorribleAccounting() is a function that became one simply because it contains common code for the cutting
   functions mentioned above.  If it turns out to be about as fast as the inline code, I'll functionize it in
    NonDegenCut() as well */

void HorribleAccounting(char *sequence, unsigned int Key, int mode, int *hashtable[4096], long CSP /*, int already[ALRDY+1]*/)
{

    long rwc;
    int m, Cur_RE, Cur_REi, pal;
    for (m=1; m<=hashtable[Key][0]; m++) {
        Cur_RE = hashtable[Key][m]; /* For clarity's sake!! */
        Cur_REi = abs(Cur_RE);  /* Cur_RE used as the index to any arrays */

        /* Also writes REAL WORLD coords to DD->Dat, so Frag_Siz will be easier */
        /* Large # of recog sites are <= 6 base pairs so the last part of the following test should be
           executed only very rarely ... NOT - actually it is gotten to a surprising # of times.. */

        /* Bit below takes care of determining if, in the degenerate cases, if there are multiple Identical
           names at the same position.  Since the overlapping degeneracy can often generate identical sites from
           degenerate siblings.  Put behind the 'if' to save time when doing nondegenerate cutting (which
           will probably be most of the time  */

        /* only the real bozo degenerates (-d/4, -D/3) have to put up with this - this is skipped if we're doing
           either the default NonDegenCutting or -i/2 cutting, but if we have to do it, we're checking whether the
           current RE has already hit at this base (assuming that we're not at the end of the seq) */


        /* for the line below, both i and already[ALRDY] should always be 0 for F.Degen
           if it did find a hit, means another degeneracy allowed a previously noted enzyme to map to the
           same base, making 2 hits of the same enz at the same base, not really nec to note. */

        /* Is this ending condition OK - think so... */
        pal = 1;  /* if it's a palindrome; pal for clarity instead of RE[Cur_RE].E_pal */
        if (Cur_RE < 0) pal = 0; /* if Cur_RE is neg, then it's the doppel that got us here */
        // debug code to check if this will peek past the sequence buffer
        //if (sequence + CSP + RE[Cur_REi].WW[pal] < sequence){fprintf(stderr, ".\n"); }

        // 1st part of following test checks if offset can force a read before startof sequence
        if ((sequence + CSP + RE[Cur_REi].WW[pal] >= sequence) && (((RE[Cur_REi].E_len < 7) && (F.Degen < 3)) ||  ((CSP + RE[Cur_REi].E_len <= DD->EOS) &&
                (Degen_Cmp(RE[Cur_REi].E_wsit[pal], sequence+CSP, RE[Cur_REi].E_len, mode, RE[Cur_REi].WW[pal]) == 1)))) {
            /* if we made it here, it's a hit - mark the enzyme as already having cut at this position */
            /* using the tail posit'n of already[] as the index */
            /* rwc = real world coordinates; has to be calculated differently for circ and linear, but looks OK for both ERROR
              and normal code - calculated with the ~tcut and ~WW specific for each enz, even if credit is given to the proto */

            /* this is the bastard that calcs the correct rwc to enter in to DD.  */
            if (F.Verbose >1) {
                fprintf(stderr, "for pal=%d, Cur_REi=%d: CSP  -  BASE_OVERLAP  +  RE[Cur_REi].E_tcut[pal]  +  RE[Cur_REi].WW[pal]\n", pal, Cur_REi);
                fprintf(stderr, "                        %ld         %d                  %d                           %d\n", CSP, BASE_OVERLAP, RE[Cur_REi].E_tcut[pal], RE[Cur_REi].WW[pal]);
            }
            rwc = CSP-BASE_OVERLAP+RE[Cur_REi].E_tcut[pal]+RE[Cur_REi].WW[pal];  /* temp var to save adds below */
            if (F.Verbose >1) {
                fprintf(stderr, "RWC=%ld  Cur_REi = %d  DD.EOS=%ld  for: <%s  >%s\n", rwc, Cur_REi, DD->EOS, RE[Cur_REi].E_wsit[0], RE[Cur_REi].E_wsit[1]);
            }

            if (rwc > 0 && rwc < DD->EOS)  {   /* if it's within the sequence */
                /* Have to fiddle with what RE index gets the credit here - if it was one of the ERROR templates that got us here,
                have to give the credit to the prototype */
                //fprintf(stderr, "hit passed @ DD->Dat[0]=%ld\n", DD->Dat[0]);
                /* verbose debugging - when -> multiple levels of V, this one should be higher */


                RE[RE[Cur_REi].proto].E_Ncuts++;      /* use the proto index (both the proto and his ERROR progeny point to himself */
                DD->Dat[(DD->Dat[0])++] = rwc;   /* now should be set to real world coordinates - correct for ERROR as well */

                /* It's OK to just use the proto pointer -  because the nonpals are
                integrated into a single RE entry, BUT we NOW NEED the neg/positive indicators here
                to allow orientation to be indicated on the map; the 'rwc' has already been
                calc'ed and entered, so this should just be  the pointer to the name of the Enz
                Change from ignoring the orientation 1-16-99 */

                if (Cur_RE > 0) {
                    DD->Dat[(DD->Dat[0])++] = RE[Cur_REi].proto;
                } else {  /* negate the PATTERN index (not the site)  */
                    DD->Dat[(DD->Dat[0])++] = RE[Cur_REi].proto * -1;
                }
            }
            DD_Mem_Check();
        }
    }
}


/* lil nothing fn() to check that DD->Dat isn't going to blow */
int DD_Mem_Check(void)
{
    if (DD->Dat[1] - DD->Dat[0] < 10) {  /* Now check to see if we need more more memory */
        DD->Dat[1] = DD->Dat[1] + 5000;
        DD->Dat = (long *) realloc (DD->Dat, sizeof(long *) * (DD->Dat[1])); /* realloc the size 5000 bigger */
        if (DD->Dat== NULL) BadMem("Cutting.c: realloc failed (DD->Dat mem)\n",1);
    }
    return 0;
}


/****************************  Function Degen_Cmp  **********************************
   Degen_Cmp compares 2 strings, the 1st a degenerate sequence (typically a
   restriction enz or transcription factor binding site) and the 2nd a
   longer, target seq (degenerate (mode=0) or NON degenerate (mode=1)
   typically the sequence that is being analyzed.  The pointers passed to
   the function point to the initial base.
******************************************************************************************/
int  Degen_Cmp (char *d, char *s, int len, int mode, int WW)
{
    /* use mode 2 for the ERRORs fn() - exact matching on degenerate sequence */
    /* needs to be modified to be able to extend the search both ways if indicated by RE[].WW
       AND needs to rematch the initial key hexamer if F.Degen = 3 (-D flag) because the generating
       function for the degeneracies (hash()) will match like -d (find ALL degeneracies, NOT find only
       definite matches in degenerate pattern/target) */
    /* The indicators for the above things to be added are available already via RE[] and F */
    /* d     pointer to the degenerate sequence that is going to be used as a probe (1st char in the
             key hexamer.
       s     pointer to the longer NON/degenerate seq that is going to be used as the target or database
             into which the probe sequence might map.  THIS FUNCTION EXPECTS 's' TO ALREADY BE ADVANCED TO
             THE OFFSET ie PRECALCULATED!!
       len   length of the pattern seq or how many steps have to be checked
       mode  whether the target is degenerate or not
       WW    Which Way - the sequence has to be checked - backward or forward - depends on where in the Enz
             pattern the best/most significant hexamer is and is in normal cases, either 0 or negative.
             For the ERRORs, it's set to a positive number so as to do the comparison in sync */

    /* Remember:  it is NOW this way....
    if F.Degen = 0  FORCEs exclusion of IUPAC and other characters in sequence)
                      (Degen_Cmp mode 0)
    if F.Degen = 1  default - cut as nondegenerate unless IUPAC characters found; then cut as '-D3'
                      (Degen_Cmp mode 0)
    if F.Degen = 2  allow IUPAC characters; ignore them in the KEY hexamer, but match outside of KEY
                      (Degen_Cmp mode 2)
    if F.Degen = 3  allow IUPAC characters; find only exact matches
                      (Degen_Cmp mode 2)
    if F.Degen = 4  allow IUPAC characters; find ALL POSSIBLE matches
                      (Degen_Cmp mode 1)
    */

    int same = 1;  /* same = 1 if dgen seq can map onto the pure seq, 0 otherwise */
    int m,         /* counter for the pattern sequence */
        start=0,   /* default here is to check everything, from start to finish; required for the -D flag */
        t = WW;  /* negate sequence counter index to adjust for starting before the default start -
                this will always be either 0 or negative, but its sign is used as indicator */

    if (WW == 0) {    /*do it the normal way */
        start = 6;  /* start checking AFTER 1' hex, cuz that's how we got here.. */
        t = 6;      /* and sequence cmp starts in the same place as well */
    }
    /*    if (WW > 0) t = 0;  for the ERRORs bit, where RE patterns are matched against slight variants
                            so the paterns have to start in the same place */

    /* that's it - the default is to check the whole seq from start to end, base by base:
       pattern goes from       start = 0               to len-of-pattern
       sequence goes from      s + RE[CurE].WW (- #)   to len-of-pattern (counted by m)
    */
    /*  	fprintf(stderr, "-D %d ", F.Degen);  */
    if (mode == 0) {  /* mode = 0; degenerate probe, NONdegenerate sequence (default/original way)*/
        for (m=start; m<len && same==1; m++,t++) { /* start at 6 - we know that it is a match to 6 letters */
            switch (d[m]) { /* or we wouldn't be here as soon as 'same' == 0, it fails, so don't bother checking further */
            case 'a':
                if (s[t] != 'a') same = 0;
                break;
            case 'c':
                if (s[t] != 'c') same = 0;
                break;
            case 'g':
                if (s[t] != 'g') same = 0;
                break;
            case 't':
                if (s[t] != 't') same = 0;
                break;
            case 'y':
                if (s[t] != 'c' && s[t] != 't') same = 0;
                break;
            case 'r':
                if (s[t] != 'a' && s[t] != 'g') same = 0;
                break;
            case 'm':
                if (s[t] != 'a' && s[t] != 'c') same = 0;
                break;
            case 'k':
                if (s[t] != 'g' && s[t] != 't') same = 0;
                break;
            case 'w':
                if (s[t] != 'a' && s[t] != 't') same = 0;
                break;
            case 's':
                if (s[t] != 'c' && s[t] != 'g') same = 0;
                break;
            case 'b':
                if (s[t] == 'a') same = 0;
                break;
            case 'd':
                if (s[t] == 'c') same = 0;
                break;
            case 'h':
                if (s[t] == 'g') same = 0;
                break;
            case 'v':
                if (s[t] == 't') same = 0;
                break;
            case 'n':
                break;   /* n = anything */
            }
        }
    } else if (mode == 1) { /* mode = 1; degenerate probe, degenerate sequence (-D 4 = all possible hits) */
        for (m=start; m<len && same==1; m++,t++) { /* start at 6 - we know that it is a match to 6 letters */
            switch (d[m]) { /* or we wouldn't be here as soon as 'same' == 0, it fails, so don't bother checking further */
            case 'a':
                if (s[t]!='a' && s[t]!='w' && s[t]!='n' && s[t]!='r' && s[t]!='k' &&
                        s[t]!='d' && s[t]!='h' && s[t]!='v') same = 0;
                break;

            case 'c':
                if (s[t]!='c' && s[t]!='s' && s[t]!='n' && s[t]!='y' && s[t]!='m' &&
                        s[t]!='b' && s[t]!='h' && s[t]!='v') same = 0;
                break;

            case 'g':
                if (s[t]!='g' && s[t]!='s' && s[t]!='n' && s[t]!='r' && s[t]!='m' &&
                        s[t]!='b' && s[t]!='d' && s[t]!='v') same = 0;
                break;

            case 't':
                if (s[t]!='t' && s[t]!='w' && s[t]!='n' && s[t]!='y' && s[t]!='k' &&
                        s[t]!='b' && s[t]!='d' && s[t]!='h') same = 0;
                break;

            case 'y':
                if (s[t]=='g' || s[t] =='a') same = 0;
                break;
            case 'r':
                if (s[t]=='c' || s[t] =='t') same = 0;
                break;
            case 'm':
                if (s[t]=='g' || s[t] =='t') same = 0;
                break;
            case 'k':
                if (s[t]=='a' || s[t] =='c') same = 0;
                break;
            case 'w':
                if (s[t]=='g' || s[t] =='c') same = 0;
                break;
            case 's':
                if (s[t]=='a' || s[t] =='t') same = 0;
                break;

            case 'b':
                if (s[t] == 'a') same = 0;
                break;
            case 'd':
                if (s[t] == 'c') same = 0;
                break;
            case 'h':
                if (s[t] == 'g') same = 0;
                break;
            case 'v':
                if (s[t] == 't') same = 0;
                break;
            case 'n':
                break;   /* n = anything */
            }
        }
    } else if (mode == 2) { /* mode = 2; degenerate probe AND sequence, but exact match to each other for -D2 -D3 flag  */
        for (m=start; m<len && same==1; m++,t++) { /* start at 6 - we know that it is a match to 6 letters */
            switch (d[m]) { /* or we wouldn't be here as soon as 'same' == 0, it fails, so don't bother checking further */
            case 'a':
                if (s[t] != 'a') same = 0;
                break;
            case 'c':
                if (s[t] != 'c') same = 0;
                break;
            case 'g':
                if (s[t] != 'g') same = 0;
                break;
            case 't':
                if (s[t] != 't') same = 0;
                break;

            case 'y':
                switch (s[t]) {
                case 'y':
                case 't':
                case 'c':
                    break; /*these will match */
                default:
                    same=0;
                    break;
                } /* anything else breaks the 'sameness' */
                break;

            case 'r':
                switch (s[t]) {
                case 'r':
                case 'a':
                case 'g':
                    break; /*these will match */
                default:
                    same=0;
                    break;
                } /* anything else breaks the 'sameness' */
                break;

            case 'm':
                switch (s[t]) {
                case 'm':
                case 'c':
                case 'a':
                    break; /*these will match */
                default:
                    same=0;
                    break;
                } /* anything else breaks the 'sameness' */
                break;

            case 'k':
                switch (s[t]) {
                case 'k':
                case 't':
                case 'g':
                    break; /*these will match */
                default:
                    same=0;
                    break;
                } /* anything else breaks the 'sameness' */
                break;

            case 'w':
                switch (s[t]) {
                case 'w':
                case 'a':
                case 't':
                    break; /*these will match */
                default:
                    same=0;
                    break;
                } /* anything else breaks the 'sameness' */
                break;

            case 's':
                switch (s[t]) {
                case 's':
                case 'g':
                case 'c':
                    break; /*these will match */
                default:
                    same=0;
                    break;
                } /* anything else breaks the 'sameness' */
                break;

            case 'b':
                switch (s[t]) {
                case 'b':
                case 'c':
                case 'g':
                case 't':
                case 'k':
                case 's':
                case 'y':
                    break; /*these will match */
                default:
                    same=0;
                    break;
                } /* anything else breaks the 'sameness' */
                break;

            case 'd':
                switch (s[t]) { /* ditto */
                case 'd':
                case 'a':
                case 'g':
                case 't':
                case 'k':
                case 'w':
                case 'r':
                    break;
                default:
                    same=0;
                    break;
                } /* ditto */
                break;

            case 'h':
                switch (s[t]) { /* ditto */
                case 'h':
                case 'a':
                case 'c':
                case 't':
                case 'm':
                case 'w':
                case 'y':
                    break;
                default:
                    same=0;
                    break;
                } /* ditto */
                break;

            case 'v':
                switch (s[t]) { /* ditto */
                case 'v':
                case 'a':
                case 'c':
                case 'g':
                case 'm':
                case 's':
                case 'r':
                    break;
                default:
                    same=0;
                    break;
                } /* ditto */
                break;

            case 'n':
                break; /* n can still = anything */
            }
        }
    }
    return same; /* if same = 1, the degen maps onto the pure seq */
}  /* if same = 0 , the degen cannot map onto the pure seq - ie they're different */
