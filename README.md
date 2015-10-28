tacg: a grep for Nucleic acids.
=================================
by Harry Mangalam <harry.mangalam@uci.edu>
v4.60, August 20th, 2012

tacg is an Open Source compiled C application for searching for patterns in 
nucleic acids. (see <http://www.biomedcentral.com/1471-2105/3/8>)

You can think of it sort of as the guts of the beautiful Macintosh 'DNA Strider'
but with text output (tho it can also produce circular plasmid maps).
It can search for:

- explicit patterns of nucleotides
- generate restriction digests and sort fragments by in different ways.
- patterns using IUPAC degeneracies
- regular expressions
- complex rules of patterns
- can tolerate degeneracies in the pattern to match as well as the DNA target.
- can translate DNA into protein in any combination of frame via several translation tables
- in both forward and reverse strands
- and many other functions.


It comes with man page, and is built via the standard './configure; make; make install' and
generally builds on modern Linux distro's and can usually be convinced to compile
on Macs.  Windows?  Yeah, right....
It also comes with the makings of a rudimentary (and badly color-coordinated) web
page & cgi (see the 'tacgi4' subdir).

Once built, 'tacg -h' should start the help system.

If it doesn't build or doesn't behave, please give me a shout.

Harry Mangalam <hjmangalam@gmail.com>



