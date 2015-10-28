#!/usr/bin/perl -w
# Please skip down to line ~35 to read the config comments, set them appropriately,
# then SET $CONFIGURED to 1 and run the script.

use Getopt::Long;
use strict;

use vars qw($WWWHOST $WEB_DOC_ROOT $WEB_CGI_ROOT $TACGLIB $TACGBIN
$CONTACT_EMAIL %SUBST $WWWUSER $WWWGROUP $tmp $CGITMPDIR
$SYSTMP $SYSTMP_DEF $TACGLIB $TACGLIB_DEF $TACG_DIR_DEF $TACG_DIR @Files
$HTMLDIR $CGI_DIR $THIS_DIR $USER $HTML_LINK $CGI_LINK
);
# $file $input_file

# following are options if you wanted to set everything from the command
# line in an installation script
&GetOptions(
   "wwwhost=s"        => \$WWWHOST,       # fully qualified domain name
   "htmldir=s"        => \$HTMLDIR,       # DocumentRoot
   "cgidir=s"         => \$CGI_DIR,       # Scriptalias
   "cgitmpdir=s"      => \$CGITMPDIR,     # CGI tmp dir (linked to systmp)
   "tacglib=s"        => \$TACGLIB,       # the TACGLIB env variable
   "tacgbin=s"        => \$TACGBIN,       # the full path name to tacg binary
   "contact_email=s"  => \$CONTACT_EMAIL, # person to email in case of trouble
   "wwwuser=s"        => \$WWWUSER,       # the www user id
   "wwwgroup=s"       => \$WWWGROUP,      # the www group id
   "systmp=s"         => \$SYSTMP,        # where the system /tmp is
   "tacgdir=s"        => \$TACG_DIR,      # what the name of the tacg-related
                                          #   dirs will be (under tmp, html, etc
);
$THIS_DIR = `pwd`;
chomp $THIS_DIR;

###########################################################################
#        Edit the following variables to suit your machine/setup          #
#       then set $CONFIGURED to 1 (below)  and then re-run the script     #
#       The values shown are examples only, of course; please use them    #
#       as templates to set the quoting & slashing correctly              #
#                                                                         #
#       !!!  This is an alpha version!!!                                  #
#       !!!  Please let me know how it fails on your system  !!!          #
###########################################################################
#   The variables below are configured for a vanilla Debian system with   #
#   user hjm, and group hjm, using standard installation dirs.  On RH and #
#   other *nix systems, you'll have to correct them yourself.  There      #
#   should be a std web configurator and there probably is, but I haven't #
#   found it.                                                             #
###########################################################################

# NB:  if configuring to have this installed in your own public_html dir,
# please use fully qualified paths:
# ie: /home/you/public_html, not ~you/public_html

# Also, You must allow your webserver to allow FollowSymLinks.

$WWWHOST   = $SUBST{WWWHOST}   = "localhost";  # your fully qualified hostname.domain (and :port, if not port 80).
$HTMLDIR   = $SUBST{HTMLDIR}   = "/var/www";     # the DocumentRoot dir from your httpd config
$CGI_DIR   = $SUBST{CGI_DIR}   = "/usr/lib/cgi-bin";    # the ScriptAlias dir from your httpd config
$TACG_DIR  = $SUBST{TACG_DIR}  = "tacg4";    # the name of the subdir that will be created in both
                                               # the system tmp dir & the DocumentRoot dir to
                                               # consolidate the tacg-related files
   # ie /tmp/tacgdir & /home/www/htdocs/tacgdir, if the settings shown were used

$TACGBIN   = $SUBST{TACGBIN}   = "/usr/local/bin/tacg"; # the tacg executable you wish to use
$TACGLIB   = $SUBST{TACGLIB}   = "/usr/local/lib/tacg"; # the TACGLIB directory containing the
                                                        # default tacg data
$SYSTMP    = $SUBST{SYSTMP}    = "/tmp"; # the system tmp directory where the cgi can write
                                         # temporary files.  Usually /tmp, /usr/tmp, or /var/tmp
$CGITMPDIR = $SUBST{CGITMPDIR} = "tmp";  # the name of the tmp dir in the cgi-bin dir that
                                         # gets linked to the system tmp dir.
$WWWUSER   = $SUBST{WWWUSER}   = "www-data";  # the user id you want the tacg CGI to execute as
$WWWGROUP  = $SUBST{WWWGROUP}  = "www-data"; # the group id you want the tacg CGI to execute as

my $CONFIGURED = 0; # set this to 1 after you've finished setting the configs.


# Assign all the vars to their SUBST equivalents

if ($CONFIGURED == 0) {

   print << "INTIAL_WARNING";

      You need to read and edit this file to configure the setting for
      your installation.  After you do this, setting \$CONFIGURED to 1
      (line 74) and re-run the script, the installation should run
      correctly, if you have the right permissions.

INTIAL_WARNING
   exit();
}


print << "WELCOME";

	OK - it looks like you've set the variables as advised.  Variable
   substitution and installation will continue now.  Here we go.


WELCOME

if ($ENV{USER} !~ "root") {
    $USER = "~" . $ENV{USER};
#    $TACG_DIR  = $SUBST{TACG_DIR} = "$USER/$TACG_DIR";
    $HTML_LINK = $SUBST{HTML_LINK} = "$TACG_DIR";
    $CGI_LINK  = $SUBST{CGI_LINK}  = "cgi-bin/$TACG_DIR";
#    $HTML_LINK = $SUBST{HTML_LINK} = "~" . "$ENV{USER}/$TACG_DIR";
#    $CGI_LINK  = $SUBST{CGI_LINK}  = "~" . "$ENV{USER}/cgi-bin/$TACG_DIR";

    print <<"NOTROOT";

    Hmmmm ... I see that you're not installing this as root - be aware
    that depending  on where you're trying to install and for whom,
    this may fail, as it's been designed to install as root. Check the
    error messages carefully to see what has failed.  However, I'll
    try to install it in your ($USER) public_html dir.

NOTROOT

} else {
    $USER = "";
    $HTML_LINK = $SUBST{HTML_LINK} = $TACG_DIR;
    $CGI_LINK  = $SUBST{CGI_LINK} =  "cgi-bin/$TACG_DIR";

}



# The only files that need to be substituted are:
@Files = ( "tacgi4-cgi.pl",
           "form4.html",
           "index.html"
          );

# go forth and substitute...
# my ($file, $input_file);

foreach my $file (@Files) {
   my $input_file = $file . ".in";
   open(IN, $input_file) or die "Couldn't open $input_file for reading";
   open(OUT,">$file") or die "Couldn't open $file for writing";
   print STDERR "\nCreating $file from $input_file \n\n";
  while (<IN>) {
    next unless /\%\%[^\%]+\%\%/;
    foreach my $pattern (keys %SUBST) {
      if ($_ =~ /\%\%$pattern\%\%/) {
         print STDERR ".";
         $_ =~  s/\%\%$pattern\%\%/$SUBST{$pattern}/g;
      }
    }
  } continue {
    print OUT;
  }
  print STDERR "\n";
  close(OUT);
  close(IN);
}

# now install all the files where they need to be to run correctly.

# ==============  create the required directories. =========================
# -------------- 1st the SYSTMP
system("mkdir -p -m 755 $SYSTMP") unless -d $SYSTMP;


# -------------- then the HTML tacg dir
if (! -d $HTMLDIR) {
	print << "WEB";

   Looks like your configuration edits don't match your httpd config
   (I can't find a pre-existing $HTMLDIR)
   you need a working web server which has been configured to use
   $HTMLDIR as a valid DocumentRoot.

   But I'll try to make it anyway...

WEB
   system("mkdir -p -m 755 $HTMLDIR") unless -d $HTMLDIR;

}
if (! -d "$HTMLDIR/$TACG_DIR") {
	print << "NOTACGDIR";

   Looks like your TACG subdirectory doesn't exist yet as well.
   (I can't find a pre-existing $HTMLDIR/$TACG_DIR)

   But I'll try to make it anyway...

NOTACGDIR
   system("mkdir -p -m 755 '$HTMLDIR/$TACG_DIR'") unless -d '$HTMLDIR\/$TACG_DIR';

#	print "\nCreating $HTMLDIR/$TACG_DIR ..\n\n";
#	system("mkdir -p -m 755 $HTMLDIR/$TACG_DIR")	 unless -d "$HTMLDIR/$TACG_DIR";
}

# -------------- then the CGI tacg dir

if (! -d $CGI_DIR) {
	print << "WEB";

   Looks like your configuration edits don't match your httpd config -
   you need a working web server which has been configured to use
   $CGI_DIR as it's ScriptAlias..
   But I'll try to create it on the fly anyway...

WEB
   system("mkdir -p -m 755 $CGI_DIR")	 unless -d "$CGI_DIR";
}

	print "\nChecking/ Creating $CGI_DIR/$TACG_DIR ..\n\n";
	system("mkdir -p -m 755 $CGI_DIR/$TACG_DIR")	 unless -d "$CGI_DIR/$TACG_DIR";

# -------------- make sure that the cgi tmp -> /tmp link is there
if (! (-d "$CGI_DIR/tmp" && -x "$CGI_DIR/tmp")) {
	print "\nCreating the symlink: $CGI_DIR/tmp -> $SYSTMP..\n\n";
	print "\nCreating the symlink: $HTMLDIR/tmp -> $SYSTMP..\n\n";
	system("ln -s  $SYSTMP  $CGI_DIR/tmp; ln -s  $SYSTMP  $HTMLDIR/tmp");
}


# ======================= then copy the files
print "\nCopying files..\n\n";
system("cp tacgi4-cgi.pl $CGI_DIR/$TACG_DIR");
#system("cd $THIS_DIR ; cp ../Docs/* $HTMLDIR/$TACG_DIR; cp nedit.macros $HTMLDIR/$TACG_DIR");
system("cp tacg_Map.pdf microplasmid.gif tacg-4.3-banner.gif form4.html  index.html $HTMLDIR/$TACG_DIR");
system("cp ../Data/rebase.data $HTMLDIR/$TACG_DIR");
system("cp ../Docs/tacg-4.3-manual.html ../Docs/tacg4.xman.html ../Docs/tacg4.man.html ../Docs/tacg4.0.main.html ../Docs/nedit.macros     $HTMLDIR/$TACG_DIR");


# ======================= Check that the TACGLIB dir exists
if ( -d $TACGLIB && -x $TACGLIB && -r $TACGLIB) {
	print "\nYour $TACGLIB looks to be installed OK..\n\n";
   if (-r "$TACGLIB/codon.data" && -r "$TACGLIB/rebase.data") {
   	print "\nAt least the minimum common data looks to be installed OK..\n\n";
   } else {
   	die "\nERROR: Minimum data 'codon.data' & 'rebase.data' aren't installed. Check!!\n\n";
   }
} else {
	die "\nERROR: Minimum tacg data is missing from $TACGLIB.  It looks like you\n",
       " haven't finished the basic tacg compilation and installation!\n";
}

# ======================= Check that the tacg binary exists and is executable
if ( -x $TACGBIN) {
	print "\nLooks like you've installed tacg ok\n\n";
} else {
	print "\nERROR: Looks like you didn't install the tacg binary (no executable \n",
   "$TACGBIN) - did you finish the install?\n\n";
   exit();
}

# then set the ownership & mode
print "\nSetting Ownership and mode from your settings..\n\n";
system("chown -R $WWWUSER:$WWWGROUP $HTMLDIR/$TACG_DIR; chmod -R a+rX $HTMLDIR/$TACG_DIR");
system("cd $CGI_DIR; chown -R $WWWUSER:$WWWGROUP $TACG_DIR; chmod -R a+rx $TACG_DIR;");

# -------- then a final message to check it
print << "FINITO";

   It looks like the installation went OK.  Please point a browser at:

     http://$WWWHOST/$HTML_LINK/index.html

   and see if it looks right.  If the install worked, the first-listed tacg
   server should be the one you just installed and clicking on the links
   should bring up the correct form and it SHOULD work (touch wood).

   If it doesn't, check your configuration settings & if you change anything,
   run the script again. If it still doesn't work, let me know.
   I can be reached at hjm\@tacgi.com.

   Thanks for trying it out.

FINITO
