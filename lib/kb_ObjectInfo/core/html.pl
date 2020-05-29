#!/usr/bin/perl
$ora_debug = shift if $ARGV[0] =~ /^-#/;

#############################################
#
#name           kegg_fun_cat
#
#description    Create KEGG functional categories web pages
#
#               "Unassigned" means the KEGG top hit is not assigned to a
#               category.  If no KEGG hit, the gene does not make
#               the list.
#
#
#syntax         See die_debug below
#
#prerequisite   KEGG blast must be run and summarized
#
#infile         /analysis/btab_vs_kegg.btab
#               /auto/GAT/db/kegg/redundant
#
#outfile        /analysis/redundant_categories
#               /analysis/redundant2
#               web page /kegg_fc.html
#               web pages /fc_$cat.html (1 per category)
#               
#author         Miriam Land
#
#history        Origin circa 1999
# 
#   10/31/02    Add more comments to the code.
#               Rename from redundant_categories
#               Add a check for medusa to minimize problems

$x = `uname -n`;
unless ($x =~ /maple/ )
{
   print "\a\n\aI'm very sorry to inform you that this script will only \n";
   print "complete correctly if you first log into MAPLE.\n";
   print "The web directories are not accessible from other machines.\n\a\n";
   exit;
}

#*****************************************************
#  Input
#*****************************************************
$org = shift || &die_usage();
$chr = shift || &die_usage();
#$analysis_name = shift || &die_usage();

$root_dir = "/auto/microbe/$org/chromosome/$chr"; 

unless (-e $root_dir)
{
	die "Did not find $root_dir";
}

#*****************************************************
#  Make sure there is a place to output the web files
#*****************************************************

$root_www = "/srv/www/devel/htdocs/microbial/$org";
unless (-e $root_www)
{
   mkdir("$root_www", 0775) || die "cannot mkdir $root_www: $!";
   &perm($root_www,"microbial");
}
$root_www = "/srv/www/devel/htdocs/microbial/$org/$chr";
unless (-e $root_www)
{
   mkdir("$root_www", 0775) || die "cannot mkdir $root_www: $!";
   &perm($root_www,"microbial");
}


#-------------------------------------------------------------------
# Connect to Oracle
#-------------------------------------------------------------------
require ("/auto/microbe/bin/Oracle_connection.pm");
$err = &def_oracle('prod');
&error($err) if ($err gt ' ');


#*****************************************************
#  Look up the full name
#*****************************************************

$query_str = "select distinct species from organism where organism = '$org'";
local $sth = $dbh->prepare("$query_str");
#print "$query_str;\n";
$sth->execute;

$species = $sth->fetchrow_array ;


$btab = "$root_dir/analysis/btab_vs_kegg.btab";
$kegg_red = "/auto/GAT/db/kegg/redundant";
$tmp_red_cat = "/usr/tmp/redundant_categories_$$";

&create_redundant ($btab,$kegg_red,$tmp_red_cat); # Create the list of redundant categories

$redcat = "$root_dir/analysis/redundant_categories";

&sort_file($tmp_red_cat, $redcat);     # Sort the list

$red_cat2 = "$root_dir/analysis/redundant2";
$kegg_fc = "$root_www/kegg_fc.html";

unlink $red_cat2 if (-e $red_cat2);
unlink $btab if (-e $btab);

#
#	From the sorted file, make the master web page and
#	the file which has no known purpose.  We now have a cgi for this.
#
# OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE
#
#&summarize ($redcat,$red_cat2,$kegg_fc);           
# Summarize the results


#-----------------------------------------------------------------
# Create the list of redundant categories
# Use the btab file and the "redundant" list
#  from KEGG.  It has one line for each category
#-----------------------------------------------------------------

sub create_redundant
{
	($btab,$kegg_red,$tmp_red_cat) = @_;

  open (TMP,">$tmp_red_cat");

  $number_to_check = 5;
  $num_genes = 0;

  undef %hash;
  $count = 0;
  $last = 'None';

#-----------------------------------------------------------------
#  Read the btab file and find the KEGG genes which are blast hits
#
#  Make %hash which has all the KEGG genes which hit this organism
#  Make %ornl_gene which has a list of all the KEGG genes which
#    correspond to an ornl gene
#  Option to stop looking and the $number_to_check
#-----------------------------------------------------------------
  open (BTAB,"<$btab") || die "Did not find the btab file $btab";
  while ($buf = <BTAB>)
  {
	next if ($buf =~ /^#/);
#      next unless ($buf =~ /Contig760/);

	@fields = split(/	/,$buf);
	if ($fields[0] ne $last)
	{
		$count = 0;
	}
	$last = $fields[0];
	$count++;
	next if ($count > $number_to_check);

	$g_name = $fields[5];

	if ($fields[0] =~ /\./)
	{
		@name_parts = split(/\./,$fields[0]);
		$or_number = $name_parts[-1];
		$contig = $name_parts[-2];
	}
	else
	{
		@name_parts = split(/ /,$fields[0]);
		$or_number = $name_parts[0];
		$contig = $name_parts[3];
	}

	if ($fields[6] eq "-1")
	{
#		print TMP "Unassigned\tNot similar to any KEGG gene\t$contig $or_number\t\t\n";
		next ;
	}
	$hash{$g_name} = 'N'; # Has it been categoried, default
	
#	print "DEBUG: ORNUM=$or_number CONTIG=$contig HIT=$fields[6]\n";
#      print "DEBUG: added $g_name\t$contig\t$or_number\n";

	$ornl_gene{"$contig $or_number"} .= "$g_name ";

	$num_genes++;
#      last if ($num_genes > 10);
  }

  close BTAB;

#-----------------------------------------------------------------
#  Find the CLASSIFICATION
#
#  Read the redundant KEGG file - keep genes which are found
#    in the organism at hand (the entire file is too large)
#-----------------------------------------------------------------
  open (KEGG,"<$kegg_red");
  while ($kegg_line = <KEGG>)
  {
     ($index,$gene_name,$def,$class) = split(/	/,$kegg_line);
#     next if ($class =~ 'Human Disease');
#     next if ($kegg_line =~ /^rpa/);
     chomp $class;
#
#	Test to see if this gene is needed
#	The oder of the variables is for sorting order
#
     if (exists $hash{$index})
     {
       $class{$index} .= "$class\t$def\t$index\t$gene_name\t";  
#       print "DEBUG:  $index - $class\n";
       if ($class =~ /Unassigned/ && $hash{$index} eq 'N')
       {
         $hash{$index} = 'U';
       }
       else
       {
         $hash{$index} = 'Y';
       }
     }
  }
  close KEGG;

#-----------------------------------------------------------------
#  Match up genes and classification
#  Keep looking through the KEGG IDs until you find one that is
#   classified.  Otherwise, take the first one.
#
#  Print out one line for each category assignment
#-----------------------------------------------------------------
foreach $gene (keys (%ornl_gene))
{
	($contig,$or_number) = split(/ /,$gene);


#    print "GENE=$or_number KEYS=$ornl_gene{$gene}\n";

	@ary = split(/ /,$ornl_gene{$gene}); # list of genes
#
#	Find the first gene with a classification
#	Look through the blast hits until you find one
#
	$first = 'None';
	foreach $index (@ary)
	{
		$classified = $hash{$index}; # Did it get classified Y/N/U
		next if ($classified eq 'U' || $classified eq 'N');
		$first = $index;
		last;
	}

	if ($first eq 'None')  
	{
		print TMP "Unassigned\tNot similar to any KEGG gene\t$contig $or_number\t\t\n";
#		print  "Unassigned\t$contig $or_number\t$ornl_gene{$gene}\n";
		next;
	}

	$line = $class{$first};
	chop $line;  # Eliminate last tab
	@lines = split(/	/,$line);

	while ($#lines > 1)
	{

		($class,$def,$index,$gene_name) = splice(@lines,0,4);
		print TMP "$class\t$def\t$contig $or_number\t$index\t$gene_name\n";  
#		print "DEBUG: $class\t$def\t$contig $or_number\t$index\t$gene_name\n"  ;  
	}
}
close TMP;

}


#********************************************
#      Sort the results
#********************************************

sub sort_file
{
  local ($tmp_red_cat, $redcat) = @_;
  $cmd = "sort $tmp_red_cat > $redcat";
  system ($cmd);
  &perm($redcat,"microbial");
  unlink $tmp_red_cat;
}

#********************************************
#      Summarize the outcome
#
#	Read the sorted, redundant categories and
#	make web pages for each category
#********************************************

#
# OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE
#
sub summarize
{

	($redcat,$red_cat2,$kegg_fc) = @_;
	open (IN,"<$redcat") || die ("Could not open $redcat\n");
	@color = ("#FFFFFF","#FFFFCC","#CCFFFF");

	open (SORT,">$red_cat2") ;
	open (FC,">$kegg_fc");
	&header3;

	@heads = ('Contig', 'Gene','KEGG Hit','KEGG Gene Name','KEGG Desc');

	$web_page = 'None';
	$page_cnt = 0;
	$last_cat = 'None';

	$table = 'None';
	$table_cnt = 0;

while ($buf = <IN>)
{
	($test,$def,@fields) = split(/	/,$buf);
	next if ($test lt 'A');
	push(@fields,$def);
#    if ($test =~ /^Metabolism;/ || $test =~ /^Unassigned/)  ## OLD STYLE
#    {
	($cat1,$cat2,$cat3) = split(/;/,$test);
	next if ($cat3 lt '     ');

#    }
#    else                                                    ## NEW STYLE
#    {
#      ($cat2,$dummy,$cat3) = split(/;/,$test);
#    }
	$name = "$cat2";
	$name = "$cat1 -- $cat2" if ($cat1 =~ /^Un/);
	$cat2 =~ s/^ //;
	$cat2 =~ s/ /_/g;
	$cat2 = 'Unassigned' if ($cat2 lt '     ');

    if ($cat1 ne $last_cat)
    {
       print FC "</blockquote>" if ($web_cnt > 0);
       print FC "<h3>$cat1</h3><blockquote>";
    }

    if ($web_page ne $cat2)
    {
       print SORT "$test\n";

       if ($web_cnt > 0)
       {
           &endTable; 
           &footer; 
           close CAT;
           &perm($cat_html,"microbial");
       }

       $cat_html = "$root_www/fc_$cat2.html";
       open (CAT,">$cat_html");
       &header($cat1,$name);   

       print FC "<a href=\"fc_$cat2.html\">$name</a><br>";
 
       $web_cnt++;
       $table_cnt = 0;

       $web_page = $cat2;
       $table = 'None';
    }
    if ($table ne $cat3)
    {
       if ($table_cnt > 0)
       {  &endTable; }
       print SORT "$test\n";
       print CAT "<h2>$cat3</h2>\n";
       &startTable("2",@heads);
       $table = $cat3;

    }
    print SORT "\t@fields";

    ($contig,$gid) = split(/ /,$fields[0]);
    $gene_number = $gid;
    $CHR=$chr;
    if ($chr eq 'combined')
    {
       $CHR = $contig;
       $CHR =~ s/Contig//;
    }
    $gid = "<a href=\"/cgi-bin/JGI_microbial/gene_viewer.cgi?org=$org&chr=$CHR&contig=$contig&gene=$gid\">$gid</a>";

#
#  Don't have a way to do this yet because multiple chromosomes
#    have the same sets of numbers and the hash $contig_link{}
#    has them mixed together. If dealing with 'combined', don't link
#

   $link = "http://maple.lsd.ornl.gov/cgi-bin/JGI_microbial/contig_viewer.cgi?org=$org&chr=$CHR&contig=$contig";
   $contig = "<a href=\"$link\">$contig</a>";
#   if (exists($contig_link{$contig}) && $chr ne 'combined')
#   {
#      $link = $contig_link{$contig};
#      $contig = "<a href=\"$link\">$contig</a>";
#   }
#   elsif (exists($contig_link{$gene_number}) && $chr ne 'combined')
#   {
#      $link = $contig_link{$gene_number};
#      $contig = "<a href=\"$link\">$contig</a>";
#   }

    $fields[0] = "$contig</td><td>$gid";

    $fields[2] = "&nbsp" if ($fields[2] le '   ');
    $color = $color[$table_cnt % 2];
    printTableRow($color,@fields);
       $table_cnt++;

    $last_cat = $cat1;

}  ## End of while loop reading redundant categories

print FC "</blockquote>";

&endTable;
&footer;
&footer3;


close SORT;
close CAT;
close IN;
close FC;

&perm($red_cat2,"microbial");
&perm($cat_html,"microbial");
&perm($kegg_fc,"microbial");

}

#********************************************
#      FOOTER
#********************************************
sub footer
{
  $foot = "/srv/www/prod/htdocs/footer.htmlf";
  open (FOOT,"$foot");
  @foot_text = <FOOT>;
  close FOOT;

  print CAT "@foot_text";

}
#********************************************
#      FOOTER
#********************************************
sub footer3
{
  $foot = "/srv/www/prod/htdocs/footer.htmlf";
  open (FOOT,"$foot");
  @foot_text = <FOOT>;
  close FOOT;

  print FC "@foot_text";

}


#********************************************
#      HEADER
#********************************************
sub header
{
  local ($cat1,$name) = @_;
  $title = "Functional Categories";
  $head = "/srv/www/prod/htdocs/GCat/header2.htmlf";
  open (HEAD,"$head");
  @head_text = <HEAD>;
  close HEAD;


  print CAT qq!
<HTML>
<HEAD>
<meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
<TITLE>$title</TITLE>
</HEAD>
<BODY BGCOLOR="#FFFFFF" BACKGROUND="">
@head_text
<H1><CENTER>Category Assignments of top KEGG Hits</CENTER></H1>
<h2><center>$cat1</center></h2>
<h2><center>$name</center></h2>
<h3><center>$species</center></h3>
!;

}

#********************************************
#      HEADER
#********************************************
sub header3
{
  $title = "Functional Categories";
  $head = "/srv/www/prod/htdocs/GCat/header2.htmlf";
  open (HEAD,"$head");
  @head_text = <HEAD>;
  close HEAD;


  print FC qq!
<HTML>
<HEAD>
<meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
<TITLE>$title</TITLE>
</HEAD>
<BODY BGCOLOR="#FFFFFF" BACKGROUND="">
@head_text
<H1><CENTER>Category Assignments of top KEGG Hits</CENTER></H1>
<h2><center>$species</center></h2>
!;

}


#********************************************
#      End a Table
#********************************************
sub endTable
{
  print CAT "</table>\n";
}

#*****************************************************
#  Print a row for a table
#*****************************************************

sub printTableRow 
{
  my ($color,@cols) = @_;
  $\ = "\n";

  print CAT "<tr bgcolor=$color>";
  local $count = 0;
  foreach (@cols)
  {
    print CAT "<td align=$align[$count]>$_</td>";
    $count++;
  }
  print CAT "</tr>";
}

#*****************************************************
#  Start a table
#*****************************************************

sub startTable
{
  my ($cell,@cols) = @_;
  $\ = "\n";

  print CAT "<table cellspacing='$cell' cellpadding='$cell' border=1>";
  print CAT "<tr>";
  foreach (@cols)
  {
    print CAT "<th>$_</th>";
  }
  print CAT "</tr>";

}
#--------------------------------------------------
#  Set permissions 
#  The file must exist and the user must be the owner
#--------------------------------------------------
sub perm
{
  my ($locfile,$group) = @_;
  $locfile =~ s/\>//g;   # Many files have >> at the front
  return unless(-e $locfile && -o $locfile);
  system ("chmod 775 $locfile");
  system ("chgrp $group $locfile");
}



sub die_usage()
{
  die "Usage:        $0 <org>  <chrom>\nPrerequisite: Summarized btab_vs_kegg.btab file";
}
