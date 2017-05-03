


== INSTALATION ==

On UNIX-like operating system (i.e., Mac OS X  and Linux), type the following commands in a terminal
window:

aclocal;
autoheader;
autoconf -f;
automake -f --add-missing;
./configure;
make clean;
make;

In order  to install the program  'coltree' that will help  you visualise the results  of a fitmodel
analysis, type the following commands:

aclocal;
autoheader;
autoconf -f;
automake -f --add-missing;
./configure  --enable-coltree;
make clean;
make;



== DESCRIPTION ==

'Fitmodel' estimates  the parameters of various  codon-based models of
substitution, including those described  in Guindon, Rodrigo, Dyer and
Huelsenbeck  (2004).  These  models  are  especially  useful  as  they
accommodate site-specific switches between selection regimes without a
priori  knowledge  of the  positions  in  the  tree where  changes  of
selection regimes occurred.

The program will ask  for two input files: a tree  file and a sequence
file.  The  tree should  be unrooted and  in NEWICK  format.  fitmodel
will crash if  you feed him with a rooted  tree.  The sequences should
be in PHYLIP interleaved or sequential format.  If you are planning to
use codon-based models, the sequence length should be a multiple of 3.
The program provides  four types of codon models: M1,  M2, M2a, and M3
(see PAML  manual).  Moreover,  M2, M2a  and M3  can be  combined with
'switching'   models  (option   'M').    Two   switching  models   are
implemented: S1 and  S2.  S1 constraints the rates  of changes between
dN/dS  values  to be  uniform  (e.g.,  the  rates of  changes  between
negative and positive  selection is constrained to be the  same as the
rate of  change between  neutrality and  positive selection)  while S2
allows for differents rates of change between the different classes of
dN/dS values.

If you are using a 'switching' model, 'fitmodel' will output file with
the        following        names:        your_sequence_file_trees_w1,
your_sequence_file_trees_w2,      your_sequence_file_trees_w3      and
your_sequence_file_trees_wbest.   The w1,  w2  and w3  files give  the
estimated tree  with probabilities  of w1, w2,  and w3  (three maximum
likelihood dN/dS ratio estimates) calculated  on each edge of the tree
and  for each  site.  Hence,  the  first tree  in one  of these  files
reports  the  probabilities  calculated  at  the  first  site  of  the
alignment.  Instead of  probabilities, the  wbest file  allows you  to
identify which  of the  tree dN/dS  is the most  probable on  any give
edge, at any given site. A branch  with label 0.0 means that w1 is the
most probable class, 0.5 indicates the w2 is the most probable and 1.0
means that w3 has the highest posterior probability.

Processing those  results manually can  be very tedious.   The program
'coltree' reads 'fitmodel' output file and generates a postscript file
which displays the phylogeny at  each codon site with different colors
on edges  depending on the  posterior probabilities estimated  on each
branch, at each site (see file p1.pdf).  This is very useful to detect
site specific changes  of selection regimes.  To  compile 'coltree' on
UNIX-like systems,  simply type the  following commands in  a terminal
window:


This  will  generate   the  binary  for  'coltree'.   'coltree'  is  a
command-line program  that takes  as first argument  'fitmodel' output
tree file.  The second  argument is  the name  of the  postscript file
where the trees will be written.


== EXAMPLE ==

The directory  example/ contains a  sequence file (p1.nxs) and  a tree
file  (p1.tree)  that  can  be   processed  by  fitmodel.   The  files
p1.nxs_fitmodel_stats,      p1.nxs_fitmodel_tree,     p1.nxs_trees_w1,
p1.nxs_trees_w2, p1.nxs_trees_w3 and p1.nxs_trees_wbest were generated
by 'fitmodel' using the M2a+S1 model.   The command line that was used
in order to generate those files is the following one:

./fitmodel < param

where 'param' is  a text file that contains the  series of letters one
have to enter in order to  select the M2a+S1 model (and other options)
through the PHYLIP-like interface.

The file  'tree.ps' is a  postscript file that displays  the posterior
probability of w1 on each edge at  each code site of the alignment. It
was generated using the following command-line:

./coltree p1.nxs_trees_w1 tree.ps



====

Need more help? Contact me: s.guindon@auckland.ac.nz



BIBLIOGRAPHY

"Modeling the site-specific variation of selection patterns along lineages"
Stephane Guindon, Allen G Rodrigo, Kelly A Dyer, John P Huelsenbeck
Proceedings of the National Academy of Sciences of the United States of America
2004

"Modelling the variability of evolutionary processes"
Olivier Gascuel, Stephane Guindon
Reconstructing evolution. New mathematical and computational advances.
Oxford University Press
2007
