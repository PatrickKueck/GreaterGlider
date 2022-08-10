#!/usr/bin/perl

use strict				;
use warnings			;
use File::Copy			;


############################################# General Script Description ###################################################
############################################################################################################################
# greaterGlider.pl is a perl written software pipeline for linux operating systems, allowing sequence simulation (INDelible)
# combined with ML tree reconstruction (IQTree) and subsequent conducted branch support estimation using the site-concordance 
# (sCF) measure (IQTree2). To enable all three stages of analysis (simulation, ML recosntruction, and branch support measure),
# allowing (among other things) a re-analysis of the original sCF study of Kück et al 2022, the pipeline consists a reduced 
# source code of the originally used Appetite-pipeline. Processing greaterGlider.pl in default mode (perl greaterGlider <enter>)
# conducts a re-analysis of the ten tree-setups /T1 to T10) as published in Kück et al 2022 with alpha=0.5. For another alpha 
# value used by Kück et al. 2022, change the second value of the %indelible_parameter 'rates' and re-start the script.
#
# Additionally to the simulations of the original study of Kück et al. 2022, greaterGlider.pl allows also other tree analyses 
# in the frame of pipeline specified parameter space (see pipeline parameter descriptions below). To change simulation underlying
# tree shapes and branch length conditions, simply substitute the newick strings in %tree_parameter. To perform a ML
# reconstruction without sCF measure, simply comment out the processing of the sCF measure (&concordance) in the 'Pipe Start'
# block, and so on...
#
# Further descriptions of greaterGlider.pl handling as well as input and output data are given at the beginning of each of 
# the three source code blocks below (see block 'INDELIBLE Sequence Simulation', 'IQTREE Tree Reconstruction', and
# 'IQTREE2 Concordance Measure').
#
# !IMPORTANT!: To use greaterGlider.pl following software must be pre-installed:
#	:: INDELible v1.03 
#		:: Fletcher W. and Yang Z. 2009. INDELible: a flexible simulator of biological sequence evolution. Mol. Biol. and Evol. 2009 26(8):1879-1888
#	:: IQTree v1.6.12 
#		:: Nguyen et al. 2015. IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Mol. Biol. and Evol. 2015 32(1):268-274
#	:: IQTree v2.2.0
#		:: Minh et al. 2020. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. and Evol. 2015 37(5):1530-1534
#	:: pdfunite
#		:: part of the TeX Live package (https://tug.org/texlive/) 
############################################################################################################################
############################################################################################################################


############################################# pipeline parameter ###########################################################
############################################################################################################################

#############################################
#	- %tree_parameter:
#		:: input tree(s) for INDELible sequence simulation, with or without '\$BL3\$' coded branch elongations
#		:: hashkey (e.g. 'T1') used as tree code
#		:: hashvalue: tree code corresponding newick string ('\$BL3\$' coded branches are stepwise elongated with a single fasta file for each branch elongation)
#		:: for setting of '\$BL3\$' coded branch elongations modify %indelible_parameter: 'minBL3', 'maxBL3', and 'stpBL3'
#		:: the ten different tree settings below represent the ten tree simulations used in the sCF study of Kück et al 2022
my %tree_parameter		=	(
	
	'T1'						=>	"((((A1:0.001,A2:0.001):\$BL3\$,(B1:0.001,B2:0.001):\$BL3\$):0.01,((C1:0.001,C2:0.001):0.1,(D1:0.001,D2:0.001):0.1):0.01):0.05,(O1:0.001,O2:0.001):0.05);"	,
	'T2'						=>	"(((A1:0.001,A2:0.001):\$BL3\$,((B1:0.001,B2:0.001):0.1,((C1:0.001,C2:0.001):0.1,(D1:0.001,D2:0.001):0.1):0.01):0.01):0.05,(O1:0.001,O2:0.001):\$BL3\$);"	,
	'T3'						=>	"(((D1:0.001,D2:0.001):0.1,((C1:0.001,C2:0.001):0.1,((A1:0.001,A2:0.001):\$BL3\$,(B1:0.001,B2:0.001):\$BL3\$):0.01):0.01):0.05,(O1:0.001,O2:0.001):0.05);"	,
	'T4'						=>	"((((A1:0.001,A2:0.001):\$BL3\$,(B1:0.001,B2:0.001):0.1):0.01,((C1:0.001,C2:0.001):\$BL3\$,(D1:0.001,D2:0.001):0.1):0.01):0.05,(O1:0.001,O2:0.001):0.05);"	,
	'T5'						=>	"(((A1:0.001,A2:0.001):0.1,((B1:0.001,B2:0.001):0.1,((C1:0.001,C2:0.001):\$BL3\$,(D1:0.001,D2:0.001):0.1):0.01):0.01):0.05,(O1:0.001,O2:0.001):\$BL3\$);"	,
	'T6'						=>	"(((A1:0.001,A2:0.001):\$BL3\$,((B1:0.001,B2:0.001):0.1,((C1:0.001,C2:0.001):\$BL3\$,(D1:0.001,D2:0.001):0.1):0.01):0.01):0.05,(O1:0.001,O2:0.001):0.05);"	,
	'T7'						=>	"((((A1:0.001,A2:0.001):0.1,(B1:0.001,B2:0.001):\$BL3\$):0.01,((C1:0.001,C2:0.001):0.1,(D1:0.001,D2:0.001):0.1):0.01):0.05,(O1:0.001,O2:0.001):\$BL3\$);"	,
	'T8'						=>	"(((A1:0.001,A2:0.001):0.1,((B1:0.001,B2:0.001):\$BL3\$,((C1:0.001,C2:0.001):0.1,(D1:0.001,D2:0.001):0.1):0.01):0.01):0.05,(O1:0.001,O2:0.001):\$BL3\$);"	,
	'T9'						=>	"(((A1:0.001,A2:0.001):\$BL3\$,((B1:0.001,B2:0.001):\$BL3\$,((C1:0.001,C2:0.001):0.1,(D1:0.001,D2:0.001):0.1):0.01):0.01):0.05,(O1:0.001,O2:0.001):0.05);"	,
	'T10'						=>	"(((D1:0.001,D2:0.001):0.1,((C1:0.001,C2:0.001):\$BL3\$,((B1:0.001,B2:0.001):\$BL3\$,(A1:0.001,A2:0.001):0.1):0.01):0.01):0.05,(O1:0.001,O2:0.001):0.05);"	,
) ;
#############################################

#############################################
#	- %indelible_parameter
#		:: parameter setting of INDELible sequence simulation following each %tree_parameter specified newick tree setup 
#		:: to perform sequence simulation, INDELible must be preinstalled  
#			:: change parameter 'scr' in respect of your INDELible terminal executable command
#		:: to change the alpha value of sequence simulation, adapt the second value of 'rates'
my %indelible_parameter		=	( 
	
	'scr'		  				=>	'indelible'								,	# INDElible parameter: script name
	'TYPE'					 	=>	'NUCLEOTIDE 1'							,	# INDElible parameter: 'NUCLEOTIDE 1', 'NUCLEOTIDE 2', 'AMINOACID 1', 'AMINOACID 2', 'CODON 1', 'CODON 2' -> simulation under nucleotide, amino-acid or codon models and which algorithm to use
	'ancestralprint'			=>	'NEW'									,	# INDElible parameter: ancestral sequences printed in new file ('NEW'), same file as recent sequences ('SAME') or not printed ('FALSE')
	'output'					=>	'FASTA'									,	# INDElible parameter: TRUE alignments either in 'FASTA', 'NEXUS', or 'PHYLIP' format, Unaligned sequences are always output in FASTA format as e.g. filename.fas
	'randomseed'				=>	12345									,	# INDElible parameter: start random-seed number
	'printrates'				=>	'TRUE' 									,	# INDElible parameter: 'TRUE' or 'FALSE'. 'TRUE' will print out a detailed output file about the site-classes or relative rates of substitution. FALSE will not print any rates information.
	'insertaslowercase'			=>	'TRUE'									,	# INDElible parameter: 'TRUE' or 'FALSE'. 'TRUE' will output inserted bases/residues in lower case letters for easy identification
	'markdeletedinsertions'		=>	'FALSE'									,	# INDElible parameter: 'TRUE' or 'FALSE'. 'TRUE' will output inserted bases/residues that have been subsequently been deleted as * rather than - for easy identification.
	'printcodonsasaminoacids'	=>	'FALSE'									,	# INDElible parameter: 'TRUE' or 'FALSE'. 'TRUE' will output codon datasets as amino-acids
	'fileperrep'				=>	'TRUE'									,	# INDElible parameter: 'TRUE' or 'FALSE'. 'TRUE will output each replicate dataset in a separate file.
	'MODEL'						=>	'mymodelname'							,	# INDElible parameter: specified model name, can be anything of alphanumeric signs
	'submodel'					=>	'GTR 0.2 0.4 0.6 0.8 1.2'				,	# INDElible parameter: For evolutionary model specifications visit indelible online manual
	'indelmodel'				=>	'LAV 1.7 20'							,	# INDElible parameter: processed only if defined. specifies the indel length distribution
	'indelrate'					=>	''	 									,	# INDElible parameter: processed only if defined. float below 1 (e.g. 0.1), rates of insertion and deletion are both 0.1
	'geneticcode'				=>	''										,	# INDElible parameter: processed only if defined. only used in CODON simulations. The value should be an integer 1 to 6, 9 to 16, or 21 to 24, corresponding to the genetic codes listed on Genbank.
	'rates'						=>	'0.3 0.5 0'								,	# INDElible parameter: pinv (float between 0 and 1) alpha (positive float number) ngamcat (integer, 0 -> continous gamma distribution)
	'statefreq' 				=>	'0.25 0.25 0.25 0.25'					,	# INDElible parameter: frequencies for T C A G, must sum up to 1
	'minBL3'					=>	0.1										,	# INDElible parameter: start length bl3, must be greater 0 if defined in newick
	'maxBL3'					=>	1.5 									,	# INDElible parameter: end length bl3, must be greater as start length if defined in newick
	'stpBL3'					=>	0.2										,	# INDElible parameter: steps of length increase bl3, must be greater 0 if defined in newick
	'rootlength' 				=>	250000									,	# INDElible parameter: sequence length (integer)
	'replicates'				=>	1										,	# INDElible parameter: number replicates of each tree (default: 1) given same random seed
	'Nloops'					=>	100										,	# INDElible parameter: original 100 ## number replicates of each tree (default: 1) given different random seeds
	'log'						=>	'indelible_logfile'						,	# INDElible parameter: storing of indelible terminal output
	'plog'						=>	1										,	# INDElible parameter: 1 -> Print log file, 0 -> no logfile
);
#############################################

#############################################
#	- %iqtree_parameter
#		:: parameter setting of IQTree ML reconstruction
#		:: to perform an IQTree ML analysis, IQTree must be preinstalled
#			:: to ensure that all parameter options are treated reliable, we recommend to pre-install IQTree v1.6.12
#			:: change parameter 'scr' in respect of your IQTree terminal executable command
my %iqtree_parameter		=	(

	'scr'						=>	'iqtree'								,	# IQTree parameter: scriptname (iqtree or iqtree2)
	'smod'						=>	'GTR'									,	# IQTree parameter: subst. model (GTR def for nuc) other models, see manual description -m <model_name> DNA: HKY (default), JC, F81, K2P, K3P, K81uf, TN/TrN, TNef, TIM, TIMef, TVM, TVMef, SYM, GTR, or Protein: LG (default), Poisson, cpREV, mtREV, Dayhoff, mtMAM, JTT, WAG, mtART, mtZOA, VT, rtREV, DCMut, PMB, HIVb, HIVw, JTTDCMut, FLU, Blosum62, GTR20, MODEL-FINDER: '-m TESTONLY' -> Standard model selection (like jModelTest, ProtTest) '-m TEST' -> Standard model selection followed by tree inference '-m MF' -> Extended model selection with FreeRate heterogeneity '-m MFP' Extended model selection followed by tree inference  '-m TESTMERGEONLY' -> Find best partition scheme (like PartitionFinder) '-m TESTMERGE' Find best partition scheme followed by tree inference '-m MF+MERGE' -> Find best partition scheme incl. FreeRate heterogeneity  '-m MFP+MERGE' -> Like -m MF+MERGE followed by tree inference
	'inva'						=>	1										,	# IQTree parameter: 1 -> model + prop. of invaraible sites (I), 0 -> without I proportion
	'pinv'						=>	0										,	# IQTree parameter: 1 -> '-i <float>' Proportion of invariable sites, 0 -> default: estimate
	'gama'						=>	1										,	# IQTree parameter: 1 -> '+G' otherwise no rate heterogeneity estimation (gamma distribution)
	'alpha'						=>	0										,	# IQTree parameter: 0 -> estimate of alpha, positive float -> -> model + gamma distribution ('-a 1')
	'alrt'						=>	0										,	# IQTree parameter: 1 -> SH-aLRT test, 0 -> without SH-aLRT test
	'alrts'						=>	0										,	# IQTree parameter: '-alrt <#replicates>' -> SH-like approximate likelihood ratio test (SH-aLRT), 'alrt 0' -> Parametric aLRT test (Anisimova and Gascuel 2006)
	'mboot'						=>	'b'										,	# IQTree parameter: STANDARD NON-PARAMETRIC BOOTSTRAP:'-b <#replicates>' -> boot+ML+Cons, '-bc <#replicates>' -> boot+Cons, '-bo <#replicates>' -> boot only; ULTRAFAST BOOTSTRAP: '-bb <#replicates>'
	'wbtl'						=>	1										,	# IQTree parameter: 1 -> '-wbtl' Write bootstrap trees to .ufboot file (default: 0 -> none)
	'nmit'						=>	1000									,	# IQTree parameter: '-nm <#iterations>' Maximum number of iterations (default: 1000)
	'nstb'						=>	100										,	# IQTree parameter: '-nstep <#iterations>' #Iterations for UFBoot stopping rule (default: 100)
	'nboot'						=>	0										,	# IQTree parameter: N bootstrap replicates (integer > 0 ), 0 -> no bootstrapping
	'core'						=>	1										,	# IQTree parameter: number of cpu's to speed up computation time
	'wltr'						=>	1										,	# IQTree parameter: 1 -> Write locally optimal trees into .treels ('-wt') file, 0 do not write local trees
	'wsrc'						=>	0										,	# IQTree parameter: 1 -> Write site rates and categories to .rate file '-wsr'
	'wrsl'						=>	1										,	# IQTree parameter: 1 -> Write site log-likelihoods to .sitelh file, 0 -> don't '-wsl'
	'wslr'						=>	1										,	# IQTree parameter: 1 -> Write site log-likelihoods per rate category '-wslr'
	'part'						=>	''										,	# IQTree parameter: q -> Edge-linked partition model (file in NEXUS/RAxML format), must be specified with inputfolder -q <partition_file> or '-spp <partition_file>' Like -q option but allowing partition-specific rates
	'ftre'						=>	0										,	# IQTree parameter: 1 -> '-te <user_tree_file> Like -t but fixing user tree (no tree search performed), 0 -> don't
	'pref'						=>	''										,	# IQTree parameter: -pre <PREFIX> Prefix for all output files ('' -> default: aln/partition)
	'quie'						=>	1										,	# IQTree parameter: '-quiet' 1 -> Quiet mode, suppress printing to screen (stdout), 0 -> don't
	'mram'						=>	0										,	# IQTree parameter: '-mem RAM' 1->  Maximal RAM usage for memory saving mode, 0 -> don't, 0 -> don't
	'optg'						=>	1										,	# IQTree parameter: 1 -> '--opt-gamma-inv' More thorough estimation for +I+G model parameters
	'mh'						=>	1										,	# IQTree parameter: 1 -> '-mh' Computing site-specific rates to .mhrate file using Meyer & von Haeseler (2003) method
	'log'						=>	'iqtree_logfile'						,	# IQTree parameter: storing of iqtree terminal output
	'plog'						=>	1										,	# IQTree parameter: 1 -> Print log file, 0 -> no logfile
	'root'						=>	'O1'									,	# IQTree parameter: if defined root tree with defined taxon (be aware of correct lettering)
);
#############################################

#############################################
#	- %cf_parameter
#		:: parameter setting of sCF analysis given one or multiple ML trees for each alignment
#		:: to perform a sCF analysis, IQTree2 must be preinstalled
#			:: to ensure that all parameter options are treated reliable, we recommend to pre-install IQTree v2.2.0
#			:: change parameter 'scr' in respect of your IQTree terminal executable command
my %cf_parameter			=	(
	
	'scr'						=>	'iqtree2'								,	# CF parameter: scriptname
	'cf_type'					=>	's'										,	# CF parameter: 's' -> site-concordance measure; 'g' -> gene-concordance measure (gene-concordance measure not implemented yet)
	'rtree'						=>	'a'										,	# CF parameter: Name of reference tree ('name.tre') to assign concordance factor (sCF & gCF), for more than one tree type 'a'
	'msa'						=>	'a'										,	# CF parameter: Name of the data alignment for site-concordance measure ('name.fas'), for more than one alignment type 'a'
	'procscf'					=>	'a'										,	# CF parameter: 'm' -> one site concordance measure given matching reference tree and alignment names, 'a' -> sCF measure for each ref tree in respect of all alignments
	'dftree'					=>	0										,	# CF parameter: 1 -> Write discordant trees associated with gDF1, 0 -> no printing
	'scfNUM'					=>	1000									,	# CF parameter: Number of quartets for site concordance factor (sCF)
	'pFILE'						=>	''										,	# CF parameter: Partition file or directory for site-concardance factor calculation
	'cfverbose'					=>	1										,	# CF parameter: 1 -> Write CF per tree/locus to cf.stat_tree/_loci, 0 -> no printing
	'cfquartet'					=>	1										,	# CF parameter: 1 -> Write sCF for all resampled quartets to .cf.quartet, 0 -> no printing
	'log'						=>	'iqtree_logfile'						,	# CF parameter: storing of iqtree terminal output
	'sim'						=>	1										,	# CF parameter: 1 -> extra output of substring cf support in respect of different branch lengths BL2; 0 -> no extra print (for all non-BL2 coded data) 
	'subnwk'					=>	1										,	# CF parameter: 1 -> substring newick results without further inner bracket structure e.g. (A,B,C,D):65, 0 -> with inner bracket structure e.g. ((A,B),(C,D)):65  
	'procsys'					=>	1										,	# CF parameter: 1 -> iqtree2 process execution, 0 -> no process execution, allows re-analysis based on pre-existing CF-measure results 
);
############################################################################################################################
############################################################################################################################



############################################# PIPE START ###################################################################
############################################################################################################################
my	$indelible_resultfolder	=	"001_indelible" ;

for my $ptree	( sort keys %tree_parameter ){
	
	my $msa_indata			=	$ptree."/".$indelible_resultfolder	;
	my $iqtree_resultfolder	=	$ptree."/002_iqtree"				;
	my $cf_resultfolder		=	$ptree."/003_sCFmeasure"			;
	
	&indelible		( \%indelible_parameter, \$ptree, \$tree_parameter{$ptree}, \$indelible_resultfolder );
	&iqtree			( \%iqtree_parameter, \$ptree, \$iqtree_resultfolder, \$msa_indata ) ;
	&concordance	( \%cf_parameter, \$ptree, \$cf_resultfolder, \$msa_indata, \$iqtree_resultfolder ) ;
	
	exit
}
############################################################################################################################
############################################################################################################################



############################################# INDELIBLE Sequence Simulation ################################################
############################################################################################################################

############################################# Brief Description &indelible #################################################
############################################################################################################################
#
# aINDELible.pm: a module to simulate sequence data due to different branch length combinations of a given topology
#
#	- reading infile data
#---------------------------------------------------------------------------------------------------------------------------
#		: newick branch length coded tree
#						:: branch lengths must be either defined directly in the tree or by branch length codes ($BL3$)  
#						:: e.g.: (((((A1:0.1,A2:0.1):$BL3$,(B1:0.1,B2:0.1):$BL3$):$BL3$,C1:0.1):$BL3$,C2:0.1):$BL3$,(D1:0.1,D2:0.1):$BL3$);
#
#	- process steps and output
#---------------------------------------------------------------------------------------------------------------------------
#		: sequence simulation based on given tree and parameter specifications
#				
#				...sequence data simulation of defined branch length combinations individually for each replication step
#
#				...output
#						:: main-output-folder: 001_indelible												
#						:: additional info files are stored in corresponding subfolders
#						:: main-output folder contains simulated raw sequences with true alignments in a separate folder 
#				
############################################################################################################################
############################################################################################################################

sub indelible{
	
	my $href_indel_setup	=	$_[0] 			; # indelible parameter									In -> defined, OUT -> unchanged
	my $sref_tree_key		=	$_[1] 			; # tree hashkey										In -> defined, OUT -> unchanged
	my $sref_tree_code		=	$_[2] 			; # tree newick incl. branch elongation code			In -> defined, OUT -> unchanged
	my $sref_ind_folder		=	$_[3] 			; # indelible result foldername							In -> defined, OUT -> unchanged
	
	print "\n\ngreaterGlider-INDELIBLE SIMULATION ", $$sref_tree_key, "..." ;
	
	############################## Result Folders and original newick tree #####################
	############################################################################################
	## generate parameter-independent indelible output folder & output specific subfolders
	my	$dir_main		=	$$sref_tree_key					;	# always
	my	$ind_main		=	$dir_main."/".$$sref_ind_folder	;	# always
	my	$ind_sub_con	=	$ind_main."/CONTROL"			;	# always
	my	$ind_sub_tre	=	$ind_main."/TREES"				;	# always
	my	$ind_sub_rat	=	$ind_main."/RATES"				;	# parameter specified, opt. later generated
	my	$ind_sub_log	=	$ind_main."/LOG"				;	# always
	my	$ind_sub_tru	=	$ind_main."/TRUE"				;	# always
	my	$ind_sub_anc	=	$ind_main."/ANCESTRAL"			;	# parameter specified, opt. later generated
	
	for my $folder ( $dir_main, $ind_main, $ind_sub_con, $ind_sub_tre, $ind_sub_log, $ind_sub_tru ){ mkdir $folder }
	##
	
	## assign specified newick tree for simulation
	my		$newick_string	=	$$sref_tree_code		;
	############################################################################################
	############################################################################################
	
	
	############################## Branch length Variables ($BL3$) #############################
	############################################################################################
	## check original newick tree for branch length variables ($BL3$)
	my	$bl3	=	qr/\$BL3\$/	;
	
	## generate hashlist if branch extension steps are coded in original newick string
	## %branchlengths_of_bl_code -> key: branch code, e.g. BL3, value: list of generated branch lengths, e.g. 0.01, 0.02...
	my	@branchlengths_of_bl_code ;
	if	(	$newick_string	=~ /$bl3/	){
		
		## assign indelible paramter settings to string variables
		my	$val_min_bl	=	$href_indel_setup->{minBL3}	;
		my	$val_max_bl	=	$href_indel_setup->{maxBL3}	;
		my	$val_stp_bl	=	$href_indel_setup->{stpBL3}	;
		
		## stop pipeline if branch length parameters are invalid
		## check start length (min) is greater zero
		## check start length is smaller as end length (max)
		## check elongation steps (stp) are greater zero
		if		( ( $val_min_bl <= 0 ) || ( $val_min_bl > $val_max_bl ) || ( $val_stp_bl <= 0 ) ){ print "\n\tParameter-ERROR in module file aINDELiblle.pm\n" and exit }
		
		## generate branch elongation numbers in @branchlengths_of_bl_code
		else	{
			
			for	( my $i = $val_min_bl;	$i <= $val_max_bl;	$i += $val_stp_bl ){ push @branchlengths_of_bl_code, $i }
		}
	}
	############################################################################################
	############################################################################################
	
	############################## write control.txt ###########################################
	############################################################################################
	# for each random seed replicate write and execute a control.txt file, considering
	# all branch length conditions
	my	$random_seed	=	$href_indel_setup->{randomseed}	- 1 ;
	my	$N_replicates	=	$href_indel_setup->{Nloops}			;
	
	for my $replicate ( 1 .. $N_replicates ){
		
		$random_seed++ ;
		
		############################################################## write control.txt
		
		##################
		## open indelible inputfile 'control.txt'
		## inputfile name 'control.txt' is mandatory
		## otherwise indelible is not running
		my	$control_file	=	"control.txt" ;
		open OUT,	">$control_file" or die "Cannot create $control_file : $!\n";
		##################
		
		##################
		## print control header to control.txt
		## not necessary needed but nice to watch
		print OUT	"////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n",
					"//                                                                                                                //\n",
					"//  INDELible V1.03 control file - control.txt                                                                    //\n",
					"//                                                                                                                //\n",
					"//                                                                                                                //\n",
					"//                                                                                                                //\n",
					"//                                                                                                                //\n",
					"////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n\n";
		##################
		
		##################
		## print indelible block parameter to control.txt
		## Block parameter have to be in the order:
		## TYPE, SETTINGS, MODEL, TREE, PARTITIONS, EVOLVE
		
		##############################################
		## Block 'TYPE'
		print OUT			"[TYPE] ", $href_indel_setup->{TYPE}, "\n" ;
		##############################################
		
		##############################################
		## Block 'SETTINGS'
		print OUT			"\n[SETTINGS]\n"	,
							"\t[randomseed]\t"	, $random_seed , "\n" ;
		
		for my $parameter	(
							'ancestralprint'			,
							'output'					,
							'printrates'				,
							'insertaslowercase'			,
							'markdeletedinsertions'		,
							'printcodonsasaminoacids'	,
							'fileperrep'				,
		)
		{ print OUT	"\t[",	$parameter,	"]\t",	$href_indel_setup->{$parameter}	,"\n" }
		##############################################
		
		##############################################
		## BLOCK 'MODEL'
		print OUT			"\n[MODEL]\t"	,	$href_indel_setup->{MODEL}		, "\n"	,
							"\t[submodel]\t",	$href_indel_setup->{submodel}	, "\n"	;
		
		for my $parameter	(
							'indelmodel'	,
							'indelrate'		,
							'geneticcode'	,
							'rates'			,
							'statefreq'		,
		)
		{ print OUT	"\t[",	$parameter,	"]\t",	$href_indel_setup->{$parameter}	,"\n" if $href_indel_setup->{$parameter} }
		##############################################
		
		##############################################
		## BLOCK 'TREE'
		my	%tree_of_treename ;
		
		## generate trees and assign tree names given defined lengths bl1, bl2, bl3
		if ( @branchlengths_of_bl_code ){
			
			## replace placeholders ($BL3$) in original newick string 
			## by BL corresponding branch length. Built treename from random seed and corresponding
			## branch length ($length), e.g. 12345_0.01 (randomseed_$length)
			## store each tree as value assigned to treename as key value in %tree_of_treename
			for	my $length (@branchlengths_of_bl_code ){
				
				(	my	$tree		=	$newick_string	)	=~	s/\$BL3\$/$length/g	;
					my	$tree_name	=	$random_seed."_".$length								;
				
				$tree_of_treename{$tree_name}	=	$tree										;
			}
			
			## print TREE block after generating trees including all coded BL length combinations
			for my $tree_name_complete	( sort keys %tree_of_treename ){
				
				print OUT	"\n[TREE]\t", $tree_name_complete, "\t", $tree_of_treename{$tree_name_complete}
			}
		}
		##
		
		## print original tree if no branch length codes are defined 
		## (trees without $bl1$, $bl2$, $bl3$ codes)
		else{
			
				my	$tree		=	$newick_string									;
			(	my	$tree_name	=	$href_indel_setup->{newick}	)	=~	s/.tre$//	;
					$tree_name	=	"tree_".$tree_name."_".$random_seed				;
			
			$tree_of_treename{$tree_name}	=	$tree ;
			print OUT	"\n[TREE]\t", $tree_name, "\t", $tree
		}
		##############################################
		
		##############################################
		## BLOCK 'PARTITIONS'
		my	$counter_partition	=	1 ;
		my	@name_partions		;
		for	my $tree_name ( sort keys %tree_of_treename ){
			
			my		$partition	=	"P".$counter_partition	; $counter_partition++ ;
			print	OUT	"\n[PARTITIONS] ", $partition, "\t[", $tree_name, "\t", $href_indel_setup->{MODEL}, "\t", $href_indel_setup->{rootlength}, "]"	;
			
			push	@name_partions, $partition
		}
		##############################################
		
		##############################################
		## BLOCK 'EVOLVE'
		print OUT	"\n\n[EVOLVE]\n"		;
		
		for	my $tree_name ( sort keys %tree_of_treename ){
			
			my	$partition	=	shift @name_partions	;
			print OUT	"\t", $partition, "\t", $href_indel_setup->{replicates} , "\t", $tree_name, "\n"	;
		}
		##############################################
		
		##################
		# close print of control.txt
		close OUT ;
		##################
		
		############################################################## execute control.txt
		
		##################
		# generate log file
		my	$logfile		=	$random_seed."_".$href_indel_setup->{log}.".txt"	;
		##################
		
		##################
		# generate system call
		my	$system_call	=	$href_indel_setup->{scr}." control.txt >".$logfile ;
		system ( $system_call ) ; 
		##################
		
		############################################################## output handling
		
		##################
		# delete space signs in fasta files
		for my $fasta ( <*.fas>  ){
			
			open	INfas, "<$fasta" or die "Cannot read file $fasta : $!\n"	;
			my		@infas	=	<INfas> ;
			close	INfas;
			unlink	$fasta ;
			
			open	OUTfas,	">$ind_main/$fasta" or die "Cannot open file $ind_main/$fasta : $!\n"	;
			for my $line ( @infas ){
				
					chomp $line ;
					
					$line =~ s/\s+//g; 
					print OUTfas $line, "\n" 
			}
			close OUTfas;
		}
		##################
		
		################## CONTROL (always)
		# move 'control' file to subfolder CONTROL
		for my $file (	<control.txt>	){ 
			
			my	$file_new	=	$ind_sub_con."/".${random_seed}."_".$file	 ;
			move ( $file, $file_new )
		}
		##################
		
		################## TREES (always)
		# move 'trees' file to subfolder TREES
		for my $file (	<trees.txt>	){ 
			
			my	$file_new	=	$ind_sub_tre."/".${random_seed}."_".$file	 ;
			move ( $file, $file_new )
		}
		##################
		
		################## RATES (optionally)
		# move 'rates' files to subfolder RATES
		# if print of sequence rates defined
		unless	(	$href_indel_setup->{ancestralprint}	eq	"FALSE"	){
			
			# generate subfolder RATES
			mkdir $ind_sub_rat	;
			
			for my $file (	<*_RATES.txt>	){
				
				my	$file_new	=	$ind_sub_rat."/".$file	 ;
				move ( $file, $file_new )
			}
		}
		##################
		
		################## LOG (optionally)
		# if defined move 'log' file to new subfolder LOG
		if	(	$href_indel_setup->{plog}	==	1	){
			
			for my $file (	<LOG.txt>	){ 
				
				my	$file_new	=	$ind_sub_log."/".${random_seed}."_".$file	 ;
				move ( $file, $file_new )
			}
		}
		
		for my $file (	<*_indelible_logfile.txt> ){
			
			my	$file_new	=	$ind_sub_log."/".$file	 ;
			move ( $file, $file_new )
		}
		##################
		
		################## ANCESTRAL (optionally)
		# if print of ancestral sequences defined
		# move files to new subfolder ANCESTRAL
		unless	(	$href_indel_setup->{ancestralprint}	eq	"FALSE"	){
			
			# generate subfolder ANCESTRAL
			mkdir $ind_sub_anc	;
			
			# move 'ancestral' files to subfolder ANCESTRAL
			for my $file (	<$ind_main/*_ANCESTRAL_*>	){ 
				
				( my	$file_new	=	$file ) =~ s/$ind_main/$ind_sub_anc/ ;
				move ( $file, $file_new )
			}
		}
		##################
		
		################## TRUE (always)
		# move 'true' files to subfolder true
		for my $file (	<$ind_main/*_TRUE_*>	){
			
			( my	$file_new	=	$file ) =~ s/$ind_main/$ind_sub_tru/ ;
			move ( $file, $file_new )
		}
		##################
	}
	############################################################################################
	############################################################################################
}
############################################################################################################################
############################################################################################################################



############################################# IQTREE Tree Reconstruction ###################################################
############################################################################################################################

############################################# Brief Description &iqtree ####################################################
############################################################################################################################
# aIQtree::iqtree : a module for ML tree reconstruction with IQTree, allowing various IQtree1 and IQtree2 parameter options
#
#	- infile data
#---------------------------------------------------------------------------------------------------------------------------
#		: multiple sequence alignment(s):	
#						:: *.fas (first choice)										
#		: partition file (yet further processing not implemented)					
#		: start tree(s) *.tre														
#							
#	- process steps and output
#---------------------------------------------------------------------------------------------------------------------------
#		: making outfile depending subresult folders for each type of output file of each iqtree analysis			
#		: best supported treefile of each analysis remains in the main outputfolder 
############################################################################################################################
############################################################################################################################
sub iqtree{
	
	my $href_iqtree_setting		=	$_[0] ; # iqtree parameter									In -> defined, OUT -> unchanged
	my $sref_tree_key			=	$_[1] ; # tree hashkey										In -> defined, OUT -> unchanged
	my $sref_iqtree_folder		=	$_[2] ; # iqtree result foldername							In -> defined, OUT -> unchanged
	my $sref_infile_msa_folder	=	$_[3] ; # infile alignment data								In -> defined, OUT -> unchanged
	
	print "\n\ngreaterGlider-APPETITE-IQTree ML reconstruction ", $$sref_tree_key, "..." ;
	
	############ IQtree system call preparation #################################
	#############################################################################
	## generate system call line based on header defined options
	my	$system_call	=	"./".$href_iqtree_setting->{scr}	; # iqtree script-name
	
	#################################
	## Substitution model and model parameter setup
	unless		(  $href_iqtree_setting->{smod} =~	/TEST|MERGE|MF|MFP/){ # without J-model test


																   $system_call	.=	" -m ".$href_iqtree_setting->{smod}		;	# substitution model
					if	( $href_iqtree_setting->{gama}	==	1	){ $system_call	.=	"+G"									}	# 1 -> +G, 0 -> without +G
					if	( $href_iqtree_setting->{inva}	==	1	){ $system_call	.=	"+I"									}	# 1 -> +I, 0 -> without +I
					if	( $href_iqtree_setting->{alpha}	>	0	){ $system_call	.=	" -a ".$href_iqtree_setting->{alpha}	}	# 0 -> alpha shape parameter value is estimated, otherwise a fixed positive float value is used
					if	( $href_iqtree_setting->{pinv}	>	0	){ $system_call	.=	" -i ".$href_iqtree_setting->{pinv}		}	# 0 -> proportion of inv sites are estimated, otherwise a fixed positive float value is used
					if	( $href_iqtree_setting->{optg}	==	1	){ $system_call	.=	" --opt-gamma-inv"						}	# 1 -> More thorough estimation for +I+G model parameters
	}
	else		{	}
	#################################
	
	#################################
	## Bootstrapping
					if	(   $href_iqtree_setting->{nboot}	>=	1															){ $system_call	.=	" -".$href_iqtree_setting->{mboot}." ".$href_iqtree_setting->{nboot}	}	# bootstrap method: b (bs+ML+Con) |bc (bs+Con)|bo (bs)|bb (ultra fast bs
					if	( ( $href_iqtree_setting->{mboot}	=~	/bb/	) && ( $href_iqtree_setting->{wbtl}	==	1	)		){ $system_call	.=	" -wbtl"																}	# Only for ultra-fast bs (-bb): Write bootstrap trees to .ufboot file (default: 0 -> none)
					if	( ( $href_iqtree_setting->{mboot}	=~	/bb/	) && ( $href_iqtree_setting->{nstb}	>=	1	)		){ $system_call	.=	" --nstep ".$href_iqtree_setting->{nstb}								}	# Only for ultra-fast bs (-bb): Iterations for UFBoot stopping rule (default: 100)
	#################################
	
	#################################
	## info output
					if	( $href_iqtree_setting->{wsrc}	==	1	){ $system_call	.=	" -wsr"									}	# 1 -> Write site rates to .rate file
					if	( $href_iqtree_setting->{wrsl}	==	1	){ $system_call	.=	" -wsl"									}	# 1 -> Write site log-likelihoods to .sitelh file
					if	( $href_iqtree_setting->{wslr}	==	1	){ $system_call	.=	" -wslr"								}	# 1 -> Write site log-likelihoods per rate category '-wslr'
					if	( $href_iqtree_setting->{wltr}	==	1	){ $system_call	.=	" -wt"									}	# 1 -> Write locally optimal trees into .treels file '-wt'
	#################################
	
	#################################
	## other options
					if	( $href_iqtree_setting->{alrt}	==	1	){ $system_call	.=	" -alrt ".$href_iqtree_setting->{alrts}	}	# 1 -> Perform SH-aLRT test: $href_iqtree_setting->{alrts}	== 0 -> Parametric aLRT test (Anisimova and Gascuel 2006); $href_iqtree_setting->{alrts}	>= 1 -> SH-like approximate likelihood ratio test (SH-aLRT) with N replicates
					if	( $href_iqtree_setting->{core}	>=	2	){ $system_call	.=	" -nt ".$href_iqtree_setting->{core}	}	# Number of cores/threads to use (REQUIRED)
					if	( $href_iqtree_setting->{root}			){ $system_call	.=	" -o ".$href_iqtree_setting->{root}		}	# outgroup taxon '-o <taxonname>', if undef -> no outgroup
					if	( $href_iqtree_setting->{pref}			){ $system_call	.=	" -pre ".$href_iqtree_setting->{pref}	}	# -pre <PREFIX> Prefix for all output files, if undef -> default: aln/partition)
					if	( $href_iqtree_setting->{quie}	==	1	){ $system_call	.=	" -quiet"								}	# Quiet mode, suppress printing to screen (stdout)
					if	( $href_iqtree_setting->{mram}	>	0	){ $system_call	.=	" -mem ".$href_iqtree_setting->{mram}	}	# Maximal RAM usage for memory saving mode
	#################################
	
	#################################
	## read in alignment and (optionally) starttree infiles
	my	( @infiles ) ;
	
	for my $indata ( <$$sref_infile_msa_folder/*.fas> ){ push @infiles, $indata }
	#################################
	
	#################################
	## prepare final system command lines and process iqtree ML reconstruction
	for my $msa_infile ( <$$sref_infile_msa_folder/*.fas> ){

		my $system_call_final = $system_call ;

		## assign msa infile to system command
		$system_call_final	.=	" -s ". $msa_infile ;
		##

		## assign logfile to system command
		(	my $logfile		 =	$msa_infile	)	=~	s/.fas$// ;
			$logfile		.=	"_".$href_iqtree_setting->{log}.".txt" ;
		if	( $href_iqtree_setting->{plog} == 1 ){	$system_call_final	.=	" >".$logfile	}
		##

		## store system call line
		print	"\n...", $msa_infile ;
		system ( $system_call_final ) ; 
		##
	}
	#############################################################################
	#############################################################################


	############## Output file handling #########################################
	#############################################################################
	## generate iqtree result folder
															mkdir $$sref_iqtree_folder				;
															mkdir $$sref_iqtree_folder."/mldist"	;
															mkdir $$sref_iqtree_folder."/log"		;
															mkdir $$sref_iqtree_folder."/info"		;
															mkdir $$sref_iqtree_folder."/contree"	;
	if	( $href_iqtree_setting->{wsrc}	==	1			){	mkdir $$sref_iqtree_folder."/rate"		}
	if	( $href_iqtree_setting->{nboot}	>=	1			){	mkdir $$sref_iqtree_folder."/boot"		}
	if	( $href_iqtree_setting->{wltr}	==	1			){	mkdir $$sref_iqtree_folder."/treels"	}
	if	( $href_iqtree_setting->{wrsl}	==	1			){	mkdir $$sref_iqtree_folder."/sitelh"	}
	#################################
	
	#################################
	## define list variables for r-ggtree plots
	my (
			@ml_trees_best	,	# store best phyml trees for graphical output preparation
			@ml_trees_boot	,	# store bootstrapped phyml trees for graphical output preparation
			@ml_trees_bmsa	,	# store best tree corresponding alignment
			@ml_start_tree	,	# store bionj generated start trees
		) ;
	#################################
	
	#################################
	## iqtree result file handling
	for	my $msa_infile	(	@infiles	){
		
		#################################
		## move iqtree results from alignment input folder to iqtree output result folder
		(	my	$outfile_new	= $msa_infile	)	=~	s/\.fas// ;	push @ml_trees_bmsa, $msa_infile ;

		(	my	$new_path_tree	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder/			;
		(	my	$new_path_cotr	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/contree/	;
		(	my	$new_path_sttr	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/bionj/		;
		(	my	$new_path_iqtr	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder/			;
		(	my	$new_path_boot	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/boot/		;
		(	my	$new_path_rate	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/rate/		;
		(	my	$new_path_trls	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/treels/	;
		(	my	$new_path_silh	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/sitelh/	;
		(	my	$new_path_mldi	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/mldist/	;
		(	my	$new_path_iqlo	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/log/		;
		(	my	$new_path_logs	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/log/		;
		(	my	$new_path_info	= $outfile_new	)	=~	s/$$sref_infile_msa_folder/$$sref_iqtree_folder\/info/		;

			my	$outfile_logs					=	$outfile_new."_".$href_iqtree_setting->{log}.".txt"		;
			my	$outfile_boot					=	$msa_infile.".boottrees"								;
			my	$outfile_rate					=	$msa_infile.".rate"										;
			my	$outfile_trls					=	$msa_infile.".treels"									;
			my	$outfile_silh					=	$msa_infile.".sitelh"									;
			my	$outfile_sttr					=	$msa_infile.".bionj"									;
			my	$outfile_iqtr					=	$msa_infile.".iqtree"									;
			my	$outfile_cotr					=	$msa_infile.".contree"									;
			my	$outfile_tree					=	$msa_infile.".treefile"									;
			my	$outfile_mldi					=	$msa_infile.".mldist"									;
			my	$outfile_info					=	$msa_infile.".ckp"										;
			my	$outfile_zipf					=	$msa_infile.".ckp.gz"									;
			my	$outfile_iqlo					=	$msa_infile.".log"										;

				$new_path_logs					=	$new_path_logs."_".$href_iqtree_setting->{log}.".txt"	;
				$new_path_boot					=	$new_path_boot."_iqtree_boottrees.tre"					;	push @ml_trees_boot	,	$new_path_boot	;
				$new_path_rate					=	$new_path_rate."_iqtree_rate.txt"						;
				$new_path_trls					=	$new_path_trls."_iqtree_treels.txt"						;
				$new_path_silh					=	$new_path_silh."_iqtree_sitelh.txt"						;
				$new_path_sttr					=	$new_path_sttr."_iqtree_bionj.tre"						;	push @ml_start_tree	,	$new_path_sttr	;
				$new_path_iqtr					=	$new_path_iqtr."_iqtree_results.txt"					;
				$new_path_cotr					=	$new_path_cotr."_iqtree_contree.tre"					;
				$new_path_tree					=	$new_path_tree."_iqtree_treefile.tre"					;	push @ml_trees_best	,	$new_path_tree	;
				$new_path_mldi					=	$new_path_mldi."_iqtree_mldistance_matrix.txt"			;
				$new_path_info					=	$new_path_info."_iqtree_info.txt"						;
				$new_path_iqlo					=	$new_path_iqlo."_iqtree_generated_log.txt"				;

		move	(	$outfile_iqtr	,	$new_path_iqtr	) ;
		move	(	$outfile_cotr	,	$new_path_cotr	) ;
		move	(	$outfile_tree	,	$new_path_tree	) ;
		move	(	$outfile_sttr	,	$new_path_sttr	) ;
		move	(	$outfile_mldi	,	$new_path_mldi	) ;
		move	(	$outfile_info	,	$new_path_info	) ;
		move	(	$outfile_iqlo	,	$new_path_iqlo	) ;
		
		############
		# Try to avoid error message Copy failed: No such file or directory at" -> possible error -> too fast
		# first, testing if file exists and include a short loop for sleep and testing again
		if	( $href_iqtree_setting->{nboot}	>=	1	){
			
			&bool (
					\$outfile_boot	,	# old filename
					\$new_path_boot	,	# new filename
			)
		}
		
		if	( $href_iqtree_setting->{wsrc}	>=	1	){
			
			&bool (
					\$outfile_rate	,	# old filename
					\$new_path_rate	,	# new filename
			)
		}
		
		if	( $href_iqtree_setting->{wltr}	==	1	){
			
			&bool (
					\$outfile_trls	,	# old filename
					\$new_path_trls	,	# new filename
			)
		}
		
		if	( $href_iqtree_setting->{wrsl}	==	1	){
			
			&bool (
					\$outfile_silh	,	# old filename
					\$new_path_silh	,	# new filename
			)
		}
		
		if	( $href_iqtree_setting->{plog}	==	1	){
			
			&bool (
					\$outfile_logs	,	# old filename
					\$new_path_logs	,	# new filename
			)
		}
		
		sub bool{
			
			my $sref_old_path	= $_[0] ; # old filename	In -> defined, OUT -> unchanged
			my $sref_new_path	= $_[1] ; # new filename	In -> defined, OUT -> unchanged
			
			my $bool	= -e $$sref_old_path ;
			my $counter	= 0 ;
			
			while( $bool == 0 ){
				
				sleep(1) ; ++$counter		;
				$bool = -e $$sref_old_path	;
				last if $counter > 15		;
			}
			if ( -e $$sref_old_path ){ move ( $$sref_old_path, $$sref_new_path	) or die "Copy failed: $!" }
		}
		#################################
		
		#################################
		# delete remaining output files
		unlink	$outfile_zipf
		#################################
	}
	#############################################################################
	#############################################################################
}
############################################################################################################################
############################################################################################################################



############################################# IQTREE2 Concordance Measure ##################################################
############################################################################################################################

############################################# Brief Description &concordance ###############################################
############################################################################################################################
# a subroutine for concordance measure given an alignment(s) and a reference tree
#
#	- Description & Commands:
#---------------------------------------------------------------------------------------------------------------------------
#		This module allows concordance measures with IQtree2 given a single or multiple reference tree(s) and an alignment. 
#		Currently, only the site-concordance measure ('cf_type' => 's') is implemented (yet, not enough time to do the gene-concordance stuff 'cf_type' => 'g')
#		Different options are possible:
#
#			: a single alignment (raw, simulated) and a single reference tree (in best case in iqtree newick format -> allows further data processings due to individual substring support), cf_parameter setting:
#				'msa'		=> 'name.fas'
#				'rtree'		=> 'name.tre'
#				'procscf'	=> 'a'
#			: a set of alignments with each alignment compared to a filename matching reference tree (can also be part of a set of reference trees), cf_parameter setting:
#				'msa'		=> 'a'
#				'rtree'		=> 'a'
#				'procscf'	=> 'm'
#			: a set of alignments with each alignment compared to a set of reference trees, cf_parameter setting:
#				'msa'		=> 'a'
#				'rtree'		=> 'a'
#				'procscf'	=> 'a'
#
#		General output files ('sim' = 1 & 'sim' = 0)
#			: <cf_subtree_summary.txt>								: list of internal branch individual cf support
#			: <cf_tree_classified_summary.txt>						: list of cf supported reference tree(s) without branch lengths
#			: <cf_tree_individual_summary.txt>						: list of cf supported reference tree(s) with branch lengths
#
#		Additional (branch length BL3 related) output. If 'sim' = 1, print of alignment individual 
#		and (BL3 depending) mean identified CF support of internal branch relationships related to different lengths BL3
#
#			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! currently only tested for IQtree formatted newick strings !
#			: <mean_complTree_cf_support.tsv>	: a tsv table file listing for each full tree and length bl3 the corresponding average CF support
#				'sim' = 1
#			: <mean_subtree_cf_support.tsv>		: a tsv table file listing for each subtree and length bl3 the corresponding average CF support
#				'sim' = 1
#			: <mean_subtree_cf_support.txt>		: a txt table file listing for each subtree and length bl3 the corresponding average CF support
#				'sim' = 1
#			: <mean_complTree_cf_support_table.tex> and <mean_complTree_cf_support_table.pdf>	: listing for each full tree and length bl3 the corresponding average CF support
#			  !! Newick substrings are listed in the first column, be aware that latex table can overflow if newickstrings are long!!
#				'sim' = 1
#			: <mean_subtree_cf_support_table.tex> and <mean_subtree_cf_support_table.pdf>		: listing for each subtree and length bl3 the corresponding average CF support
#			  !! Newick substrings are listed in the first column, be aware that latex table can overflow if newickstrings are long!!
#				'sim' = 1
#			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#	- infile data
#---------------------------------------------------------------------------------------------------------------------------
#		: multiple sequence alignment(s):
#						:: *.fas
#		: reference tree(s) in newick format (optionally '$BL3$' coded):
#						:: *.tre
############################################################################################################################
############################################################################################################################
sub concordance{
	
	my $href_concordance_setting		=	$_[0] ; # iqtree parameter									In -> defined, OUT -> unchanged
	my $sref_tree_key					=	$_[1] ; # tree hashkey										In -> defined, OUT -> unchanged
	my $sref_concordance_folder			=	$_[2] ; # resultfolder of concordance measure				In -> defined, OUT -> unchanged
	my $sref_infile_data_folder			=	$_[3] ; # infile alignment or gene trees					In -> defined, OUT -> unchanged
	my $sref_infile_rtree_folder		=	$_[4] ; # infile reference tree data						In -> defined, OUT -> unchanged
	
	print "\n\ngreaterGlider-APPETITE-Concordance Measure ", $$sref_tree_key, "..." ;
	
	####################################
	####################################
	# print appetite parameter specifications
	# for my $parameter ( sort keys %$href_concordance_setting ){ print "\nparameter: ", $parameter, "\tvalue: ", $href_concordance_setting->{$parameter} } exit ;
	####################################
	
	####################################
	####################################
	# generate main directory and subfolder due to parameter speicifactions
	mkdir $$sref_concordance_folder ;
	####################################
	
	
	################# Read Reference Tree(s) ####################################
	#############################################################################
	# sampling of reference tree(s) in %ref_trees
	my %ref_trees	; # key1: reftree filename (with path); value: number of occurence
	
	#############
	# if multiple reference trees have to be analysed
	if		( $href_concordance_setting->{rtree} =~ /a/ ){
		
		for	my $ref_tree ( <$$sref_infile_rtree_folder/*.tre> ){ $ref_trees{$ref_tree}++ }
	}
	#############
	
	#############
	# if a single reference tree is specified
	elsif 	( $href_concordance_setting->{rtree} =~ /.tre$/ ){
		
		my	$single_reftree = $$sref_infile_rtree_folder."/".$href_concordance_setting->{rtree} ;
			$ref_trees{$single_reftree}++
	}
	#############
	
	#############
	else{ die "\nError in aIQtree::concordance, no reference tree (*.tre) in ", $$sref_infile_rtree_folder, "/\n" }
	#for my $ref_tree ( sort keys %ref_trees ){ print "\nref tree: ", $ref_tree, "\t", $ref_trees{$ref_tree} } exit;
	#############################################################################
	#############################################################################
	
	
	################ Concordance measure system call preparation & execution ####
	############################################################################# 
	## generate system call line based on header defined options
	my	$system_call	=	"./".$href_concordance_setting->{scr}	; # iqtree script-name
	if	( $href_concordance_setting->{dftree}		==	1	){ $system_call	.=	" --df-tree"	}		# Write discordant trees associated with gDF1
	if	( $href_concordance_setting->{cfverbose}	==	1	){ $system_call	.=	" --cf-verbose"	}		# Write CF per tree/locus to cf.stat_tree/_loci
	
	my	@system_calls	;
	my	%infile_of_call	;
	
	## Parameter check of defined CF measure ################
	#########################################################
	unless ( $href_concordance_setting->{cf_type} =~ /^s$|^g$/i ) { die "\nERROR: Unknown parameter for concordance measure 'cf_type: ", $href_concordance_setting->{cf_type},"'\n\n" }
	#########################################################
	
	################ site concordance measure (sCF) #########
	#########################################################
	# Site-concordance measure (sCF) is based on alignment and reference tree comparison
	if	(( $href_concordance_setting->{cf_type} =~	/^s$/i ) && ( $href_concordance_setting->{procsys} == 1 ) ){
		
		#########################################################
		# extend system command to general scf depending commands
		if	( $href_concordance_setting->{cfquartet}	==	1	){ $system_call	.=	" --cf-quartet"									}	# Write sCF for all resampled quartets to .cf.quartet
		if	( $href_concordance_setting->{scfNUM}		>=	1	){ $system_call	.=	" --scf ".$href_concordance_setting->{scfNUM}	}	# Number of quartets for site concordance factor (sCF)
		#########################################################
		
		############ Read alignment(s) ##########################
		#########################################################
		# assign input alignments to @alignments
		my	@alignments ; 
		
		######################
		# if $href_concordance_setting->{msa} =~ /a/: assign all inpath fasta (.fas) alignments
		if ( $href_concordance_setting->{msa} =~ /a/ ){
			
			for my	$fas ( <$$sref_infile_data_folder/*.fas> ){ push @alignments, $fas	}
		}
		######################
		
		######################
		# Otherwise assign the single defined fasta (.fas) alignment of inpath to @alignments
		elsif ( $href_concordance_setting->{msa} =~ /.fas$/ ){
			
			$alignments[0] = $$sref_infile_data_folder."/".$href_concordance_setting->{msa} 
		}
		######################
		
		######################
		else { die "\nError in aIQtree::concordance, no alignment in ", $$sref_infile_data_folder, "/\n" }
		#########################################################
			
		############# For each alignment... #####################
		#########################################################
		MSA: for my	$msa ( sort @alignments ){
				
			####################################
			####################################
			# complete system call due to each alignment individually ($system_call_final)
			my	$system_call_final = $system_call ;
			####################################
			
			
			####################################
			####################################
			# substitute inpath prefix and file suffix (.fas) for name comparison between $msa_raw_name and pot. ref trees
			( my $msa_raw_name	=	$msa )	=~	s/^$$sref_infile_data_folder\/|.fas$//g ;	# print "\ncompl: ", $msa, "\tno path: ", $msa_raw_name ;
			####################################
			
			
			###### msa comparison with msa-name matching ref tree ###
			#########################################################
			# analyse each alignment once given a name matching reftree (one CF analysis for each alignment specified by parameter procsf='m')
			# if match has been found, process system call and output handling
			# Afterwards proceed with next alignment
			if		( $href_concordance_setting->{procscf}	=~ /m/	){ #print "\nwrong!!!!\n"; exit;
					
				#################################### assign filename match between msa & reftree & process execution
				# check name match between $msa_raw_name 
				# and potential ref. trees 
				for my $rtree ( %ref_trees ){
					
					#print "\n\nreftree: $rtree\nmsa: $msa_raw_name\n\n"; exit;
					
					if ( $rtree =~ /$msa_raw_name/ ){
						
						$system_call_final 	.= " -t ".$rtree." -s ".$msa 		; #print "\nmatch: ", $rtree, "\t", $msa_raw_name ;
						
						system ($system_call_final )
					}
				}
				####################################
				
				#################################### Output handling
				# rename output files by adding msa filename inbetween reftree filename (prefix) and cf output codes (suffix)
				# adding the msa filename avoids file overwriting, because cf outputs are only based on the reftree filename
				# thus multiple alignment tests with the same reftree are not overwritten
				my	$split_pattern		=	"tre.cf"					;	# pattern on which the output string is splitted
				my	$new_string_pattern	=	"tre.".$msa_raw_name."_msa.cf"	;	# new pattern string, replacing $split_pattern in cf output filenames
				
				&rename_cf_output(
					
					\$$sref_infile_rtree_folder	,	# path to located cf output files, usually the reftree input path	In -> defined, OUT -> unchanged
					\$split_pattern				,	# filename pattern replaced by $new_string_pattern					In -> defined, OUT -> unchanged
					\$new_string_pattern		,	# string to add to cf output files									In -> defined, OUT -> unchanged
				) ;
				####################################
			}
			#########################################################
				
			
			######## msa comparison with each reftree ###############	
			#########################################################
			# analyse each alignment with each reference tree (multiple CF analyses for each msa if number ref trees > 1)
			elsif	( $href_concordance_setting->{procscf}	=~ /a/ ){ #print "\ncorrect!!!!\n"; exit;
				
				####################################
				# check name match between $msa_raw_name and potential ref. trees 
				for my $rtree ( keys %ref_trees ){ print "\nrtree (a): ", $rtree,"\n";
						
					$system_call_final .= " -t ".$rtree." -s ".$msa			;
					
					system ($system_call_final)  ;
					
					$system_call_final = $system_call ;
				}
				####################################
				#exit;
				####################################
				# rename output files by adding msa filename inbetween reftree filename (prefix) and cf output codes (suffix)
				# adding the msa filename avoids file overwriting, because cf outputs are only based on the reftree filename
				# thus multiple alignment tests with the same reftree are not overwritten
				my	$split_pattern		=	"tre.cf"						;	# pattern on which the output string is splitted
				my	$new_string_pattern	=	"tre.".$msa_raw_name."_msa.cf"	;	# new pattern string, replacing $split_pattern in cf output filenames
				
				&rename_cf_output(
					
					\$$sref_infile_rtree_folder	,	# path to located cf output files, usually the reftree input path	In -> defined, OUT -> unchanged
					\$split_pattern				,	# filename pattern replaced by $new_string_pattern					In -> defined, OUT -> unchanged
					\$new_string_pattern		,	# string to add to cf output files									In -> defined, OUT -> unchanged
				) ;
			}
			#########################################################
		} 
		#########################################################
	}
	#########################################################
	#########################################################
	
	
	################ gene concordance measure (gCF) #########
	#########################################################
	# Gene-concordance measure (gCF) is based on multiple intrees and reference tree comparison
	elsif( $href_concordance_setting->{cf_type} =~	/^g$/i ){ die "\nERROR: gene-concordance measure is currently not implemented in aIQtree::concordance\n\n" }
	############################################################################# 
	############################################################################# 
	
	
	################ Sorting of general output files ############################
	############################################################################# 
	# generate outputpath related subfolder for each type of result, which is 
	# printed by iqtree2 independently of sCF or gCF calculation and move 
	# corresponding result file from reftree inpath to new subfolder
	my $subfolder_cf_tree = $$sref_concordance_folder."/cf_tree"	; mkdir $subfolder_cf_tree ;
			
	for my $cftree (<$$sref_infile_rtree_folder/*.cf.tree>){
		
		( my $cftree_path = $cftree ) =~ s/$$sref_infile_rtree_folder/$subfolder_cf_tree/ ;
		move ( $cftree,	$cftree_path	) ;
	}
	#########################################
	
	#########################################
	my $subfolder = $$sref_concordance_folder."/cf_tree_nex" 		; mkdir $subfolder ;			# includes...
	
	for my $cftreenex (<$$sref_infile_rtree_folder/*.cf.tree.nex>){
				
		( my $cftreenex_path = $cftreenex ) =~ s/$$sref_infile_rtree_folder/$subfolder/ ;
		move ( $cftreenex,	$cftreenex_path	)
	}
	#########################################
	
	#########################################
	$subfolder = $$sref_concordance_folder."/cf_stat" 				; mkdir $subfolder ;			# includes...
	
	for my $cfstat (<$$sref_infile_rtree_folder/*.cf.stat>){
				
		( my $cfstat_path = $cfstat ) =~ s/$$sref_infile_rtree_folder/$subfolder/ ;
		move ( $cfstat,	$cfstat_path	)
	}
	#########################################
	
	#########################################
	$subfolder = $$sref_concordance_folder."/cf_branch" 			; mkdir $subfolder ;			# includes...
	
	for my $cfbranch (<$$sref_infile_rtree_folder/*.cf.branch>){
				
		( my $cfbranch_path = $cfbranch ) =~ s/$$sref_infile_rtree_folder/$subfolder/ ;
		move ( $cfbranch,	$cfbranch_path	)
	}
	#########################################
	
	#########################################
	$subfolder = $$sref_concordance_folder."/log" 					; mkdir $subfolder ;			# includes...
	
	for my $log (<$$sref_infile_rtree_folder/*.tre.log>){
				
		( my $log_path = $log ) =~ s/$$sref_infile_rtree_folder/$subfolder/ ;
		move ( $log,	$log_path	)
	}
	#############################################################################
	#############################################################################
	
	
	################ Sorting of sCF specific output files #######################
	#############################################################################
	# generate outputpath related subfolder for each type of result, which is 
	# printed by iqtree2 based on sCF calculation and move 
	# corresponding result file from reftree inpath to new subfolder
	if	( $href_concordance_setting->{cf_type} =~	/^s$/i ){
		
		#########################################
		if	( $href_concordance_setting->{cfverbose}	==	1	){
			
			my $subfolder = $$sref_concordance_folder."/cf_stat_tree" ; mkdir $subfolder ;		# includes CF per tree/locus to cf.stat_tree/_loci
			
			for my $verbose (<$$sref_infile_rtree_folder/*.stat_tree>){
				
				( my $verbose_path = $verbose ) =~ s/$$sref_infile_rtree_folder/$subfolder/ ;
				move ( $verbose,	$verbose_path	)
			}
		}
		#########################################
		
		#########################################
		if	( $href_concordance_setting->{cfquartet}	==	1	){
			
			my $subfolder = $$sref_concordance_folder."/cf_quartet"	; mkdir $subfolder ;		# includes sCF for all resampled quartets to .cf.quartet
			
			for my $cfquartet (<$$sref_infile_rtree_folder/*.quartet>){
				
				( my $cfquartet_path = $cfquartet ) =~ s/$$sref_infile_rtree_folder/$subfolder/ ;
				move ( $cfquartet,	$cfquartet_path	)
			}
		}
		#########################################
	}
	#############################################################################
	#############################################################################
	
	
	################ Sorting of gCF specific output files #######################
	#############################################################################
	# gCF MEASURE NOT IMPLEMENTED YET
	#############################################################################
	#############################################################################
	
	
	################ Further processing of scf/gcf support results ##############
	#############################################################################
	# summarization and further processing of cf evaluated reference trees
	# output includes tables in .tex, .txt, and .tsv format
	# as well as optionally .pdf R plotting
	&summarize_cf_trees(
		
		\$subfolder_cf_tree						,	# path to cf supported reftrees														In -> defined, OUT -> unchanged
		\$$sref_concordance_folder				,	# path to new cf-tree summarized result prints										In -> defined, OUT -> unchanged
		\$href_concordance_setting->{sim}		,	# flag for further result print in respect of different lengths bl2					In -> defined, OUT -> unchanged
		\$href_concordance_setting->{cf_type}	,	# specified cf type (scf or gcf measure) impedes further filehandling 				In -> defined, OUT -> unchanged
		\$href_concordance_setting->{rgr}		,	# if 1 -> R plot, otherwise no R plot												In -> defined, OUT -> unchanged
		\$href_concordance_setting->{subnwk}	,	# if 1 -> newick substrings without inner brackets, 0 -> with inner brackets		In -> defined, OUT -> unchanged
		\$href_concordance_setting->{prtnwkcf}	,	# defined cf-tree layout for ggtree print											In -> defined, OUT -> unchanged
		\$href_concordance_setting->{prtnwkraw}	,	# defined raw-tree layout for ggtree print											In -> defined, OUT -> unchanged
	) ;
	#############################################################################
	#############################################################################
	
	################# subroutines in &concordance ############################################################################################################
	##########################################################################################################################################################
	sub summarize_cf_trees{
		
		############################################################################# 
		# called from &concordance (main routine)
		#
		# summarization and further processing of cf evaluated reference trees
		# output includes tables in .tex, .txt, and .tsv format
		#
		# uses following subroutines:
		#	- &structure_handling
		#	- &bl2_subtree_analysis
		############################################################################# 
		
		my	$sref_infile_path			= $_[0]	; # path to cf supported reftrees													In -> defined, OUT -> unchanged
		my	$sref_outfile_path			= $_[1]	; # path to new cf-tree summarized result prints									In -> defined, OUT -> unchanged
		my	$sref_bl2_output_flag		= $_[2]	; # flag for further result print in respect of different lengths bl2				In -> defined, OUT -> unchanged
		my	$sref_cf_type_setting		= $_[3]	; # specified cf type (scf or gcf measure) impedes further filehandling 			In -> defined, OUT -> unchanged
		my	$sref_rplot_param			= $_[4]	; # if 1 -> R plot, otherwise no R plot												In -> defined, OUT -> unchanged
		my	$sref_subnwk_parameter		= $_[5]	; # if 1 -> newick substrings without inner brackets, 0 -> with inner brackets		In -> defined, OUT -> unchanged
		my	$sref_prtnwkcf_parameter	= $_[6]	; # defined cf-tree layout for ggtree print											In -> defined, OUT -> unchanged
		my	$sref_prtnwkraw_parameter	= $_[7]	; # defined raw-tree layout for ggtree print										In -> defined, OUT -> unchanged
		
		#########################################
		# Open result summary file <cf_tree_individual_summary.txt>
		# to print 'Name_cf_tree_support_resultfile \t Newick_string_with_support_and_branch_lengths'
		# for each result file in a separate line
		my		$tree_with_cf_summary_complete	= $$sref_outfile_path."/cf_tree_individual_summary.txt" ;
		open	OUTi, 	">", $tree_with_cf_summary_complete ;
		print	OUTi	"Data\tIndividual Node Support (sCF)\tAverage Tree Support (sCF)\n" ;
		#########################################
		
		#########################################
		# print each CF analysed reference tree...
		#  1) without branch lengths but with CF scores
		#  2) without branch lengths and without cf score
		# tabstop delimited wth corresponding filename (all together in a new line)
		# to new main output path file cf_tree_classified_summary.txt
		my	$tree_with_cf_summary_reduced	= $$sref_outfile_path."/cf_tree_classified_summary.txt" ;
		open 	OUTc,	">", $tree_with_cf_summary_reduced ;
		print	OUTc	"Data\tIndividual Node Support (sCF)\tRaw Newick\tAverage Tree Support (sCF)\n" ;
		#########################################
		
		#########################################
		# for each scf analysed tree file, extract...
		# 1) ...and print the full output newick string (incl branch lengths and cf scores): cf_tree_individual_summary.txt
		# 2) ...and print the full output newick string without branch lengths, but with cf scores: cf_tree_classified_summary.txt
		# 3) ...and print the full output newick string without branch lengths and without cf scores: cf_tree_classified_summary.txt
		# 4) newick string related subtrees and corresponding cf score (currently only for iqtree analysed newick format)
		my	%cf_of_subtree_of_treefile	;	# key1: tree filename; key2: subtree string; value: subtree related cf score
		my	%cf_of_cpltree_of_treefile	;	# key1: tree filename; key2: raw nwk string of complete tree; value: tree related overall cf score
		my	%seen_subtree				;	# key1: subtree string; value: number of occurence given all tree files
		my	%seen_nwk_raw				;	# key1: complete tree string without CF and branch lengths; value: number of occurence given all tree files
		my	%seen_iqtree				;	# flag for iqtree formatted newick trees. a positive flag is allowing subtree analyses 
		my	@cf_trees					;	# list of cf evaluated trees for ggtree output (if specified via ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		
		for my $cf_tree (<$$sref_infile_path/*.cf.tree>){
			
			#########################################
			# sample cf supported newick files for ggtree
			push @cf_trees, $cf_tree ;
			#########################################
			
			#########################################
			# extract actual treefilename from inputpath prefix
			( my $cftree_raw = $cf_tree ) =~ s/${$sref_infile_path}\/// ;
			#########################################
			
			#########################################
			# read newick string of reference treefile (includes branch lengths and cf scores)
			open	IN, "<", $cf_tree ;
			chomp ( my $newick_line = <IN> ) ;
			close	IN ;
			#########################################
			
			#########################################
			# substitute branch lengths of original newick string
			( my $tree_without_brlengths	=	$newick_line			) =~ s/:\d\.\d+//g ;
			#########################################
			
			#########################################
			# substitute cf scores in newick without branch lengths
			( my $tree_raw					=	$tree_without_brlengths ) =~ s/\)\d+(\.\d+)?/\)/g ;
			#########################################
			
			#########################################
			# calculate average CF score of all internal node support
			( my $nwk_cf_raw				=	$tree_without_brlengths ) =~ s/\)|,\(/:::/g ;
			  my @nwk_prts					=	split ":::", $nwk_cf_raw ;
			  my $cf_total					=	0 ;
			  my $cf_counter				=	0 ;
			
			for my $nwk_prt ( @nwk_prts ){
				
				if ( $nwk_prt =~ /^\d+(\.\d+)?$/ ){ $cf_total += $nwk_prt ; $cf_counter++ }
			}
			
			  my $avrg_cf					=	sprintf("%.1f", ( $cf_total / $cf_counter ) );
				 
			$cf_of_cpltree_of_treefile{$cf_tree}{$tree_raw} = $avrg_cf ;
			$seen_nwk_raw{$tree_raw}++ ;
			#########################################
			
			#########################################
			# print tree filename and complete newick string to 
			# cf_tree_individual_summary.txt
			print OUTi $cftree_raw, "\t", $newick_line, "\t", $avrg_cf, "\n" ;
			#########################################
			
			#########################################
			# print tree filename, newick without branch lengths, and raw newick tabstop delimited
			# to cf_tree_classified_summary.txt
			print OUTc $cftree_raw, "\t", $tree_without_brlengths, "\t", $tree_raw, "\t", $avrg_cf, "\n" ;
			#########################################
			
			#########################################
			# further ćf support extraction due to individual subtree support
			# currently only tested for iqtree nested newick bracketstring formatted trees 
			if ( $cf_tree =~ /_iqtree_/ ){
				
				#########################################
				# store trefile in seen_iqtree for further bl2 depending 
				# subtree depending cf analyses (if parameter 'sim==1')
				$seen_iqtree{$cf_tree}++ ;
				#########################################
				
				#########################################
				# newick string decomposition into all internal node classified
				# subtree relationships and assignment of subtree related cf score
				# @substring_info stores final substrings and corresponding cf score as separate elements, 
				# e.g. qw/(A1,A2), 99.5, (B1,B2), 99.6, .../ 
				my @substring_info = &structure_handling( \$tree_without_brlengths ) ;
				#########################################
				
				####################
				# terminal test print
				# print "\n", $tree_without_brlengths; for ( @substring_info ){ print "\n", $_ } print "\n" ; #exit ;
				####################
				
				#########################################
				# for each treefile (key1) hash assignment of 
				# each substring (key2) and corresponding cf score (value)
				while ( @substring_info ){
					
					my $subtree	=	shift @substring_info ; 					# extract newick subtree
					my $cfscore	=	shift @substring_info ; 					# extract subtree corresponding cf score
					
					####################
					# substitution of inner bracket strutures in given subtree
					# subtree assigned CF score does not refer to inner bracket relationships
					# sorting of inner sequence names avoids multiple prints of the same set of ingroup sequences
					if ( $$sref_subnwk_parameter == 1 ){ 
						
						$subtree		=~	s/\(|\)//g	;
						my	@seqs		=	split	",", $subtree ;
							@seqs		=	sort		 @seqs ;
							$subtree	=	join	",", @seqs ;
							$subtree	=	"(".$subtree.")" ;
					}
					####################
					
					$seen_subtree{$subtree}++ ;									# count subtrees overall tree files
					$cf_of_subtree_of_treefile{$cf_tree}{$subtree} = $cfscore ; # hash assignment
				}
				#########################################
			}
			#########################################
		}
		#########################################
		
		#########################################
		# close filehandle
		close OUTi ; # cf_tree_individual_summary.txt
		close OUTc ; # cf_tree_classified_summary.txt
		#########################################		
		
		#########################################
		# print & optionally plotting of newick string related subtrees and corresponding cf score (currently only for iqtree analysed newick format)
		# for each tree file, all substrings and related cf scores are printed separately in a single line (tabstop delimited)
		# printed to new cf_subtree_summary.txt file in main output path
		if ( %seen_iqtree ){
			
			#########################################
			# for each tree file, sampling of all treefile related substrings 
			# and cf scores combined to a single string in @datalines
			my	@datalines ;
			for my $treefile ( sort keys %cf_of_subtree_of_treefile ){
				
				####################
				# begin data string with tree filename
				my	$dataline = $treefile ;
				####################
				
				####################
				# check for each of the overall treefiles identfied substrings
				# if substring is cf suported in actual tree file.
				# if true, add substring and cf score to data string
				# otherwise, add only tabstop delimiters (keep both string infos empty)
				# this guarantees that similar subtrees of different treefiles are listed in the same tabstop order
				for my $substr ( sort keys %seen_subtree ){
					
					if 	( $cf_of_subtree_of_treefile{$treefile}{$substr} )	{ $dataline .= "\t".$substr."\t".$cf_of_subtree_of_treefile{$treefile}{$substr} }
				}
				####################
				
				####################
				# store complete data line of actual treefile
				push @datalines, $dataline ;
				####################
			}
			#########################################
			
			#########################################
			# join treefile individual data strings by combining them with a newline
			my	$data_out					= join "\n", @datalines ;
			#########################################
			
			#########################################
			# print combined data-string to cf_subtree_summary.txt
			my	$general_subtree_summary	= $$sref_outfile_path."/cf_subtree_summary.txt" ;
			open  OUTd , ">", $general_subtree_summary ;
			print OUTd	 $data_out ;
			close OUTd ;
			#########################################
			
			#########################################
			# print and (optionally) plot single and mean CF support of internal branch relationships
			# related to different lengths BL3
			if ( $$sref_bl2_output_flag == 1 ){
				
				&bl2_subtree_analysis(
					
					\%seen_iqtree				, # as iqtree identified tree files												In -> defined, OUT -> unchanged
					\%seen_subtree				, # key1: tree filename; key2: subtree string; value: subtree related cf score	In -> defined, OUT -> unchanged
					\%seen_nwk_raw				, # key1: complete tree string without CF and branch lengths; value: number		In -> defined, OUT -> unchanged
					\%cf_of_subtree_of_treefile	, # key1: tree filename; key2: subtree string; value: subtree related cf score	In -> defined, OUT -> unchanged
					\%cf_of_cpltree_of_treefile	, # key1: tree filename; key2: raw nwk string; value: total cf score			In -> defined, OUT -> unchanged
					\$$sref_outfile_path		, # path to new cf-tree summarized result prints								In -> defined, OUT -> unchanged
					\$$sref_cf_type_setting		, # specified cf type (scf or gcf measure) impedes further filehandling 		In -> defined, OUT -> unchanged
					\$$sref_rplot_param			, # if 1 -> R plot, otherwise no R plot											In -> defined, OUT -> unchanged
				);
			}
			#########################################
		}
		#########################################
		
		
	}
	#####

	#####
	sub bl2_subtree_analysis{
		
		############################################################################# 
		# called from &summarize_cf_trees
		#
		# print and (optionally) plot single and mean CF support of internal branch 
		# relationships related to different lengths BL3
		#
		# uses following modules/subroutines:
		#	- aRPLOTS::r_lineplot
		#	- &latex_table_simple
		############################################################################# 
		
		my	$href_seen_iqtree					= $_[0]	; # as iqtree identified tree files												In -> defined, OUT -> unchanged
		my	$href_seen_subtree					= $_[1]	; # key1: tree filename; key2: subtree string; value: subtree related cf score	In -> defined, OUT -> unchanged
		my	$href_seen_nwk_raw					= $_[2]	; # key1: complete tree string without CF and branch lengths; value: number		In -> defined, OUT -> unchanged
		my	$href_cf_of_subtree_of_treefile		= $_[3]	; # key1: tree filename; key2: subtree string; value: subtree related cf score	In -> defined, OUT -> unchanged
		my	$href_cf_of_cpltree_of_treefile		= $_[4]	; # key1: tree filename; key2: raw nwk string; value: total cf score			In -> defined, OUT -> unchanged
		my	$sref_outpath						= $_[5]	; # path to new cf-tree summarized result prints								In -> defined, OUT -> unchanged
		my	$sref_cf_setting					= $_[6]	; # specified cf type (scf or gcf measure) impedes further filehandling 		In -> defined, OUT -> unchanged
		my	$sref_rplot_parameter				= $_[7]	; # if 1 -> R plot, otherwise no R plot											In -> defined, OUT -> unchanged
		
		######## list sampling of BL & subtree individual CF support #####################
		##################################################################################
		# Filtering of filename coded branch elongation value (BL3)
		# and BL3 individual cf sampling of treefile individual subtrees
		# in %hol_cf_values_of_subtree_of_bl
		my	%hol_cf_values_of_subtree_of_bl ; # key1: branch length BL3; key2: subtree newick; value: array of overall treefiles extracted subtree & BL corresponding cf support 
		my	%hol_cf_averag_of_cpltree_of_bl ; # key1: branch length BL3; key2: complete raw newick; value: list of tree individually averaged CF support
		
		for my $iqtree_file ( sort keys %$href_seen_iqtree ){
			
			#########################################
			# path substition in tree filename
			( my $iqtree_file_modified	=	$iqtree_file )	=~	s/^(.*\/)+// ;	#print "\n", $iqtree_file_modified ;
			#########################################
			
			#########################################
			# tree filename consists of original reftree name and alignment filename, e.g.
			# 12345_0.1_TRUE_1_iqtree_treefile.tre.12345_0.1_TRUE_1_msa.cf.tree
			# with pattern '.tre.' always as seperator
			# extraction of msa filename for further info usage, e.g.
			# 12345_0.1_TRUE_1_msa.cf.tree -> $nameprts[1]
				$iqtree_file_modified	=~	s/\.tre\./:::/ ; 
			my	@nameprts				=	split ":::", $iqtree_file_modified ;	#print "\n", $nameprts[1] ;
			#########################################
			
			#########################################
			# extraction of bl3 branch length info
			# branch length bl3 coded in alignment filename on different
			# underscore positions (depending on the indelible bl1 to bl3 simulation setting):
			# 1) RandomSeedNumber_BL1_BL2_BL3
			# 2) RandomSeedNumber_BL1_BL2_BL3
			# 3) RandomSeedNumber_BL1_BL3
			# 4) RandomSeedNumber_BL2_BL3
			# 5) RandomS$subtree$subtreeeedNumber_BL3
			# these pattern can be followed by a single number ('_1') or alphanumeric signes ('_TRUE' or '_msa')
			# '_msa' should always be part of the filename, whereas '_TRUE' and '_1' appear optional
			my	$bl_coded_filename	=	$nameprts[1] 							;	#print "\n", $bl_coded_filename	;	# e.g. 12346_1.5_TRUE_1_msa.cf.tree
				$bl_coded_filename	=~	s/_TRUE_1_|_1_/_/ 						;	#print "\n", $bl_coded_filename	;	# e.g. 12346_1.5_msa.cf.tree
			
			########################################################################################################################### MISSING SPLIT PATTERN FOR gCF measure
			########################################################################################################################### gCF has not been tested yet (28.10.21)
			my	@prts_bl ;	
			if	( $$sref_cf_setting =~ /^s$/ ){ @prts_bl = split "_msa", $bl_coded_filename }
			else{ print "\nSorry, currently only site concordance measures can be further analyzed in respect of different branch length (BL) simulations\n\n"; exit } 	
			###########################################################################################################################
			###########################################################################################################################
			
			my	@subprts_bl			=	split 	"_"		, $prts_bl[0] 			;
			my	$bl					=	pop				  @subprts_bl 			;	#print "\n", $bl 				;	# e.g. 1.5
			#########################################
			
			#########################################
			# For each treefile, sample cf support of each subtree
			# in respect of similar simulated branch lengths BL3
			for	my $subtree ( sort keys %$href_seen_subtree ){
				
				if ( $href_cf_of_subtree_of_treefile->{$iqtree_file}{$subtree} ){	#print "\n", $href_cf_of_subtree_of_treefile->{$iqtree_file}{$subtree} ;
					
					push @{$hol_cf_values_of_subtree_of_bl{$bl}{$subtree}}, $href_cf_of_subtree_of_treefile->{$iqtree_file}{$subtree}
				}
			}
			#########################################
			
			#########################################
			# For each treefile, sample averaged overall tree support of each complete tree
			# in respect of similar simulated branch lengths BL3
			for my $raw_nwk ( sort keys %$href_seen_nwk_raw ){
				
				if ( $href_cf_of_cpltree_of_treefile->{$iqtree_file}{$raw_nwk} ){
					
					push @{$hol_cf_averag_of_cpltree_of_bl{$bl}{$raw_nwk}}, $href_cf_of_cpltree_of_treefile->{$iqtree_file}{$raw_nwk}
				}
			}
			#########################################
		}
		##################################################################################
		##################################################################################
		
		
		####### BL & subtree-newick individual CF mean support & tsv table print #########
		##################################################################################
		# print .tsv table file of BL3 and subtree individually
		# identified mean-cf support for further R plotting
		# mean CF for different lengths BL3 are listed vertically
		my	 	$tsv_out = $$sref_outpath."/mean_subtree_cf_support.tsv" ;
		open 	OUTtsv, ">$tsv_out" || die "Cannot write tsv file ", $$sref_outpath/$tsv_out, ": $!" ;
		#################
		
		#################
		# print tsv header
		print	OUTtsv	"BL\tSubtree\tCFmean\n" ;
		#################
		
		#################
		# hash sampling of BL3 and subtree individually identified mean-cf support
		my	%cf_mean_of_subtree_of_bl	; # key1: branch length BL3; key2: subtree newick; value: mean subtree & BL corresponding cf support 
		#################
		
		#################
		# mean cf calculation from all subtree and BL3 individually related cf scores
		for my 	$bl ( sort keys %hol_cf_values_of_subtree_of_bl ){
			
			for	my $subtree ( sort keys %$href_seen_subtree ){
				
				if ( $hol_cf_values_of_subtree_of_bl{$bl}{$subtree} ){
					
					#################
					# collect BL3 and subtree corresponding cf values
					my	@cf_scorings = @{$hol_cf_values_of_subtree_of_bl{$bl}{$subtree}} ;
					#################
					
					#################
					# calculate mean cf support
						my	$N_values	= @cf_scorings ;
						my	$sum_cf		= 0 ;
					for my	$value ( @cf_scorings ){ $sum_cf += $value }
						my	$mean_cf	= sprintf("%.1f", ($sum_cf / $N_values) ) ;	
					
					$cf_mean_of_subtree_of_bl{$bl}{$subtree} = $mean_cf ;
					#################
					
					#################
					# terminal test print
					#my $string_values = join "_", @cf_scorings; print "\nvalues: ",$string_values  ,"\nmean ", $mean_cf ;
					#################
					
					#################
					# tsv table print
					print OUTtsv $bl, "\t", $subtree, "\t", $mean_cf, "\n" ;
					#################
				}
				
				else{print OUTtsv $bl, "\t", $subtree, "\tNA\n" ;}
				#################
			}
			#################
		}
		#################
		
		#################
		# close tsv table print 
		close	OUTtsv ;
		##################################################################################
		##################################################################################
		
		
		####### BL & complete-newick individual CF mean support & tsv table print ########
		##################################################################################
		# print .tsv table file of BL3 and complete-tree individually
		# identified mean-cf support for further R plotting
		# mean CF for different lengths BL3 are listed vertically
		my	 	$tcv_out = $$sref_outpath."/mean_complTree_cf_support.tsv" ;
		open 	OUTtcv, ">$tcv_out" || die "Cannot write tsv file ", $$sref_outpath/$tcv_out, ": $!" ;
		#################
		
		#################
		# print tsv header
		print	OUTtcv	"BL\tComplTree\tCFavrg_cmpl\n" ;
		#################
		
		#################
		# hash sampling of BL3 and tree individually identified mean-cf support
		my	%cf_mean_of_cpltree_of_bl	; # key1: branch length BL3; key2: complete-tree newick; value: mean tree & BL corresponding cf support 
		#################
		
		#################
		# mean cf calculation from all trees and BL3 individually related cf scores
		for my 	$bl ( sort keys %hol_cf_averag_of_cpltree_of_bl ){
			
			for	my $cpltree ( sort keys %$href_seen_nwk_raw ){
					
				if ( $hol_cf_averag_of_cpltree_of_bl{$bl}{$cpltree} ){
					
					#################
					# collect BL3 and complete-tree corresponding cf values
					my	@cf_scorings = @{$hol_cf_averag_of_cpltree_of_bl{$bl}{$cpltree}} ;
					#################
					
					#################
					# calculate mean cf support
						my	$N_values	= @cf_scorings ;
						my	$sum_cf		= 0 ;
					for my	$value ( @cf_scorings ){ $sum_cf += $value }
						my	$mean_cf	= sprintf("%.1f", ($sum_cf / $N_values) ) ;	
				
					$cf_mean_of_cpltree_of_bl{$bl}{$cpltree} = $mean_cf ;
					#################
					
					#################
					# terminal test print
					#my $string_values = join "_", @cf_scorings; print "\nvalues: ",$string_values  ,"\nmean ", $mean_cf ;
					#################
					
					#################
					# tsv table print
					print OUTtcv $bl, "\t", $cpltree, "\t", $mean_cf, "\n" ;
					#################
				} ;
				#################
			}
			#################
		}
		#################
		
		#################
		# close tsv table print 
		close	OUTtcv ;
		##################################################################################
		##################################################################################
		
		
		####### BL & subtree individual CF mean txt table print ##########################
		##################################################################################
		# print .txt table file of BL3 and subtree individually
		# identified mean-cf support for further processing
		# mean CF for different lengths BL3 are listed horizontally
		my	 	$txt_out = $$sref_outpath."/mean_subtree_cf_support.txt" ;
		open 	OUTtxt, ">$txt_out" || die "Cannot write txt file ", $$sref_outpath/$txt_out, ": $!" ;
		#################
		
		#################
		# print table header
		print	OUTtxt	"Subtree" ;
		for my	$length_bl ( sort {$a<=>$b} keys %cf_mean_of_subtree_of_bl ){ print OUTtxt "\tBL: ", $length_bl }
		print	OUTtxt	"\n" ;
		#########################################
		
		#################
		# data print
		for	my $subtree ( sort keys %$href_seen_subtree ){
			
			print OUTtxt $subtree ;
			
			for my $bl ( sort keys %cf_mean_of_subtree_of_bl ){
				
				if  ( $cf_mean_of_subtree_of_bl{$bl}{$subtree} ){ print OUTtxt "\t", $cf_mean_of_subtree_of_bl{$bl}{$subtree} }
				else{ print OUTtxt "\tNA" }
			}
			
			print OUTtxt "\n"
		}
		#################
		
		#################
		# close txt table print 
		close	OUTtxt ;
		##################################################################################
		##################################################################################
		
		
		####### Generate Latex table of Mean CF Support for each length bl2 and subtree ##
		##################################################################################
		my	@table_lines2 ;
		
		#################
		# do string of subtree header line 
		my		@header_elements2	=	sort {$a<=>$b} keys %cf_mean_of_subtree_of_bl ;
		unshift @header_elements2	,	"Subtree" ;
		my		$Ncols2			 	=	@header_elements2 ;
		my		$header_string2		=	join "\t", @header_elements2 ;	#print "\n", $header_string ; exit ;
		push	@table_lines2	 	,	$header_string2 ;
		#################
		
		#################
		# do subtree data lines
		for	my $subtree ( sort keys %$href_seen_subtree ){
			
			my	$dataline = $subtree ;
			
			for my $bl ( sort keys %cf_mean_of_subtree_of_bl ){
				
				if  ( $cf_mean_of_subtree_of_bl{$bl}{$subtree} ){ $dataline .= "\t".$cf_mean_of_subtree_of_bl{$bl}{$subtree} }
				else{ $dataline .= "\tNA" }
			}
			
			push @table_lines2, $dataline
		}
		#################
		
		#################
		# define caption text
		my	$table_caption2 = "Average CF support for subtree relationships given different steps of branch elongation" ;
		#################
		
		#################
		# full output filename (with output-path) 
		my	$tex_filename2 = $$sref_outpath."/mean_subtree_cf_support_table.tex" ;
		#################
		
		#################
		# do latex table
		&latex_table_simple(
			
			\@table_lines2		,	# list of correctly ordered table lines (with header line as first element and column values tab delimited in each line)	In -> defined, OUT -> unchanged
			\$Ncols2			,	# N columns																													In -> defined, OUT -> unchanged
			\$tex_filename2		,	# filename																													In -> defined, OUT -> unchanged
			\$table_caption2	,	# caption text																												In -> defined, OUT -> unchanged
			\$$sref_outpath		,	# Output path																												In -> defined, OUT -> unchanged
			\"c"				,	# table placement (l: left, c: center, r: right)																			In -> defined, OUT -> unchanged
		) ;
		##################################################################################
		##################################################################################
		
		
		## Generate Latex table of Mean CF Support for each length bl2 and complete tree #
		##################################################################################
		my	@table_lines1 ;
		
		#################
		# do string of complete-tree header line 
		my		@header_elements1	=	sort {$a<=>$b} keys %cf_mean_of_cpltree_of_bl ;
		unshift @header_elements1	,	"Tree" ;
		my		$Ncols1			 	=	@header_elements1 ;
		my		$header_string1	 	=	join "\t", @header_elements1 ;	#print "\n", $header_string1 ; exit ;
		push	@table_lines1	 	,	$header_string1 ;
		#################
		
		#################
		# do complete-tree data lines
		for	my $cpltree ( sort keys %$href_seen_nwk_raw ){
			
			my	$dataline = $cpltree ;
			
			for my $bl ( sort keys %cf_mean_of_cpltree_of_bl ){
				
				if  ( $cf_mean_of_cpltree_of_bl{$bl}{$cpltree} ){ $dataline .= "\t".$cf_mean_of_cpltree_of_bl{$bl}{$cpltree} }
				else{ $dataline .= "\tNA" }
			}
			
			push @table_lines1, $dataline
		}
		#################
		
		#################
		# define caption text
		my	$table_caption1 = "Average CF support for tree relationships given different steps of branch elongation" ;
		#################
		
		#################
		# full output filename (with output-path) 
		my	$tex_filename1 = $$sref_outpath."/mean_complTree_cf_support_table.tex" ;
		#################
		
		#################
		# do latex table
		&latex_table_simple(
			
			\@table_lines1		,	# list of correctly ordered table lines (with header line as first element and column values tab delimited in each line)	In -> defined, OUT -> unchanged
			\$Ncols1			,	# N columns																													In -> defined, OUT -> unchanged
			\$tex_filename1		,	# filename																													In -> defined, OUT -> unchanged
			\$table_caption1	,	# caption text																												In -> defined, OUT -> unchanged
			\$$sref_outpath		,	# Output path																												In -> defined, OUT -> unchanged
			\"l"				,	# table placement (l: left, c: center, r: right), c if undefined															In -> defined, OUT -> unchanged
		) ;
		##################################################################################
		##################################################################################
	}
	#####
	
	#####
	sub latex_table_simple{
		
		############################################################################# 
		# called from &bl2_subtree_analysis
		#
		# prints a latex .tex table file from a list of string coded data lines with
		# line columns tabstop delimited @$aref_table_lines
		# latex table is further converted to pdf using pdflatex
		############################################################################# 
		
		my	$aref_table_lines	= $_[0] ;	# list of correctly ordered table lines (with header line as first element and column values tab delimited in each line)	In -> defined, OUT -> unchanged
		my	$sref_Ncols			= $_[1] ;	# N columns																													In -> defined, OUT -> unchanged
		my	$sref_tex_filename	= $_[2] ;	# filename																													In -> defined, OUT -> unchanged
		my	$sref_caption		= $_[3] ;	# caption text																												In -> defined, OUT -> unchanged
		my	$sref_print_path	= $_[4] ;	# output path																												In -> defined, OUT -> unchanged
		my	$sref_tab_position	= $_[5] ;	# table placement (l: left, c: center, r: right), c if undefined															In -> defined, OUT -> unchanged
		
		###################################################
		# Header preparation
		my	$header_line		=	shift @$aref_table_lines ;
		my	@header_cols		=	split "\t", $header_line ;
		for my $hcol ( @header_cols ){ 
			
			$hcol =~ s/_/\_/g 					; # protect underscores in latex (substitute '_' to '\_') 
			$hcol =  "\\textbf{".$hcol."}"		; # textbf
		}
		$header_line = join " & ", @header_cols ;
		###################################################
		
		###################################################
		# data line preparation
		my	$N_lines = @$aref_table_lines ;
		for	my $i ( 0 .. $N_lines-1 ){
			
			my @lineprts = split "\t", $aref_table_lines->[$i] ;
			
			for my $lcol ( @lineprts ){ 
			
				$lcol =~ s/_/\_/g 							; # protect underscores in latex (substitute '_' to '\_') 
				$lcol =  "\\textcolor{black}{".$lcol."}"	; # allows subsequent changes of textcolor
			}
			
			$aref_table_lines->[$i] = join " & ", @lineprts 	;
			
			if ( $i % 2 == 0 ){ $aref_table_lines->[$i] = "\\rowcolor{cell} ".$aref_table_lines->[$i] }
		}
		my	$datalines = join "\\\\\n", @$aref_table_lines ;
			$datalines = $datalines."\\\\\n" ;
		###################################################
		
		###################################################
		# cmidrule preparation
		my	$cmidrule_string ;
		for my $col ( 1 .. $$sref_Ncols ){ $cmidrule_string .= "\\cmidrule(lr){".$col."-".$col."}" }
		###################################################
		
		###################################################
		# open tex file
		open ( my $FH, '>', $$sref_tex_filename ) or die "Cannot open $$sref_tex_filename : $!\n" ;
		###################################################
		
		###################################################
		# Define table position
		my	$table_position = "c" ;	# default is center position
		if ( ( $$sref_tab_position ) && ( $$sref_tab_position =~ /^l$|^r$/ ) ){ $table_position = $$sref_tab_position }
		###################################################
		
		###################################################
		print {$FH}	"\\documentclass[10pt,oneside,a4paper]{article}\n"							,
					"\n"																		,
					"\\usepackage{longtable}\n"													,
					"\\usepackage[table]{xcolor}\n"												,
					"\\usepackage{helvet}\n"													,
					"\\usepackage{booktabs}\n"													,
					"\\usepackage[width=1.25\\textwidth]{caption}\n"							,
					"\n"																		,
					"\\renewcommand{\\familydefault}{\\sfdefault}\n"							,
					"\n"																		,
					"\n"																		,
					"\\definecolor{cell}{RGB}{220,230,240}\n"									,
					"\\definecolor{line}{RGB}{220,20,60}\n"										,
					"\n"																		,
					"\n"																		,
					"\\begin{document}\n"														,
					"\n"																		,
					"\n"																		,
					"\\arrayrulecolor{line}\n"													,
					"\n"																		,
					"\\begin{longtable}[",$table_position,"]{", "l" x $$sref_Ncols, "}\n"		,
					"\n"																		,
					"\\caption[VARIABEL]{\\label{tab:parameter}",$$sref_caption,"}\\\\\n"		,
					"\\hline \\toprule\n"														,
					"\n"																		,
					$header_line																,
					"\\\\\n"																	, 
					$cmidrule_string															,
					"\n"																		,
					'\\endfirsthead'															,
					"\n"																		,
					"\n"																		,
					"\\multicolumn{".$$sref_Ncols."}{r}{Continues on next page(s).}\\\\"		,
					"\n"																		,
					'\\endfoot'																	,
					"\n"																		,
					'\\hline'																	,
					"\n"																		,
					"\n"																		,
					"\\multicolumn{".$$sref_Ncols."}{r}{End of table \\ref{tab:parameter}}\\\\"	,
					"\n"																		,
					'\\endlastfoot'																,
					"\n"																		,
					"\n"																		,
					$datalines																	,
					"\n"                                                                        ,
					"\\bottomrule"																,
					"\n"																		,
					"\\end{longtable}"															,
					"\n"																		,
					"\\end{document}"															;
		
		close $FH	;
		
		####################################
		# must be executed twice, otherwise table format can be incorrect
		system ("pdflatex  $$sref_tex_filename") ;
		system ("pdflatex  $$sref_tex_filename") ;
		####################################
		
		#################################### latex outfile handling
		# move latex pdf from main appetite folder
		# to specified output path
		( my $latexpdf_new = $$sref_tex_filename ) =~ s/\.tex$/\.pdf/ ;
		( my $latexpdf_old = $latexpdf_new		 ) =~ s/^$$sref_print_path\/// ;
		( my $latexout_inf = $latexpdf_old		 ) =~ s/\.pdf$// ;
		move ( $latexpdf_old, $latexpdf_new ) ;
		##
		
		##
		# delete other latex outfiles (.aux, .log) in appetite main folder
		for my $latexinfofile ( <$latexout_inf.*> ){ unlink $latexinfofile }
		####################################
	}
	#####
	
	#####
	sub structure_handling{
		
		############################################################################# 
		# called from &summarize_cf_trees
		#
		# newick string decomposition into all internal node classified
		# subtree relationships and assignment of subtree related cf score
		# @substrings stores final substrings and corresponding cf score as separate elements, 
		# e.g. qw/(A1,A2), 99.5, (B1,B2), 99.6, .../  
		# 
		# return @substrings to &summarize_cf_trees
		############################################################################# 
		
		#########################################
		# subroutine can handle only iqtree formatted newick strings, like:
		# ((((A1,A2)99.5,(B1,B2)99.6)41.8,((C1,C2)99.5,(D1,D2)99.5)41.9)99.6,O1,O2);
		my	$sref_string        = $_[0] ;	#print "\ntree ", $$sref_string ;
		#########################################
		
		#########################################
		# subroutine global variables  
		my ( 
			@substrings , # stores final substrings and corresponding cf score as separate elements, e.g. qw/(A1,A2), 99.5, (B1,B2), 99.6, .../ 
			@forward	, # stores site number of unmatching open brakets during the analysis
		);
		#########################################
		
		#########################################
		# split newick string at each position and store each element in @structures
		my @structures = split "", $$sref_string ;
		#########################################
		
		#########################################
		# extraction of internal node relationships and corresponding cf support, e.g.:
		# ((((A1,A2)99.5,(B1,B2)99.6)41.8,((C1,C2)99.5,(D1,D2)99.5)41.9)99.6,O1,O2);
		# 
		# ...separated in substrings...
		# (A1,A2)99.5
		# (B1,B2)99.6
		# ((A1,A2),(B1,B2))41.9
		# (C1,C2)99.5
		# (D1,D2)99.5
		# ((C1,C2),(D1,D2))41.9
		# (((A1,A2),(B1,B2)),((C1,C2),(D1,D2)))99.6
		#
		# ...with each substring stored in @substrings
		my $flag = 0 ;
		
		for my $pos ( 0 .. @structures-1 ){
			
			#########################################
			# if a closed bracket ')' has been identified earlier as substring-end ($flag==1), extraction and cleaning of substring :
			# i)  the position of the last open braket '(' seen in front of the closed braket marks the substring start position ($substr_start),
			# ii) the closed braket following commata ',' or open '(' or closed ')' bracket defines the end of a substring, e.g.:
			# substring '(A1,A2)99.5,' in  ((((A1,A2)99.5,(B1,B2)99.6)41.8,((C1,C2)99.5,(D1,D2)99.5)41.9)99.6,O1,O2);
			if ((	$structures[$pos]	=~ /,|\(|\)/ 	) && 
				(	$flag == 1 							) ){
				
				####################
				# substring extraction ($substring)
				my		$substr_start 	= 	pop @forward			; # substring start position (last open bracket position before)
				my		$substr_length	= 	$pos - $substr_start	; # substring length identified by start and end position (end position defined by $pos)
				my		$substring		= 	substr( $$sref_string, $substr_start, $substr_length ) ; # substring extraction
				####################
				
				####################
				# substring cleaning of further non-directly substring related cf support of more inner nodes, e.g.:
				# ((A1,A2)99.5,(B1,B2)99.6)41.9 to ((A1,A2),(B1,B2))41.9
				# in some substring cases while loop embedded pattern substitution must be done multiple times
						$substring												=~	s/\)\d+(\.\d+)?,/\),/g	 ; # first step:   ((A1,A2)99.5,(B1,B2)99.6)41.9 to ((A1,A2),(B1,B2)99.6)41.9
				while ( $substring		=~	 /\)\d+(\.\d+)?\)/	){ $substring	=~	s/\)\d+(\.\d+)?\)/\)\)/g } # second step:  ((A1,A2),(B1,B2)99.6)41.9 to ((A1,A2),(B1,B2))41.9
				####################
				
				####################
				# substring splitting, separating cf score from original substring
						$substring		=~	s/\)(\d+)/\)::$1/ 		;	# add '::' betwenn substring and cf score
				my		@subparts		=	split "::",	$substring	;	# use '::' as string separator
				#for my $part ( @subparts ){ print "\n", $part }
				####################
				
				####################
				# substring sorage in @substrings
				push	@substrings		, 	@subparts		; 
				####################
				
				####################
				# flag resetting		
				$flag = 0 ;
				####################
			}
			#########################################
			
			#########################################
			# push position of open braket as start position in @forward
			if ( 	$structures[$pos]	=~ /\(/   ){ push @forward, $pos }#; print "\nstart\t", $structures[$pos]	}
			#########################################
			
			#########################################
			# end of a substring, substring corresponding cf value follows directly afterwards
			# thus, set flag of substring end to 1
			if ( 	$structures[$pos]	=~ /\)/   ){ $flag = 1			 }#; print "\nflag!\t", $structures[$pos]	}
			#########################################
		}
		#########################################
		
		####################
		# terminal test print
		# for ( @substrings ){ print "\n", $_ } print "\n" ; #exit ;
		####################
		
		####################
		# return @substrings
		return @substrings ;
		####################
	}
	#####

	#####
	sub rename_cf_output{
		
		############################################################################# 
		# called from &concordance (main routine)
		#
		# substitutes old filenames to new filenames in specified filepath
		############################################################################# 
		
		my	$sref_filepath		=	$_[0]	;	# path to located cf outputfiles, usually the reftree input path	In -> defined, OUT -> unchanged
		my	$sref_old_string	=	$_[1]	;	# string to replace in cf output files								In -> defined, OUT -> unchanged
		my	$sref_new_string	=	$_[2]	;	# string to add to cf output files									In -> defined, OUT -> unchanged
		
		for my $file ( <$$sref_filepath/*> ){
			
			if ( $file	=~	/${$sref_old_string}/ ){
				
				#print "\n", $file ;
				
				( my $new_filename	= $file ) =~ s/${$sref_old_string}/${$sref_new_string}/ ;
				move ( $file, $new_filename ) ;
			}
		}
	}
	#############################################################################
}
############################################################################################################################
############################################################################################################################



























