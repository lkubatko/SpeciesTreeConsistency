#!/usr/bin/perl

print "\n\n\n Simulation starting ..... \n\n\n";


printf "Running tree with branch length setting 1 \n\n\n";

    foreach my $ngenes (2000, 3000, 4000, 5000, 7500, 10000) {

	# ngenes is a parameter for the number of genes to simulate
	# the theory is derived for the case where each site has its own gene tree; to simulate this, pick a large number of genes,
	# and give each length 1 (genelength parameter below) 

	system("/bin/rm -f alldat.out alldat_use.out ML.out ML_details.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps");
	
	# the entire script can be run multiple times by changing the value of NLOOP
	use constant NLOOP1 => 500;  
	$count = 0;
	while ( $count < NLOOP1)
	{
	    
	    ###########################                                                                                            
   	    ### Simulate data       ###                                                                                         
	    ###########################        
	    
	    # number of sites in each gene
	    $genelength = 1;
	    
	    # the code below writes the input file for the "intra" program, which will simulate the gene trees from the coalescent model
	    # the only thing that needs to be changed is the name of the file containing the species tree:  mysptree4tax1 
	    open MYINFILE, ">ourinfile";
	    select MYINFILE;
	    print "begin coal;\n";
	    print "  ntaxa = 4;\n";
	    print "  species tree file = mysptree4tax_bl1;\n";
	    print "  intra no;\n";
	    print "  theta = 2;\n";
	    print " \n"; 
	    print "  gene tree file = simulated;\n";
	    print " \n"; 
	    print "  taxa names = Species1 Species2 Species3 Species4;\n";
	    print " \n";
	    print "  blstyle none;\n"; 
	    print " \n";
	    print "  outfile = myout;\n";
	    print " \n";
	    print "  ngtrees ", $ngenes, ";\n";   
	    print "  nstrees 1;\n";
	    my $range = 100000;
	    my $random_number = int(rand($range));
	    print "  seed1 = ", $random_number, "\n"; 
	    my $random_number = int(rand($range));
	    print "  seed2 = ", $random_number, "\n";
	    print " \n";
	    print "end;\n";
	    close MYINFILE;
	    
	    $nsites = $ngenes*$genelength;
	    system("cp ourinfile infile");
	    # the command below runs the intra program
	    system("./intra > outfile");
	    
	    # get the file with the gene trees ready for seq-gen:
	    system(qq( awk '{print "[$genelength]" \$1}' simtrees.dat > simtrees2.dat));
	    # The commands below call seq-gen to simulate sites along the gene trees; theta is set using the -s flag, e.g., the code below has theta=0.01
	    
	    # JC69 model
	    system("./seq-gen < simtrees2.dat > infile -q -mHKY -s 0.01 -l $nsites -p $ngenes");
	    system(qq( awk 'NR<2{print \$0; next}{print | "sort"}' infile > data.phy ));

	    ###########################
	    ### Get SVDQ estimate   ### 
	    ###########################
	    
	    system("/bin/cp -f settings_PrintFlat settings");
	    system("./PrintFlat");
	    system("cat svdq.out species_quartets.out > temp.svdq.out");
	    system("/bin/cp -f temp.svdq.out svdq.out");
	    
	    # check that there is enough data to compute scores
	    system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile));
	    system("cat temp.scoresfile valid.svdq > temp.valid.svdq");
	    system("/bin/cp -f temp.valid.svdq valid.svdq");
	    
	    #################################                                                                                                                       
            ### Get SVDQ estimate - SNPs  ###                                                                                                                
            #################################                                                                                                                       

            system("./PrintFlatSNPs");                                                                                                
            system("cat svdq.snps.out species_quartets.out > temp.svdq.snps.out");
            system("/bin/cp -f temp.svdq.snps.out svdq.snps.out");

            # check that there is enough data to compute scores                                                      
            system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile.snps));            
	    system("cat temp.scoresfile.snps valid.svdq.snps > temp.valid.svdq.snps");                                                                             
            system("/bin/cp -f temp.valid.svdq.snps valid.svdq.snps");   

	    ###########################                                                                                                            
	    ### Get the MLE         ###                                                                                                                    
	    ###########################        
       
	    # Now run ssa to get the site pattern count for each set of 4 taxa of interest
	    # (1,2,3,4)
	    system("/bin/cp -f settings_ssa settings");
	    system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("/bin/cp -f infile_data check1");
	    system("./ssa > temp1.ssa");
       
	    # (1,3,2,4)
	    system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp2.ssa");
       
	    # (1,4,2,3)
	    system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp3.ssa");
       
	    # (2,1,3,4)
	    system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp4.ssa"); 
       
	    # (2,3,1,4)                                                                                                            
	    system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp5.ssa");  
       
	    # (2,4,1,3)                                                                                                            
	    system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp6.ssa");     
       
	    # (3,1,2,4)                                                                                                            
	    system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp7.ssa");    
       
	    # (3,2,1,4)                                                                                                            
	    system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data)); 
	    system("./ssa > temp8.ssa");                    
       
	    # (3,4,1,2)                                                                                                            
	    system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp9.ssa");      
       
	    # (4,1,2,3)                                                                                                            
	    system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));       
	    system("./ssa > temp10.ssa"); 
       
	    # (4,2,1,3)                                                                                                            
	    system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp11.ssa"); 
       
	    # (4,3,1,2)                                                                                                            
	    system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
	    system("./ssa > temp12.ssa");

	    # Use R code to find the ML tree
	    system("R CMD BATCH ./get_lik.R");
	    system("cat temp.R ML.out > temp.ML");
	    system("/bin/cp -f temp.ML ML.out");

	    $count = $count+1;
	    
	}
	
	system("paste ML.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps > alldat.out");
	system("sed -f sed.tree.ind alldat.out > alldat_bl1_$ngenes.out");	
    }

    system("zip alldat_bl1_0.01.zip alldat_bl1_*");


printf "Running tree with branch length setting 2 \n\n\n";

    foreach my $ngenes (2000, 3000, 4000, 5000, 7500, 10000) {   

	# ngenes is a parameter for the number of genes to simulate
	# the theory is derived for the case where each site has its own gene tree; to simulate this, pick a large number of genes,
	# and give each length 1 (genelength parameter below) 
	# $ngenes = 200; 
	
	system("/bin/rm -f alldat.out alldat_use.out ML.out ML_details.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps");
	
	# the entire script can be run multiple times by changing the value of NLOOP
	use constant NLOOP1 => 500;  
	$count = 0;
	while ( $count < NLOOP1)
	{
	    
	    ###########################                                                                                                                        
	    ### Simulate data       ###                                                                                                                       
	    ###########################        
	    
	    # number of sites in each gene
	    $genelength = 1;
	    
	    # the code below writes the input file for the "intra" program, which will simulate the gene trees from the coalescent model
	    # the only thing that needs to be changed is the name of the file containing the species tree:  mysptree4tax1 
	    open MYINFILE, ">ourinfile";
	    select MYINFILE;
	    print "begin coal;\n";
	    print "  ntaxa = 4;\n";
	    print "  species tree file = mysptree4tax_bl2;\n";
	    print "  intra no;\n";
	    print "  theta = 2;\n";
	    print " \n"; 
	    print "  gene tree file = simulated;\n";
	    print " \n"; 
	    print "  taxa names = Species1 Species2 Species3 Species4;\n";
	    print " \n";
	    print "  blstyle none;\n"; 
	    print " \n";
	    print "  outfile = myout;\n";
	    print " \n";
	    print "  ngtrees ", $ngenes, ";\n";   
	    print "  nstrees 1;\n";
	    my $range = 100000;
	    my $random_number = int(rand($range));
	    print "  seed1 = ", $random_number, "\n"; 
	    my $random_number = int(rand($range));
	    print "  seed2 = ", $random_number, "\n";
	    print " \n";
	    print "end;\n";
	    close MYINFILE;
	    
	    $nsites = $ngenes*$genelength;
	    system("cp ourinfile infile");
	    # the command below runs the intra program
	    system("./intra > outfile");
	    
	    # get the file with the gene trees ready for seq-gen:
	    system(qq( awk '{print "[$genelength]" \$1}' simtrees.dat > simtrees2.dat));
	    # The commands below call seq-gen to simulate sites along the gene trees; theta is set using the -s flag, e.g., the code below has theta=0.01
	    
	    # JC69 model
	    system("./seq-gen < simtrees2.dat > infile -q -mHKY -s 0.01 -l $nsites -p $ngenes");
	    system(qq( awk 'NR<2{print \$0; next}{print | "sort"}' infile > data.phy ));

	    ###########################
	    ### Get SVDQ estimate   ### 
	    ###########################
	    
	    system("/bin/cp -f settings_PrintFlat settings");
	    system("./PrintFlat");
	    system("cat svdq.out species_quartets.out > temp.svdq.out");
	    system("/bin/cp -f temp.svdq.out svdq.out");
	    
	    # check that there is enough data to compute scores
	    system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile));
	    system("cat temp.scoresfile valid.svdq > temp.valid.svdq");
	    system("/bin/cp -f temp.valid.svdq valid.svdq");
	    
	    #################################                                                         
	    ### Get SVDQ estimate - SNPS  ###
	    #################################

            system("./PrintFlatSNPs");                                                                                                                        
            system("cat svdq.snps.out species_quartets.out > temp.svdq.snps.out");                                                                  
            system("/bin/cp -f temp.svdq.snps.out svdq.snps.out");                                                                                                   

            # check that there is enough data to compute scores                                                                                                      
            system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile.snps));                                                                        
            system("cat temp.scoresfile.snps valid.svdq.snps > temp.valid.svdq.snps");                                                                               
            system("/bin/cp -f temp.valid.svdq.snps valid.svdq.snps");     

       ###########################                                                                                                            
       ### Get the MLE         ###                                                                                                                    
       ###########################        
       
       # Now run ssa to get the site pattern count for each set of 4 taxa of interest
       # (1,2,3,4)
       system("/bin/cp -f settings_ssa settings");
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("/bin/cp -f infile_data check1");
       system("./ssa > temp1.ssa");
       
       # (1,3,2,4)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp2.ssa");
       
       # (1,4,2,3)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp3.ssa");
       
       # (2,1,3,4)
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp4.ssa"); 
       
       # (2,3,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp5.ssa");  
       
       # (2,4,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp6.ssa");     
       
       # (3,1,2,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
        system("./ssa > temp7.ssa");    
       
       # (3,2,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data)); 
       system("./ssa > temp8.ssa");                    
       
       # (3,4,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp9.ssa");      
       
       # (4,1,2,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));       
       system("./ssa > temp10.ssa"); 
       
       # (4,2,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp11.ssa"); 
       
       # (4,3,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp12.ssa");

       # Use R code to find the ML tree
       system("R CMD BATCH ./get_lik.R");
       system("cat temp.R ML.out > temp.ML");
       system("/bin/cp -f temp.ML ML.out");

	    $count = $count+1;
	    
	}
	
	system("paste ML.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps > alldat.out");
	system("sed -f sed.tree.ind alldat.out > alldat_bl2_$ngenes.out");                                                                                      
    }                                                                                                                                                          
                                                                                                                                         
    system("zip alldat_bl2_0.01.zip alldat_bl2_*");    



printf "Running tree with branch length setting 3 \n\n\n";

    foreach my $ngenes (2000, 3000, 4000, 5000, 7500, 10000) {   

	# ngenes is a parameter for the number of genes to simulate
	# the theory is derived for the case where each site has its own gene tree; to simulate this, pick a large number of genes,
	# and give each length 1 (genelength parameter below) 
	# $ngenes = 200; 
	
	system("/bin/rm -f alldat.out alldat_use.out ML.out ML_details.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps");
	
	# the entire script can be run multiple times by changing the value of NLOOP
	use constant NLOOP1 => 500;  
	$count = 0;
	while ( $count < NLOOP1)
	{
	    
	    ###########################                                                                                                                        
	    ### Simulate data       ###                                                                                                                       
	    ###########################        
	    
	    # number of sites in each gene
	    $genelength = 1;
	    
	    # the code below writes the input file for the "intra" program, which will simulate the gene trees from the coalescent model
	    # the only thing that needs to be changed is the name of the file containing the species tree:  mysptree4tax1 
	    open MYINFILE, ">ourinfile";
	    select MYINFILE;
	    print "begin coal;\n";
	    print "  ntaxa = 4;\n";
	    print "  species tree file = mysptree4tax_bl3;\n";
	    print "  intra no;\n";
	    print "  theta = 2;\n";
	    print " \n"; 
	    print "  gene tree file = simulated;\n";
	    print " \n"; 
	    print "  taxa names = Species1 Species2 Species3 Species4;\n";
	    print " \n";
	    print "  blstyle none;\n"; 
	    print " \n";
	    print "  outfile = myout;\n";
	    print " \n";
	    print "  ngtrees ", $ngenes, ";\n";   
	    print "  nstrees 1;\n";
	    my $range = 100000;
	    my $random_number = int(rand($range));
	    print "  seed1 = ", $random_number, "\n"; 
	    my $random_number = int(rand($range));
	    print "  seed2 = ", $random_number, "\n";
	    print " \n";
	    print "end;\n";
	    close MYINFILE;
	    
	    $nsites = $ngenes*$genelength;
	    system("cp ourinfile infile");
	    # the command below runs the intra program
	    system("./intra > outfile");
	    
	    # get the file with the gene trees ready for seq-gen:
	    system(qq( awk '{print "[$genelength]" \$1}' simtrees.dat > simtrees2.dat));
	    # The commands below call seq-gen to simulate sites along the gene trees; theta is set using the -s flag, e.g., the code below has theta=0.01
	    
	    # JC69 model
	    system("./seq-gen < simtrees2.dat > infile -q -mHKY -s 0.01 -l $nsites -p $ngenes");
	    system(qq( awk 'NR<2{print \$0; next}{print | "sort"}' infile > data.phy ));

	    ###########################
	    ### Get SVDQ estimate   ### 
	    ###########################
	    
	    system("/bin/cp -f settings_PrintFlat settings");
	    system("./PrintFlat");
	    system("cat svdq.out species_quartets.out > temp.svdq.out");
	    system("/bin/cp -f temp.svdq.out svdq.out");
	    
	    # check that there is enough data to compute scores
	    system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile));
	    system("cat temp.scoresfile valid.svdq > temp.valid.svdq");
	    system("/bin/cp -f temp.valid.svdq valid.svdq");
	    
	    #################################                                                                                                                                    ### Get SVDQ estimate - SNPs  ###                                                                                                                                    #################################                                                                                                                         
            system("./PrintFlatSNPs");                                                                                                                                           system("cat svdq.snps.out species_quartets.out > temp.svdq.snps.out");                                                                                               system("/bin/cp -f temp.svdq.snps.out svdq.snps.out");                                                                                                   

            # check that there is enough data to compute scores                                                                                                      

            system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile.snps));                                                                                    system("cat temp.scoresfile.snps valid.svdq.snps > temp.valid.svdq.snps");                                                                                           system("/bin/cp -f temp.valid.svdq.snps valid.svdq.snps");     

	           ###########################                                                                                                            
       ### Get the MLE         ###                                                                                                                    
       ###########################        
       
       # Now run ssa to get the site pattern count for each set of 4 taxa of interest
       # (1,2,3,4)
       system("/bin/cp -f settings_ssa settings");
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("/bin/cp -f infile_data check1");
       system("./ssa > temp1.ssa");
       
       # (1,3,2,4)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp2.ssa");
       
       # (1,4,2,3)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp3.ssa");
       
       # (2,1,3,4)
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp4.ssa"); 
       
       # (2,3,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp5.ssa");  
       
       # (2,4,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp6.ssa");     
       
       # (3,1,2,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
        system("./ssa > temp7.ssa");    
       
       # (3,2,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data)); 
       system("./ssa > temp8.ssa");                    
       
       # (3,4,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp9.ssa");      
       
       # (4,1,2,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));       
       system("./ssa > temp10.ssa"); 
       
       # (4,2,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp11.ssa"); 
       
       # (4,3,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp12.ssa");

       # Use R code to find the ML tree
       system("R CMD BATCH ./get_lik.R");
       system("cat temp.R ML.out > temp.ML");
       system("/bin/cp -f temp.ML ML.out");
       
	    $count = $count+1;
	    
	}
	
	system("paste ML.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps > alldat.out");
	system("sed -f sed.tree.ind alldat.out > alldat_bl3_$ngenes.out");  
	
    }

    system("zip alldat_bl3_0.01.zip alldat_bl3_*");  



printf "Running tree with branch length setting 4 \n\n\n";

    foreach my $ngenes (2000, 3000, 4000, 5000, 7500, 10000) {   

	# ngenes is a parameter for the number of genes to simulate
	# the theory is derived for the case where each site has its own gene tree; to simulate this, pick a large number of genes,
	# and give each length 1 (genelength parameter below) 
	# $ngenes = 200; 
	
	system("/bin/rm -f alldat.out alldat_use.out ML.out ML_details.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps");
	
	# the entire script can be run multiple times by changing the value of NLOOP
	use constant NLOOP1 => 500;  
	$count = 0;
	while ( $count < NLOOP1)
	{
	    
	    ###########################                                                                                                                        
	    ### Simulate data       ###                                                                                                                       
	    ###########################        
	    
	    # number of sites in each gene
	    $genelength = 1;
	    
	    # the code below writes the input file for the "intra" program, which will simulate the gene trees from the coalescent model
	    # the only thing that needs to be changed is the name of the file containing the species tree:  mysptree4tax1 
	    open MYINFILE, ">ourinfile";
	    select MYINFILE;
	    print "begin coal;\n";
	    print "  ntaxa = 4;\n";
	    print "  species tree file = mysptree4tax_bl4;\n";
	    print "  intra no;\n";
	    print "  theta = 2;\n";
	    print " \n"; 
	    print "  gene tree file = simulated;\n";
	    print " \n"; 
	    print "  taxa names = Species1 Species2 Species3 Species4;\n";
	    print " \n";
	    print "  blstyle none;\n"; 
	    print " \n";
	    print "  outfile = myout;\n";
	    print " \n";
	    print "  ngtrees ", $ngenes, ";\n";   
	    print "  nstrees 1;\n";
	    my $range = 100000;
	    my $random_number = int(rand($range));
	    print "  seed1 = ", $random_number, "\n"; 
	    my $random_number = int(rand($range));
	    print "  seed2 = ", $random_number, "\n";
	    print " \n";
	    print "end;\n";
	    close MYINFILE;
	    
	    $nsites = $ngenes*$genelength;
	    system("cp ourinfile infile");
	    # the command below runs the intra program
	    system("./intra > outfile");
	    
	    # get the file with the gene trees ready for seq-gen:
	    system(qq( awk '{print "[$genelength]" \$1}' simtrees.dat > simtrees2.dat));
	    # The commands below call seq-gen to simulate sites along the gene trees; theta is set using the -s flag, e.g., the code below has theta=0.01
	    
	    # JC69 model
	    system("./seq-gen < simtrees2.dat > infile -q -mHKY -s 0.01 -l $nsites -p $ngenes");
	    system(qq( awk 'NR<2{print \$0; next}{print | "sort"}' infile > data.phy ));

	    ###########################
	    ### Get SVDQ estimate   ### 
	    ###########################
	    
	    system("/bin/cp -f settings_PrintFlat settings");
	    system("./PrintFlat");
	    system("cat svdq.out species_quartets.out > temp.svdq.out");
	    system("/bin/cp -f temp.svdq.out svdq.out");
	    
	    # check that there is enough data to compute scores
	    system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile));
	    system("cat temp.scoresfile valid.svdq > temp.valid.svdq");
	    system("/bin/cp -f temp.valid.svdq valid.svdq");

	    #################################                                                                                                                       
            ### Get SVDQ estimate - SNPs  ###                                                                                                                        
            #################################                                                                                                                        
                                                                                                                                                                     
            system("./PrintFlatSNPs");                                                                                                                               
	    system("cat svdq.snps.out species_quartets.out > temp.svdq.snps.out");                                                                               
	    system("/bin/cp -f temp.svdq.snps.out svdq.snps.out");                                                                                               
                                                                                                                                                                     
            # check that there is enough data to compute scores                                                                                                      
                                                                                                                                                                     
            system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile.snps));                                                                        
	    system("cat temp.scoresfile.snps valid.svdq.snps > temp.valid.svdq.snps");                                                                          
	    system("/bin/cp -f temp.valid.svdq.snps valid.svdq.snps");                                                                                           
                             
       ###########################                                                                                                            
       ### Get the MLE         ###                                                                                                                    
       ###########################        
       
       # Now run ssa to get the site pattern count for each set of 4 taxa of interest
       # (1,2,3,4)
       system("/bin/cp -f settings_ssa settings");
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("/bin/cp -f infile_data check1");
       system("./ssa > temp1.ssa");
       
       # (1,3,2,4)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp2.ssa");
       
       # (1,4,2,3)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp3.ssa");
       
       # (2,1,3,4)
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp4.ssa"); 
       
       # (2,3,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp5.ssa");  
       
       # (2,4,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp6.ssa");     
       
       # (3,1,2,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
        system("./ssa > temp7.ssa");    
       
       # (3,2,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data)); 
       system("./ssa > temp8.ssa");                    
       
       # (3,4,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp9.ssa");      
       
       # (4,1,2,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));       
       system("./ssa > temp10.ssa"); 
       
       # (4,2,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp11.ssa"); 
       
       # (4,3,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp12.ssa");

       # Use R code to find the ML tree
       system("R CMD BATCH ./get_lik.R");
       system("cat temp.R ML.out > temp.ML");
       system("/bin/cp -f temp.ML ML.out");
                                                                                                                                        
            $count = $count+1;                                                                                                                                       
                                                                                                                                                                    
        }                                                                                                                                                            
                                                                                                                                                                  
        system("paste ML.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps > alldat.out");     
	system("sed -f sed.tree.ind alldat.out > alldat_bl4_$ngenes.out");

    }

    system("zip alldat_bl4_0.01.zip alldat_bl4_*");  



printf "Running tree with branch length setting 5 \n\n\n";

    foreach my $ngenes (2000, 3000, 4000, 5000, 7500, 10000) {   

	# ngenes is a parameter for the number of genes to simulate
	# the theory is derived for the case where each site has its own gene tree; to simulate this, pick a large number of genes,
	# and give each length 1 (genelength parameter below) 
	# $ngenes = 200; 
	
	system("/bin/rm -f alldat.out alldat_use.out ML.out ML_details.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps");
	
	# the entire script can be run multiple times by changing the value of NLOOP
	use constant NLOOP1 => 500;  
	$count = 0;
	while ( $count < NLOOP1)
	{
	    
	    ###########################                                                                                                                        
	    ### Simulate data       ###                                                                                                                       
	    ###########################        
	    
	    # number of sites in each gene
	    $genelength = 1;
	    
	    # the code below writes the input file for the "intra" program, which will simulate the gene trees from the coalescent model
	    # the only thing that needs to be changed is the name of the file containing the species tree:  mysptree4tax1 
	    open MYINFILE, ">ourinfile";
	    select MYINFILE;
	    print "begin coal;\n";
	    print "  ntaxa = 4;\n";
	    print "  species tree file = mysptree4tax_bl5;\n";
	    print "  intra no;\n";
	    print "  theta = 2;\n";
	    print " \n"; 
	    print "  gene tree file = simulated;\n";
	    print " \n"; 
	    print "  taxa names = Species1 Species2 Species3 Species4;\n";
	    print " \n";
	    print "  blstyle none;\n"; 
	    print " \n";
	    print "  outfile = myout;\n";
	    print " \n";
	    print "  ngtrees ", $ngenes, ";\n";   
	    print "  nstrees 1;\n";
	    my $range = 100000;
	    my $random_number = int(rand($range));
	    print "  seed1 = ", $random_number, "\n"; 
	    my $random_number = int(rand($range));
	    print "  seed2 = ", $random_number, "\n";
	    print " \n";
	    print "end;\n";
	    close MYINFILE;
	    
	    $nsites = $ngenes*$genelength;
	    system("cp ourinfile infile");
	    # the command below runs the intra program
	    system("./intra > outfile");
	    
	    # get the file with the gene trees ready for seq-gen:
	    system(qq( awk '{print "[$genelength]" \$1}' simtrees.dat > simtrees2.dat));
	    # The commands below call seq-gen to simulate sites along the gene trees; theta is set using the -s flag, e.g., the code below has theta=0.01
	    
	    # JC69 model
	    system("./seq-gen < simtrees2.dat > infile -q -mHKY -s 0.01 -l $nsites -p $ngenes");
	    system(qq( awk 'NR<2{print \$0; next}{print | "sort"}' infile > data.phy ));

	    ###########################
	    ### Get SVDQ estimate   ### 
	    ###########################
	    
	    system("/bin/cp -f settings_PrintFlat settings");
	    system("./PrintFlat");
	    system("cat svdq.out species_quartets.out > temp.svdq.out");
	    system("/bin/cp -f temp.svdq.out svdq.out");
	    
	    # check that there is enough data to compute scores
	    system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile));
	    system("cat temp.scoresfile valid.svdq > temp.valid.svdq");
	    system("/bin/cp -f temp.valid.svdq valid.svdq");
	    
	    #################################                                                                                                                        
            ### Get SVDQ estimate - SNPs  ###                                                                                                                        
            #################################                                                                                                                        
                                                                                                                                                                     
            system("./PrintFlatSNPs");                                                                                                                               
            system("cat svdq.snps.out species_quartets.out > temp.svdq.snps.out");                                                                                   
            system("/bin/cp -f temp.svdq.snps.out svdq.snps.out");                                                                                                   
                                                                                                                                                                     
            # check that there is enough data to compute scores                                                                                                      
                                                                                                                                                                     
            system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile.snps));                                                                        
            system("cat temp.scoresfile.snps valid.svdq.snps > temp.valid.svdq.snps");                                                                               
            system("/bin/cp -f temp.valid.svdq.snps valid.svdq.snps");      

       ###########################                                                                                                            
       ### Get the MLE         ###                                                                                                                    
       ###########################        
       
       # Now run ssa to get the site pattern count for each set of 4 taxa of interest
       # (1,2,3,4)
       system("/bin/cp -f settings_ssa settings");
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("/bin/cp -f infile_data check1");
       system("./ssa > temp1.ssa");
       
       # (1,3,2,4)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp2.ssa");
       
       # (1,4,2,3)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp3.ssa");
       
       # (2,1,3,4)
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp4.ssa"); 
       
       # (2,3,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp5.ssa");  
       
       # (2,4,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp6.ssa");     
       
       # (3,1,2,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
        system("./ssa > temp7.ssa");    
       
       # (3,2,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data)); 
       system("./ssa > temp8.ssa");                    
       
       # (3,4,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp9.ssa");      
       
       # (4,1,2,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));       
       system("./ssa > temp10.ssa"); 
       
       # (4,2,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp11.ssa"); 
       
       # (4,3,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp12.ssa");

       # Use R code to find the ML tree
       system("R CMD BATCH ./get_lik.R");
       system("cat temp.R ML.out > temp.ML");
       system("/bin/cp -f temp.ML ML.out");
       
	    $count = $count+1;
	    
	}
	
	system("paste ML.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps > alldat.out");
	system("sed -f sed.tree.ind alldat.out > alldat_bl5_$ngenes.out");  
	
    }

    system("zip alldat_bl5_0.01.zip alldat_bl5_*");  



printf "Running tree with branch length setting 6 \n\n\n";

    foreach my $ngenes (2000, 3000, 4000, 5000, 7500, 10000) {   

	# ngenes is a parameter for the number of genes to simulate
	# the theory is derived for the case where each site has its own gene tree; to simulate this, pick a large number of genes,
	# and give each length 1 (genelength parameter below) 
	# $ngenes = 200; 
	
	system("/bin/rm -f alldat.out alldat_use.out ML.out ML_details.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps");
	
	# the entire script can be run multiple times by changing the value of NLOOP
	use constant NLOOP1 => 500;  
	$count = 0;
	while ( $count < NLOOP1)
	{
	    
	    ###########################                                                                                                                        
	    ### Simulate data       ###                                                                                                                       
	    ###########################        
	    
	    # number of sites in each gene
	    $genelength = 1;
	    
	    # the code below writes the input file for the "intra" program, which will simulate the gene trees from the coalescent model
	    # the only thing that needs to be changed is the name of the file containing the species tree:  mysptree4tax1 
	    open MYINFILE, ">ourinfile";
	    select MYINFILE;
	    print "begin coal;\n";
	    print "  ntaxa = 4;\n";
	    print "  species tree file = mysptree4tax_bl6;\n";
	    print "  intra no;\n";
	    print "  theta = 2;\n";
	    print " \n"; 
	    print "  gene tree file = simulated;\n";
	    print " \n"; 
	    print "  taxa names = Species1 Species2 Species3 Species4;\n";
	    print " \n";
	    print "  blstyle none;\n"; 
	    print " \n";
	    print "  outfile = myout;\n";
	    print " \n";
	    print "  ngtrees ", $ngenes, ";\n";   
	    print "  nstrees 1;\n";
	    my $range = 100000;
	    my $random_number = int(rand($range));
	    print "  seed1 = ", $random_number, "\n"; 
	    my $random_number = int(rand($range));
	    print "  seed2 = ", $random_number, "\n";
	    print " \n";
	    print "end;\n";
	    close MYINFILE;
	    
	    $nsites = $ngenes*$genelength;
	    system("cp ourinfile infile");
	    # the command below runs the intra program
	    system("./intra > outfile");
	    
	    # get the file with the gene trees ready for seq-gen:
	    system(qq( awk '{print "[$genelength]" \$1}' simtrees.dat > simtrees2.dat));
	    # The commands below call seq-gen to simulate sites along the gene trees; theta is set using the -s flag, e.g., the code below has theta=0.01
	    
	    # JC69 model
	    system("./seq-gen < simtrees2.dat > infile -q -mHKY -s 0.01 -l $nsites -p $ngenes");
	    system(qq( awk 'NR<2{print \$0; next}{print | "sort"}' infile > data.phy ));

	    ###########################
	    ### Get SVDQ estimate   ### 
	    ###########################
	    
	    system("/bin/cp -f settings_PrintFlat settings");
	    system("./PrintFlat");
	    system("cat svdq.out species_quartets.out > temp.svdq.out");
	    system("/bin/cp -f temp.svdq.out svdq.out");
	    
	    # check that there is enough data to compute scores
	    system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile));
	    system("cat temp.scoresfile valid.svdq > temp.valid.svdq");
	    system("/bin/cp -f temp.valid.svdq valid.svdq");
	    
	    #################################                                                                                                                        
            ### Get SVDQ estimate - SNPs  ###                                                                                                                        
            #################################                                                                                                                        
                                                                                                                                                                    
            system("./PrintFlatSNPs");                                                                                                                               
            system("cat svdq.snps.out species_quartets.out > temp.svdq.snps.out");                                                                                   
            system("/bin/cp -f temp.svdq.snps.out svdq.snps.out");                                                                                                   
                                                                                                                                                                     
            # check that there is enough data to compute scores                                                                                                      
                                                                                                                                                                    
            system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile.snps));                                                                        
            system("cat temp.scoresfile.snps valid.svdq.snps > temp.valid.svdq.snps");                                                                               
            system("/bin/cp -f temp.valid.svdq.snps valid.svdq.snps");      

       ###########################                                                                                                            
       ### Get the MLE         ###                                                                                                                    
       ###########################        
       
       # Now run ssa to get the site pattern count for each set of 4 taxa of interest
       # (1,2,3,4)
       system("/bin/cp -f settings_ssa settings");
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("/bin/cp -f infile_data check1");
       system("./ssa > temp1.ssa");
       
       # (1,3,2,4)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp2.ssa");
       
       # (1,4,2,3)
       system(qq( awk 'FNR==NR && /Species1/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp3.ssa");
       
       # (2,1,3,4)
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species3/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp4.ssa"); 
       
       # (2,3,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp5.ssa");  
       
       # (2,4,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species2/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp6.ssa");     
       
       # (3,1,2,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
        system("./ssa > temp7.ssa");    
       
       # (3,2,1,4)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species4/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data)); 
       system("./ssa > temp8.ssa");                    
       
       # (3,4,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species3/{a[0]=\$0; next}/Species4/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp9.ssa");      
       
       # (4,1,2,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species1/{a[1]=\$0; next}/Species2/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));       
       system("./ssa > temp10.ssa"); 
       
       # (4,2,1,3)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species2/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species3/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp11.ssa"); 
       
       # (4,3,1,2)                                                                                                            
       system(qq( awk 'FNR==NR && /Species4/{a[0]=\$0; next}/Species3/{a[1]=\$0; next}/Species1/{a[2]=\$0; next}/Species2/{a[3]=\$0 }END{printf "4 $nsites \\n%s\\n%s\\n%s\\n%s\\n",a[0],a[1],a[2],a[3]}' infile > infile_data));
       system("./ssa > temp12.ssa");

       # Use R code to find the ML tree
       system("R CMD BATCH ./get_lik.R");
       system("cat temp.R ML.out > temp.ML");
       system("/bin/cp -f temp.ML ML.out");

	    $count = $count+1;
	    
	}
	
	system("paste ML.out svdq.out valid.svdq svdq.snps.out valid.svdq.snps > alldat.out");
	system("sed -f sed.tree.ind alldat.out > alldat_bl6_$ngenes.out");  
	
}

    system("zip alldat_bl6_0.01.zip alldat_bl6_*");  


print "\n\nDone.\n";



