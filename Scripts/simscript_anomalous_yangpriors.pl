#!/usr/bin/perl

print "\n\n\n Simulation starting ..... \n\n\n";

foreach my $ngenes(50, 100, 150, 200, 300, 400) { 

   # ngenes is a parameter for the number of genes to simulate
   # the theory is derived for the case where each site has its own gene tree; to simulate this, pick a large number of genes,
   # and give each length 1 (genelength parameter below) 

   system("/bin/rm -f alldat.out alldat_use.out bpp_reps.out bppd_reps.out svdq.out valid.svdq x*");

   # the entire script can be run multiple times by changing the value of NLOOP
   use constant NLOOP1 => 100;  
   $count = 0;
   while ( $count < NLOOP1)
   {

       ###########################                                                                                                                                                        
       ### Simulate data       ###                                                                                                                                                        
       ###########################        
       
       # number of sites in each gene
       $genelength = 200;
       
       # the code below writes the input file for the "intra" program, which will simulate the gene trees from the coalescent model
       # the only thing that needs to be changed is the name of the file containing the species tree:  mysptree4tax1 
       open MYINFILE, ">ourinfile";
       select MYINFILE;
       print "begin coal;\n";
       print "  ntaxa = 4;\n";
       print "  species tree file = mysptree4tax_a1;\n";
       print "  intra no;\n";
       print "  theta = 2;\n";
       print " \n"; 
       print "  gene tree file = simulated;\n";
       print " \n"; 
       print "  taxa names = ^Species1 ^Species4 ^Species2 ^Species3;\n";
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
       system("./seq-gen < simtrees2.dat > infile -q -mHKY -s 0.025 -l $genelength -p $ngenes");

       ###########################     
       ### Run bpp             ###
       ###########################


       #### Yang priors #### 
       system(qq( sed 's/nloci = XX/nloci = $ngenes/g' a01_yang.ctl > a01.ctl));
       system("./bpp --c a01.ctl");
       system(qq(tail -n 1 infile.out.txt | grep "(A, B)" && echo "1" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(C, D)" && echo "1" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(A, C)" && echo "2" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(B, D)" && echo "2" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(A, D)" && echo "3" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(B, C)" && echo "3" > bpp.out));
       system("cat bpp.out bpp_reps.out > temp.bpp.out"); 
       system("/bin/cp -f temp.bpp.out bpp_reps.out");    

       #### Default priors ####
       system(qq( sed 's/nloci = XX/nloci = $ngenes/g' a01_master.ctl > a01.ctl));                                                        
       system("./bpp --c a01.ctl");                                                                                                     
       system(qq(tail -n 1 infile.out.txt | grep "(A, B)" && echo "1" > bppd.out));                                                      
       system(qq(tail -n 1 infile.out.txt | grep "(C, D)" && echo "1" > bppd.out));                                                      
       system(qq(tail -n 1 infile.out.txt | grep "(A, C)" && echo "2" > bppd.out));                                                      
       system(qq(tail -n 1 infile.out.txt | grep "(B, D)" && echo "2" > bppd.out));                                                      
       system(qq(tail -n 1 infile.out.txt | grep "(A, D)" && echo "3" > bppd.out));                                                      
       system(qq(tail -n 1 infile.out.txt | grep "(B, C)" && echo "3" > bppd.out));                                                      
       system("cat bppd.out bppd_reps.out > temp.bppd.out");                                                                               
       system("/bin/cp -f temp.bppd.out bppd_reps.out");        



       ###########################
       ### Get SVDQ estimate   ### 
       ###########################
       
       system("/bin/rm -f x*");
       system("/bin/rm -f infile2 infile3");
       system(qq(awk '/Species1/ {print \$2}' infile > species1_v1));
       system(qq(awk '/Species2/ {print \$2}' infile > species2_v1));
       system(qq(awk '/Species3/ {print \$2}' infile > species3_v1)); 
       system(qq(awk '/Species4/ {print \$2}' infile > species4_v1)); 
       system(qq(tr -d '\n' < species1_v1 > species1_v2));
       system(qq(tr -d '\n' < species2_v1 > species2_v2));
       system(qq(tr -d '\n' < species3_v1 > species3_v2));
       system(qq(tr -d '\n' < species4_v1 > species4_v2));
       system(qq(awk 'BEGIN {print "^Species1 "}' > species1));
       system("paste species1 species1_v2 > species1_v3");
       system(qq(awk 'BEGIN {print "^Species2 "}' > species2));
       system("paste species2 species2_v2 > species2_v3");   
       system(qq(awk 'BEGIN {print "^Species3 "}' > species3)); 
       system("paste species3 species3_v2 > species3_v3");   
       system(qq(awk 'BEGIN {print "^Species4 "}' > species4));
       system("paste species4 species4_v2 > species4_v3");   
       system("cat species*_v3 > infile2");
       system("/bin/rm -f data.phy");
       $nsites = $ngenes*$genelength;
       system(qq( awk 'BEGIN {print "4 ", $nsites}' > header));
       system("cat header infile2 > data.phy");
      
       system("/bin/cp -f settings_PrintFlat settings");
       system("./PrintFlat");
       system("cat svdq.out species_quartets.out > temp.svdq.out");
       system("/bin/cp -f temp.svdq.out svdq.out");


       # check that there is enough data to compute scores
       system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile));
       system("cat temp.scoresfile valid.svdq > temp.valid.svdq");
       system("/bin/cp -f temp.valid.svdq valid.svdq");
       
       $count = $count+1;
       
   }

   system("paste bpp_reps.out bppd_reps.out svdq.out valid.svdq > alldat.out");
   system("sed -f sed.tree.ind alldat.out > alldat_a1_025_200bp_yangpriors_$ngenes.out");
  
}


foreach my $ngenes(50, 100, 150, 200, 300, 400) { 

   # ngenes is a parameter for the number of genes to simulate
   # the theory is derived for the case where each site has its own gene tree; to simulate this, pick a large number of genes,
   # and give each length 1 (genelength parameter below) 

   system("/bin/rm -f alldat.out alldat_use.out bpp_reps.out bppd_reps.out svdq.out valid.svdq x*");

   # the entire script can be run multiple times by changing the value of NLOOP
   use constant NLOOP1 => 100;  
   $count = 0;
   while ( $count < NLOOP1)
   {

       ###########################                                                                                                                                                        
       ### Simulate data       ###                                                                                                                                                        
       ###########################        
       
       # number of sites in each gene
       $genelength = 1000;
       
       # the code below writes the input file for the "intra" program, which will simulate the gene trees from the coalescent model
       # the only thing that needs to be changed is the name of the file containing the species tree:  mysptree4tax1 
       open MYINFILE, ">ourinfile";
       select MYINFILE;
       print "begin coal;\n";
       print "  ntaxa = 4;\n";
       print "  species tree file = mysptree4tax_a1;\n";
       print "  intra no;\n";
       print "  theta = 2;\n";
       print " \n"; 
       print "  gene tree file = simulated;\n";
       print " \n"; 
       print "  taxa names = ^Species1 ^Species2 ^Species3 ^Species4;\n";
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
       system("./seq-gen < simtrees2.dat > infile -q -mHKY -s 0.025 -l $genelength -p $ngenes");

       ###########################     
       ### Run bpp             ###
       ###########################

       #### Yang priors ####
       system(qq( sed 's/nloci = XX/nloci = $ngenes/g' a01_yang.ctl > a01.ctl));
       system("./bpp --c a01.ctl");
       system(qq(tail -n 1 infile.out.txt | grep "(A, B)" && echo "1" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(C, D)" && echo "1" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(A, C)" && echo "2" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(B, D)" && echo "2" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(A, D)" && echo "3" > bpp.out));
       system(qq(tail -n 1 infile.out.txt | grep "(B, C)" && echo "3" > bpp.out));
       system("cat bpp.out bpp_reps.out > temp.bpp.out"); 
       system("/bin/cp -f temp.bpp.out bpp_reps.out");    

       #### Default priors ####
       system(qq( sed 's/nloci = XX/nloci = $ngenes/g' a01_master.ctl > a01.ctl));                                                                      
       system("./bpp --c a01.ctl");                                                                                                                   
       system(qq(tail -n 1 infile.out.txt | grep "(A, B)" && echo "1" > bppd.out));                                                                    
       system(qq(tail -n 1 infile.out.txt | grep "(C, D)" && echo "1" > bppd.out));                                                                    
       system(qq(tail -n 1 infile.out.txt | grep "(A, C)" && echo "2" > bppd.out));                                                                    
       system(qq(tail -n 1 infile.out.txt | grep "(B, D)" && echo "2" > bppd.out));                                                                    
       system(qq(tail -n 1 infile.out.txt | grep "(A, D)" && echo "3" > bppd.out));                                                                    
       system(qq(tail -n 1 infile.out.txt | grep "(B, C)" && echo "3" > bppd.out));                                                                    
       system("cat bppd.out bppd_reps.out > temp.bppd.out");                                                                                             
       system("/bin/cp -f temp.bppd.out bppd_reps.out");       


       ###########################
       ### Get SVDQ estimate   ### 
       ###########################
       
       system("/bin/rm -f x*");
       system("/bin/rm -f infile2 infile3");
       system(qq(awk '/Species1/ {print \$2}' infile > species1_v1));
       system(qq(awk '/Species2/ {print \$2}' infile > species2_v1));
       system(qq(awk '/Species3/ {print \$2}' infile > species3_v1)); 
       system(qq(awk '/Species4/ {print \$2}' infile > species4_v1)); 
       system(qq(tr -d '\n' < species1_v1 > species1_v2));
       system(qq(tr -d '\n' < species2_v1 > species2_v2));
       system(qq(tr -d '\n' < species3_v1 > species3_v2));
       system(qq(tr -d '\n' < species4_v1 > species4_v2));
       system(qq(awk 'BEGIN {print "^Species1 "}' > species1));
       system("paste species1 species1_v2 > species1_v3");
       system(qq(awk 'BEGIN {print "^Species2 "}' > species2));
       system("paste species2 species2_v2 > species2_v3");   
       system(qq(awk 'BEGIN {print "^Species3 "}' > species3)); 
       system("paste species3 species3_v2 > species3_v3");   
       system(qq(awk 'BEGIN {print "^Species4 "}' > species4));
       system("paste species4 species4_v2 > species4_v3");   
       system("cat species*_v3 > infile2");
       system("/bin/rm -f data.phy");
       $nsites = $ngenes*$genelength;
       system(qq( awk 'BEGIN {print "4 ", $nsites}' > header));
       system("cat header infile2 > data.phy");
      
       system("/bin/cp -f settings_PrintFlat settings");
       system("./PrintFlat");
       system("cat svdq.out species_quartets.out > temp.svdq.out");
       system("/bin/cp -f temp.svdq.out svdq.out");


       # check that there is enough data to compute scores
       system( qq(awk '{sum += \$2} END{print sum}' scoresfile > temp.scoresfile));
       system("cat temp.scoresfile valid.svdq > temp.valid.svdq");
       system("/bin/cp -f temp.valid.svdq valid.svdq");
       
       $count = $count+1;
       
   }

   system("paste bpp_reps.out bppd_reps.out svdq.out valid.svdq > alldat.out");
   system("sed -f sed.tree.ind alldat.out > alldat_a1_025_1000bp_yangpriors_$ngenes.out");
  
}


print "\n\nDone.\n";



