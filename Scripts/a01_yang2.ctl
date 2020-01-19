          seed =  12345

       seqfile = infile
      Imapfile = infile.Imap.txt
       outfile = infile.out.txt
      mcmcfile = infile.mcmc.txt

  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
         speciestree = 1

*   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2:uniformSLH; 3: uniformSRooted

  species&tree = 4  A  B  C  D
                    1  1 1  1
                   ((A, B), (C, D));
                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = XX  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes,     0:no)?

    thetaprior = 3 0.02 E  # invgamma(a, b) for theta
      tauprior = 3 0.024    # invgamma(a, b) for root tau &      Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for     GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars,         Genetrees
        burnin = 400
      sampfreq = 2
       nsample = 1500
