# First things first: Get PAF latest version
Log into gridui and load root:

    source /gpfs/csic_projects/cms/sw/ROOT/current/root/bin/thisroot.sh

Now clone the git PAF repository (into your home or a suitable directory)

    git clone https://github.com/PROOF-Analysis-Framework/PROOF-Analysis-Framework

Move to the new directory and install PAF

    cd PROOF-Analysis-Framework
    cd build
    make
    make install
  
To follow the basics steps collected here no extra knowledge is required, but to understand how the latest PAF version works is strongly recommended, and essential if you want to edit the files here and make your own implementations. For this task you should check the git repository:

	https://github.com/PROOF-Analysis-Framework/
	
And the simple and well-explained tutorials there or at the site:

	http://www.hep.uniovi.es/PAF/
	
# Now get this code (and make it your own!)
If not already done, log into your github account and click "Fork" at the top-right corner of this
page to get your own copy of this repository.

Now, go back to your terminal where you logged into gridui and do the following. This is just to make git work properly, specially when submitting your own changes (it avoids using a ssh key):
  Go to another directory if you like, and set the CMS environment:

     source /cvmfs/cms.cern.ch/cmsset_default.sh

  And the architecture:

     export SCRAM_ARCH=slc6_amd64_gcc481

  Now you can choose your favorite CMSSW release.

     cmsrel CMSSW_7_2_0
     cd CMSSW_7_2_0/src
     cmsenv
     
After doing this go back to your working directory and get the material:

    git clone https://github.com/YOUR_GITHUB_USERNAME/PAF-MuonAnalyzer

Note that you need to place YOUR OWN GITHUB USERNAME in the path, because you will be cloning your fork of the repository.

# Run the code
First, go to the directory where you cloned the repository

	cd PAF-MuonAnalyzer
	
Now source the PAF setup

	source wathever_path_leading_up_to_PAF/PROOF-Analysis-Framework/PAF_setup.sh
	
Take a brief look to the file RunMuonAnalyzer.C. It's the one you will use to run the PAF Project with root CINT interpreter. Note
