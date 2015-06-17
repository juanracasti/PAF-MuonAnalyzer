# First things first: Get PAF latest version
Log into gridui and load root:

	ssh ...
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

After doing this go back to your working directory and get the material:

    git clone https://github.com/YOUR_GITHUB_USERNAME/PAF-MuonAnalyzer

Note that you need to place YOUR OWN GITHUB USERNAME in the path, because you will be cloning your fork of the repository.

# Run the code
First, go to the directory where you cloned the repository

	cd PAF-MuonAnalyzer
	
Now source the PAF setup

	source wathever_path_leading_up_to_PAF/PROOF-Analysis-Framework/PAF_setup.sh
	
Take a brief look at the file RunMuonAnalyzer.C. It's the one you will use to run the PAF Project with root CINT interpreter. Note that there is a variable called 'signal' that can take several values depending on the Trees to load. So to run the code, you just do:

	root -l
	.x RunMuonAnalyzer.C("SIGNAL")
	
where signal is one of the possible values. For example try:

	.x RunMuonAnalyzer.C("MC_DY_PU20bx25")

# Edit the code
Read the PAF documentation first! Despite that, the files are well documented and should be quite self-explanatory. Try to mantain that philosophy :)

If you want to keep up to date with the origin repository, before doing any edits, first get the latest changes in the repository, if any, to avoid merging conflicts. Always keep in mind the possible differences with your current copy of the repository.

    git pull https://github.com/juanracasti/PAF-MuonAnalyzer.git

Always check if the code compiles and produces the expected results before submitting any changes.

# Submit your changes (if you want!)
You may not want to do this if you just want to keep your own copy of the repository and then do your own stuff locally. However it is always good to maintain some sort of version control!
Always check if the code compiles and produces the expected results before submitting any changes.
Now, if you do want to commit your changes, do the following trick. This is just to make git work properly when submitting your own changes (it avoids using a ssh key):
  Go to another directory if you like, and set the CMS environment:

     source /cvmfs/cms.cern.ch/cmsset_default.sh

  And the architecture:

     export SCRAM_ARCH=slc6_amd64_gcc481

  Now you can choose your favorite CMSSW release.

     cmsrel CMSSW_7_4_0
     cd CMSSW_7_2_4/src
     cmsenv

And return to the current directory. Then commit your changes:

    git status
    git add <filepattern>
    git commit -m 'Modified'
    git push

Now proceed with caution. Only follow the steps below if you have made contributions that you think are universal for the code or that all the users could benefit from them, like a new selector or an improvement of some sort. Do not do it if what you have implemented is ONLY useful for yourself.

Open a new tab in your browser with your copy of the repository:

    https://github.com/your_github_username/MuonTreeAnalyzer

Open a new pull request with your commits and a brief description of the changes. Eventually it will be accepted and integrated into the origin repository.
