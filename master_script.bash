#bash get_genomes.bash -f K.pneumoniae &> out_record.txt
#bash Kallisto.sh -f K.pneumoniae &> kallisto_record.txt

bash get_genomes.sh -f test_dir &> test_records.txt



# Notes for getting panoct to work:
# Add .pl files to path
# Install Anaconda
# Install perl and bioperl through anaconda using this code for perl:
#   https://stackoverflow.com/questions/55227065/x86-64-conda-cos6-linux-gnu-gcc-not-found
# Install perl packages Bio::SeqIO and Term::ReadKey a la https://www.putorius.net/how-to-install-perl-modules-with-cpan.html
# Set PERL5LIB variable appropriately in .bashrc a la https://stackoverflow.com/questions/58290190/how-to-fix-perl-from-anaconda-not-installing-bioperl-bailing-out-the-installat
