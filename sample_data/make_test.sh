#! /bin/bash

wget https://s3.amazonaws.com/paladin.aligner/test.fq
wget https://s3.amazonaws.com/paladin.aligner/paladin_test.faa
../paladin index -r2 paladin_test.faa
../paladin align -t4 -u2 paladin_test.faa test.fq > paladin_uniprot_report.txt

if [ -f paladin_uniprot_report.txt ];
        then
            echo "\n\nPALADIN HAS BEEN SUCCESSFULLY INSTALLED\n\n"
        else
            echo "\n\nOOPS: SOMETHING WENT WRONG WITH INSTALLATION, OR YOU ARE NOT CONNECTED TO THE INTERNET\n\n"
fi

rm paladin_test.faa*
