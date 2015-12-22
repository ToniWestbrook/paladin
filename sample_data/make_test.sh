#! /bin/bash

curl -O https://s3.amazonaws.com/paladin.aligner/test.fq
curl -O https://s3.amazonaws.com/paladin.aligner/paladin_test.faa
../paladin index -r3 paladin_test.faa
../paladin prepare -r1 -f paladin_test.faa
../paladin align -t4 paladin_test.faa test.fq -o test

if [ -s test_uniprot.tsv ];
        then
            echo "\n\nPALADIN HAS BEEN SUCCESSFULLY INSTALLED\n\n"
        else
            echo "\n\nOOPS: SOMETHING WENT WRONG WITH INSTALLATION, OR YOU ARE NOT CONNECTED TO THE INTERNET\n\n"
fi

rm paladin_test.faa*
