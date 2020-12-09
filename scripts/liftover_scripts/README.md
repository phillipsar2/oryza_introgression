### Step 1: prepare the genomes for alignments
lastdb.sh


### Step 2: generate the LAST alignment, process the laignment, merge files, and generate a net and chain file
make_alignments.sh


### Step 3: subset chain file and merge subset to form final liftOver chain. Then liftover the file.
chain_liftover.sh
