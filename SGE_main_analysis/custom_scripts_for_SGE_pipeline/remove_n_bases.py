#remove_n_bases.py
#set up calls to this script using run_remove_n_bases.py script in directory with seqprep merged output files (.fastq.gz)
#arguments -- [1], input file (designed to be anything pre or post merger);  [2], output file (will be an unzipped fastq)

import gzip
import sys
from itertools import islice
out_file = open(sys.argv[2],'w')

with gzip.open(sys.argv[1], 'r') as my_reads:
    reads = 0
    reads_w_n = 0
    while True:
        my_fastq = list(islice(my_reads, 4))
        if not my_fastq:
            break
        else:
            reads+=1
            seq = str(my_fastq[1])
            if seq.find("N") != -1:
                reads_w_n += 1
            else:
                out_fastq = str(''.join(my_fastq))
                out_file.write(out_fastq)
    print(str(sys.argv[1])+' finished processing.')
    print(str(reads_w_n)+' reads with N base')
    print(str(reads)+' total reads')