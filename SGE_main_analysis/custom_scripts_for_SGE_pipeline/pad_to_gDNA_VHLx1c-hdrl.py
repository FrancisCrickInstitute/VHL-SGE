#cDNA_to_gDNA_SGE_pipeline.py

#arguments -- [1], input file (designed to be anything pre or post merger);  [2], output file (will be an unzipped fastq) [3] add_on seq
#run in  no_Ns folder
#python /camp/lab/findlayg/home/users/bucklem/_Scripts/_VHL/Analysis_pipelines/pad_to_gDNA_VHLx1c-hdrl.py VHLx1c2r-hdrl.merged.fastq VHLx1c2-hdrl.merged.fastq GAGAT

import gzip
import sys
from itertools import islice
out_file = open(sys.argv[2],'w')
add_on = sys.argv[3] #sequence to pad on 3' end of read

#read in the editing data for all amplicons
line_count = 0
reads = 0
with open(sys.argv[1], 'r') as my_reads:
    while True:
        my_fastq = list(islice(my_reads, 4))
        if not my_fastq:
            break
        else:
            reads+=1
            seq = my_fastq[1].strip()
            qs = my_fastq[3].strip()
            new_seq = seq+add_on
            new_qs = qs+len(add_on)*'H'
            out_file.write(''.join([my_fastq[0],new_seq,'\n',my_fastq[2],new_qs,'\n']))
        
print sys.argv[1], 'finished processing.'
print reads, 'total reads processed'





