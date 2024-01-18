# RGP

#Generating shuffled sequences 
To create shuffled sequence sets use test_hg38_seqs_196kb.fasta and biasaway (sample command below)

##Local Shuffling
`biasaway w -f /project/st-cdeboer-1/iluthra/enformer_random_DNA/enformer_data/test_hg38_seqs_196kb.fasta -k 3 -s 100 -w 200 > /project/st-cdeboer-1/iluthra/enformer_random_DNA/enformer_data/test_hg38_seqs_196kb_trinuc_local.fasta`

##Global Shuffling
`biasaway k -f /project/st-cdeboer-1/iluthra/enformer_random_DNA/enformer_data/test_hg38_seqs_196kb.fasta -k 3 > /project/st-cdeboer-1/iluthra/enformer_random_DNA/enformer_data/test_hg38_seqs_196kb_trinuc.fasta`

#References
Khan, A., Riudavets Puig, R., Boddie, P. & Mathelier, A. BiasAway: command-line and web server to generate nucleotide composition-matched DNA background sequences. Bioinformatics 37, 1607–1609 (2021).

Avsec, Ž. et al. Effective gene expression prediction from sequence by integrating long-range interactions. Nat Methods 18, 1196–1203 (2021).