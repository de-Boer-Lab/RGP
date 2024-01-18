

# Model Implementation codes:

All the modules in this repo are taken from the goggle colab [here](https://colab.research.google.com/github/deepmind/deepmind_research/blob/master/enformer/enformer-usage.ipynb). 

**Enformer.py**:

Made modifications in `EnformerScoreVariantsRaw` class, to return a list of reference and alternate allele predictions as opposed to a difference between the mean for each feature between the two alleles.

**FastaExt.py**:

Made the following modifications in `variant_centered_sequences`function:
- Added `tss` argument that takes a genomic coordinate to centers the input sequence at that position before making predictions. This is to support TSS based prediction for non coding mutations. By default, this is set to be `None` and the centering is done based on the position in the vcf file.