### Spectra/21BRCA.*.tsv:

1. Download wget ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerExtractor/Example_data/21BRCA.zip,
see [this link](https://osf.io/t6j7u/wiki/2.%20Quick%20Start%20Example/) for more details.
2. Unzip the zip file to your destination folder:
```
unzip 21BRCA.zip
```
3. Run SigProfilerMatrixGenerator to obtain multiple catalog formats. 
See [this link](https://osf.io/s93d5/wiki/6.%20Quick%20Start%20Example/) for more details:

```
$ python3

>>> ## Installation of SigProfilerMatrixGenerator
>>> pip install SigProfilerMatrixGenerator


>>> ## Installation of reference genome (skip if previously installed)
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh37')

>>> ## Generate SigProfiler-formatted catalog matrices from example 21BRCA vcf files:
>>> from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
>>> matrices = matGen.SigProfilerMatrixGeneratorFunc("vcftest", "GRCh37", "/<Working_directory>/21BRCA/21BRCA_vcf",plot=True,tsb_stat=True)
>>>

```

4. Rename vcftest.*.all into vcftest.*.tsv, copy them to "./Spectra"

5. Copy 21BRCA.SBS96.tsv and 21BRCA.DBS78.tsv into "tests/testthat/testdata/SigPro-Cat" 


## Spectra/Strelka.ID.GRCh37.s1.*.tsv:

1. Copy <ICAMS_HOME>/tests/testthat/testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf into an empty folder <FOLDER_PATH>

2. Run SigProfilerMatrixGenerator to obtain multiple indel catalog formats. 

```
$ python3

>>> ## Installation of SigProfilerMatrixGenerator
>>> pip install SigProfilerMatrixGenerator


>>> ## Installation of reference genome (skip if previously installed)
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh37')

>>> ## Generate SigProfiler-formatted catalog matrices from example 21BRCA vcf files:
>>> from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
>>> matrices = matGen.SigProfilerMatrixGeneratorFunc("Strelka.ID.GRCh37.s1", "GRCh37", "<FOLDER_PATH>/",plot=True)
>>>

```

3. Rename Strelka.ID.GRCh37.s1.*.all as Strelka.ID.GRCh37.s1.*.tsv, and copy them into 
<ICAMS_HOME>/data-raw/data/SigPro/Spectra.

4. Copy Strelka.ID.GRCh37.s1.ID83.tsv into "tests/testthat/testdata/SigPro-Cat".