![Definition of chipper](/pics/chipper.png?raw=true "Definition of chipper")

Current status of the [chipper build](https://travis-ci.org/massie/chipper)
[![Build Status](https://travis-ci.org/massie/chipper.svg?branch=master)](https://travis-ci.org/massie/chipper)

## Quick Install

Chipper is comprised for a frontend, written in C, and a backend written in a [Jupyter notebook](/py/create_chipper_model.ipynb) in Python.

The frontend runs on MacOS and Linux and requires a BLAS library, e.g. [OpenBLAS](http://www.openblas.net/) and `cmake` to build. For MacOS brew users, you can satisfy these requirements by running `brew install cmake blas`. Chipper supports both the `clang` and `gcc` compilers.

```
$ git clone https://github.com/massie/chipper.git
$ cd chipper
$ mkdir build
$ cd build
$ cmake ..
$ make all test install
```

By default, chipper will install to your `/usr/local/` directory. If you want to install chipper in a different directory, change the `cmake` commandline above, e.g.

```
$ cmake -DCMAKE_INSTALL_PREFIX:PATH=/your/alternative/prefix ..
```
For example,

```
$ cmake -DCMAKE_INSTALL_PREFIX:PATH=/tmp/example ..
$ make all test install
$ find /tmp/example
/tmp/example
/tmp/example/bin
/tmp/example/bin/chipper
/tmp/example/include
/tmp/example/include/chipper.h
/tmp/example/lib
/tmp/example/lib/libchipper.a
/tmp/example/share
/tmp/example/share/chipper
/tmp/example/share/chipper/lr.model
/tmp/example/share/chipper/svc.model
/tmp/example/share/chipper/test.fa
```

Now that you have `chipper` installed, you can use the included test FASTA file, `test.fa` to ensure that it's working, e.g.

```
$ /tmp/example/bin/chipper -i /tmp/example/share/chipper/test.fa | gunzip
@gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKKPAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKKK
+gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
!!!!!!!!!!++-++-+--++++++++----+++++----+--++++++++-+------++++++-++++-+-+-+---++++-+++++-+--+++++-------++++++++----+-+++++-++++++----+++++++++++--++-+++--++++++++-+-+++++--+++++++--+-++-+++++++++++-++++-++++!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!!-+-++++++++-+-+--+-----+-++++----++-----+-+-+-++-++-+-++++++++-----++++++----++++-+--+-++---+++++-+++++++++++++--------++++-+-+++--+--+++++++--+--++--+--++-++-+-++-++--+-+++-++++----+++-+-+-++++++-!!!!!!!!!!
```

To learn more about the chipper commandline, run `chipper -h`. To learn more about the output, continue reading.
## Chipper Backend

You do **not** need to run the `chipper` backend in order to use `chipper`. The backend generates the models and some C code that the frontend uses. It's shared in the repo for transparency and to allow you to customize the backend. You can [view the jupyter notebook](/py/create_chipper_model.ipynb) in this repository including all output. 

If you want to run the backend notebook, you will need a copy of the Swiss-Prot database at `swiss_db_path`, `blast` installed, and list of epitopes exported from IEDB.

The notebook has the following high-level steps:

1. The IEDB csv file with epitopes is loaded and only epitopes of length [`min_epitope_len`, `max_epitope_len`] are used.
2. Each epitope has a protein ID which is cross-referenced against Swiss-Prot database. Each epitope seqence and location are double-checked. 
3. Any sequences that were part of the Saxova dataset are removed in order to test the model (see *Appendix A* of the paper "[Predicting proteasomal cleavage sites: a comparison of available methods](http://intimm.oxfordjournals.org/content/15/7/781.full)")
4. For each epitope, one cleavage and one non-cleavage example each of length `generated_sample_len` (see more details below).
5. The amino acid sequences are converted into features using the amino acid properties found in the paper "[A new set of amino acid descriptors and its application in peptide QSARs](https://www.ncbi.nlm.nih.gov/pubmed/15895431)". Each amino acids is defined by 50 properties. If the `generated_sample_len` is `20`, then each sample is comprised of `20 * 50 = 1000` features. *Note that there are currently no metrics for selenocysteine. Chipper will print a warning and continue processing.*
6. The [scikit learn](http://scikit-learn.org/) [MinMaxScaler](http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html) is used to scale the training and test datasets to values from `0` to `1`. Note that later in the model creation process, the scale values determined by the `MinMaxScaler` for each amino acid properties are used in the C code generation step to create a scaled amino acid property lookup table.
7. Two [liblinear](http://www.csie.ntu.edu.tw/~cjlin/liblinear/) models are created using L2-regularized logistic regression (the default) and L2-regularized L2-loss support vector classification (primal). The optimal `C` parameter used is found using liblinear's hyperparameter tuning feature with 5-fold cross-validation.
8. The performance of the model is measured against the Saxova dataset (see the performance numbers below).


### Cleavage and Non-Cleavage Samples

The C-terminus of each epitope is considered a proteasome cleavage position. The midpoint between the N-terminus and C-terminus is considered a non-cleavage position. For example, if a protein that's 500 amino acids long has an epitope that starts at position 240 and ends at position 254. If the `generated_sample_len` is `20`, meaning that 10 flanking amino acids are taken from the protein to create a sample, then the following protein positions will be used to create the samples.


|               | Start Position (N-terminus) | End Position (C-terminus) |
| ------------- | -------------- | -------------|
| Epitope       | **240**  | **254** |
| Cleavage Sample  | 254-10=**244**  | 254+10=**264** |
| Non-Cleavage Sample | (240+254)/2-10=**237** | (240+254)/2+10=**257** |

## Chipper Frontend

To see the `chipper` commandline options, e.g.

```
chipper -h
Usage: chipper [-hvp] [-i <fasta file>] [-m <liblinear model>] [-o <fastq file>] [-c <cutoff>]
 -h, --help                Display this help and exit
 -v, --version             Display version info and exit
 -i, --input=<fasta file>  FASTA file with protein(s) (default: stdin)
 -m, --model=<liblinear model> (default:/usr/local/share/chipper/lr.model)
 -o, --output=<fastq file> FASTQ file with predictions (default: stdout)
 -c, --cutoff=<cutoff>     Cutoff in range (0,1) for probability models (default: cutoff with highest MCC)
 -p, --probs               Output probabilities to FASTQ
```

The `chipper` frontend is a command-line tool takes proteins in FASTA format and output a gzip FASTQ file with cleavage positions expressed. The quality values are expressed using the following characters: `!0123456789+-`. The character `!` means that `chipper` is unable to make a prediction about cleavage because there is not enough flanking amino acids to decide. A '0' means <10% probability of cleavage, '1' means a 10 to 20% exclusive probability, etc. with '9' meaning 90 to 100% cleavage. For classification output, a '-' means no cleavage and '+' marks a cleavage position.

There is a small FASTA file, `test.fa`, which is installed in `${CMAKE_INSTALL_PREFIX}/share/chipper/test.fa` that will be used for these examples.

### Predict cleavage positions using the default Logistic Regression (LR) model

```
$ chipper -i /usr/local/share/chipper/test.fa | gunzip
@gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKKPAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKKK
+gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
!!!!!!!!!!++-++-+--++++++++----+++++----+--++++++++-+------++++++-++++-+-+-+---++++-+++++-+--+++++-------++++++++----+-+++++-++++++----+++++++++++--++-+++--++++++++-+-+++++--+++++++--+-++-+++++++++++-++++-++++!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!!-+-++++++++-+-+--+-----+-++++----++-----+-+-+-++-++-+-++++++++-----++++++----++++-+--+-++---+++++-+++++++++++++--------++++-+-+++--+--+++++++--+--++--+--++-++-+-++-++--+-+++-++++----+++-+-+-++++++-!!!!!!!!!!
```

Since `chipper` needs ten flanking amino acids in order to predict cleavage points, the first and last ten positions of each protein has '!!!!!!!!!!`.

To see the default cutoff that `chipper` used with the logistic regression model, see the version output, e.g.

```
$ chipper -v
chipper version: 0.1.0
Best cutoff is 0.35 (MCC= 0.62)
```

#### Predict cleavage positions using the default LR model and a custom cutoff.

For example, if you want to only mark a position as a cleavage position when the probability of cleavage is `>=` 80%, e.g.

```
$ chipper -i /usr/local/share/chipper/test.fa -c 0.8 | gunzip
@gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKKPAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKKK
+gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
!!!!!!!!!!-+--+------++++++-----++++-------++-+--++-+------+---+--+-+-----------++--+--++----++--+---------+-++++-------+-+----+-------+++--+--------+-++----+-++-------+--+-----++-+--+------+-+--+--+-+----+--+!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!!------+----------------+-++-+----++-----+-+-+--+---------++++-------++-++------++-------+------++-----+--+++++---------+++--+--++-----+--++++------+-----++-+--+------------+--+++-----+----+-+++++--!!!!!!!!!!
```

When your compare the output to the previous (which used a cutoff of 0.35), you see less positions marked as cleavage positions. Using a higher probability reduces the number of false positives but also increases the number of false negatives.

#### Output probabilities of cleavage using the default Logistic Regression model

Instead of '-/+' notation, you can see the probabilities as each position, using the `-p/--probs` flag, e.g.

```
$ chipper -i /usr/local/share/chipper/test.fa **-p** | gunzip
@gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKKPAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKKK
+gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
!!!!!!!!!!5915815204798999910017999800004019949759938020020976787195960505060017885084399050089768111200066949998010040697961759776012389864957554014818962159698565162497480277398791281761349596793383954619579!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!!25135696477073412600021929978003099220019180804917707056599894000025894980100469915307248010556883644794589989400200001899708168910400846998920400481172098095391741740172769049890000597250938999841!!!!!!!!!!
```

#### Predict cleavage positions using support vector classification (SVC)

```
$ chipper -i /usr/local/share/chipper/test.fa -m /usr/local/share/chipper/svc.model  | gunzip
@gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKKPAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKKK
+gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
!!!!!!!!!!++-++-+---+++++++----+++++----+--++-+++++-+------++++++-++++-+-+-+---++++-++-++-+--+++++-------++++++++----+-+++++-++++++----++++++++++---++-+++--++++++++-+-+++-+--++-++++--+-++--+++++++--+-++-+-++++!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!!-+--+++++++-+----+-----+-++++----++-----+-+-+--+-++-+-++++++++-----+++-++-----+++-+--+-++---+++++-+--++++++++++--------++++-+-+++-----+++++++-----++--+--++-++-+-++-+---+-+++-++++----+++-+-+-+++++--!!!!!!!!!!
```

The SVC model runs just as fast as the LR model but has slightly different sensitivity, specificity, and MCC metrics (see below).

In fact, you can generate your own `liblinear` model and pass it on the command-line. 

## Performance Comparison

Chipper compares favorably against other proteasomal cleavage detection systems. Metrics reported in "[Predicting proteasomal cleavage sites: a comparison of available methods](http://intimm.oxfordjournals.org/content/15/7/781.full)" compared with chipper. Here's the performance using the same dataset using the same test dataset.

| Method | Sensitivity (recall) | Specificity | MCC    |
|--------|-------------|-------------|--------|
|[PAProC](http://www.paproc.de/)| 46.4 |64.7 |0.10 |
|[FragPredict](http://www.mpiib-berlin.mpg.de/MAPPP/fragpredict.html) | 72.1 | 41.4 | 0.12 |
|[NetChop](http://www.cbs.dtu.dk/services/NetChop/) 1.0 |34.4 |91.4 |0.31 |
|[NetChop](http://www.cbs.dtu.dk/services/NetChop/) 2.0 |57.4 |76.4 |0.32 |
|chipper (LR) | **77.4** | **85.2** | **0.62** |
|chipper (SVC) | **78.2** | **79.0** | **0.57** |

The SVC model uses a linear kernel which makes interpreting the effect of the hydrophobic, steric, and electronic properties possible, while the performance of the logistic regression model has better specificity.

The `chipper` command is also fast. For example, you can predict cleavage locations of 98,016 proteins in about 47 seconds on a Mid 2014 MacBook Pro. A rate of over 2,000 proteins per second.

```
$ time chipper -i /workspace/chipper_data/training/proteins.fasta -o /dev/null
real    0m47.451s
user    0m47.337s
sys     0m0.049s
```

The ROC curve for the chipper logistic regression model shows an area under the curve of 0.87. The chipper backend calculates the best cutoff using the Matthews correlation coefficient (MCC). You can find this value in the version information for `chipper`, e.g. `chipper -v`.

![ROC Curve for Logistic Regression](/pics/lr_roc_curve.png?raw=true "ROC Curve for Logistic Regression")


## Contributing

Contributions welcomed! If you find an issue, feel free to submit a pull request or [file an issue](https://github.com/massie/chipper/issues).

## Biological Information

[![MHC](/pics/MHC_Class_I_processing.png?raw=true "MHC")](https://en.wikipedia.org/wiki/Major_histocompatibility_complex)

From [Major histocompatibility complex class I binding predictions as a tool in epitope discovery](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2913210/)...

> Cytotoxic T lymphocytes (CTLs) are the effector cells of the adaptive immune response that deal with infected, or malfunctioning, cells. Whereas intracellular pathogens are shielded from antibodies, CTLs are endowed with the ability to recognize and destroy cells harbouring intracellular threats. This obviously requires that information on the intracellular protein metabolism (including that of any intracellular pathogen) be translocated to the outside of the cell, where the CTL reside. To this end, the immune system has created an elaborate system of antigen processing and presentation. During the initial phase of antigen processing, peptide antigens are generated from intracellular pathogens and translocated into the endoplasmic reticulum. **In here, these peptide antigens are specifically sampled by major histocompatibility complex (MHC) class I molecules and then exported to the cell surface, where they are presented as stable peptide: MHC I complexes awaiting the arrival of scrutinizing T cells.** Hence, identifying which peptides are able to induce CTLs is of general interest for our understanding of the immune system, and of particular interest for the development of vaccines and immunotherapy directed against infectious pathogens, as previously reviewed. Peptide binding to MHC molecules is the key feature in cell-mediated immunity, because it is the peptide–MHC class I complex that can be recognized by the T-cell receptor (TCR) and thereby initiate the immune response. The CTLs are CD8+ T cells, whose TCRs recognize foreign peptides in complex with MHC class I molecules. In addition to peptide binding to MHC molecules, several other events have to be considered to be able to explain why a given peptide is eventually presented at the cell surface. Generally, an immunogenic peptide is generated from proteins expressed within the presenting cell, and peptides originating from proteins with high expression rate will normally have a higher chance of being immunogenic, compared with peptides from proteins with a lower expression rate. There are, however, significant exceptions to this generalization, e.g. cross-presentation, but this will be ignored in the following. In the classical MHC class I presenting pathway (*see image on right*) proteins expressed within a cell will be degraded in the cytosol by the protease complex, named the **proteasome**. **The proteasome digests polypeptides into smaller peptides 5–25 amino acids in length and is the major protease responsible for generating peptide C termini**. Some of the peptides that survive further degradation by other cytosolic exopeptidases can be bound by the transporter associated with antigen presentation (TAP), reviewed by Schölz et al. This transporter molecule binds peptides of lengths 9–20 amino acids and transports the peptides into the endoplasmic reticulum, where partially folded MHC molecules [in humans called human leucocyte antigens (HLA)], will complete folding if the peptide is able to bind to the particular allelic MHC molecule. The latter step is furthermore facilitated by the endoplasmic-reticulum-hosted protein tapasin. Each of these steps has been characterized and their individual importance has been related to final presentation on the cell surface.

