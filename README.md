![Definition of chipper](/pics/chipper.png?raw=true "Definition of chipper")

Current status of the [chipper build](https://travis-ci.org/massie/chipper)
[![Build Status](https://travis-ci.org/massie/chipper.svg?branch=master)](https://travis-ci.org/massie/chipper)

Chipper predicts cleavage sites of the human proteasome. See the biological information at the end of this readme for more details. Chipper takes FASTA/FASTQ files for input and outputs the predictions in FASTQ or NetChop 3.0 format.

## Coordinates

Amino acids have an amino group (the N-terminus side) and a carboxylic acid group (the C-terminus side). Chipper predicts whether or not the *C-terminus* at a specific location will be cleaved by proteasomes.

## Quick Install

Chipper is comprised of a frontend, written in C, and a backend written in a [Jupyter notebook](/py/create_chipper_model.ipynb) in Python.

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
!!!!!!!!!++-++-+--++++++++----+++++----+--++++++++-+------++++++-++++-+-+-+---++++-+++++-+--+++++-------++++++++----+-+++++-++++++----+++++++++++--++-+++--++++++++-+-+++++--+++++++--+-++-+++++++++++-++++-++++-!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!-+-++++++++-+-+--+-----+-++++----++-----+-+-+-++-++-+-++++++++-----++++++----++++-+--+-++---+++++-+++++++++++++--------++++-+-+++--+--+++++++--+--++--+--++-++-+-++-++--+-+++-++++----+++-+-+-++++++--!!!!!!!!!!
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
Usage: chipper [-hvnsp] [-i <fasta file>] [-m <liblinear model>] [-o <file>] [-c <cutoff>]
 -h, --help                Display this help and exit
 -v, --version             Display version info and exit
 -i, --input=<fasta file>  FASTA file with protein(s) (default: stdin)
 -m, --model=<liblinear model> (default:/usr/local/share/chipper/lr.model)
 -o, --output=<file>       File to write output (default: stdout)
 -n, --netchop             Output NetChop 3.0 format instead of FASTQ
 -s, --netchop_short       Output NetChop 3.0 short format
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
!!!!!!!!!++-++-+--++++++++----+++++----+--++++++++-+------++++++-++++-+-+-+---++++-+++++-+--+++++-------++++++++----+-+++++-++++++----+++++++++++--++-+++--++++++++-+-+++++--+++++++--+-++-+++++++++++-++++-++++-!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!-+-++++++++-+-+--+-----+-++++----++-----+-+-+-++-++-+-++++++++-----++++++----++++-+--+-++---+++++-+++++++++++++--------++++-+-+++--+--+++++++--+--++--+--++-++-+-++-++--+-+++-++++----+++-+-+-++++++--!!!!!!!!!!
```

Since `chipper` needs ten flanking amino acids in order to predict cleavage points, the first nine and last ten positions of each protein has '!`.

To see the default cutoff that `chipper` used with the logistic regression model, see the version output, e.g.

```
$ chipper -v
chipper version: 0.2.0
Best cutoff is 0.35 (MCC= 0.62)
Trained with 25655 publicly available MHC class I ligands
```

#### Predict cleavage positions using the default LR model and a custom cutoff.

For example, if you want to only mark a position as a cleavage position when the probability of cleavage is `>=` 80%, e.g.

```
$ chipper -i /usr/local/share/chipper/test.fa -c 0.8 | gunzip
@gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKKPAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKKK
+gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
!!!!!!!!!-+--+------++++++-----++++-------++-+--++-+------+---+--+-+-----------++--+--++----++--+---------+-++++-------+-+----+-------+++--+--------+-++----+-++-------+--+-----++-+--+------+-+--+--+-+----+--+-!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!------+----------------+-++-+----++-----+-+-+--+---------++++-------++-++------++-------+------++-----+--+++++---------+++--+--++-----+--++++------+-----++-+--+------------+--+++-----+----+-+++++---!!!!!!!!!!
```

When your compare the output to the previous (which used a cutoff of 0.35), you see less positions marked as cleavage positions. Using a higher probability reduces the number of false positives but also increases the number of false negatives.

#### Output probabilities of cleavage using the default Logistic Regression model

Instead of '-/+' notation, you can see the probabilities as each position, using the `-p/--probs` flag, e.g.

```
$ chipper -i /usr/local/share/chipper/test.fa **-p** | gunzip
@gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKKPAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKKK
+gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
!!!!!!!!!59158152047989999100179998000040199497599380200209767871959605050600178850843990500897681112000669499980100406979617597760123898649575540148189621596985651624974802773987912817613495967933839546195790!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!251356964770734126000219299780030992200191808049177070565998940000258949801004699153072480105568836447945899894002000018997081689104008469989204004811720980953917417401727690498900005972509389998411!!!!!!!!!!
```

#### Predict cleavage positions using support vector classification (SVC)

```
$ chipper -i /usr/local/share/chipper/test.fa -m /usr/local/share/chipper/svc.model  | gunzip
@gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRIKLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKKPAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKKK
+gi|121919|sp|P10412.2|H14_HUMAN RecName: Full=Histone H1.4; AltName: Full=Histone H1b; AltName: Full=Histone H1s-4
!!!!!!!!!++-++-+---+++++++----+++++----+--++-+++++-+------++++++-++++-+-+-+---++++-++-++-+--+++++-------++++++++----+-+++++-++++++----++++++++++---++-+++--++++++++-+-+++-+--++-++++--+-++--+++++++--+-++-+-++++-!!!!!!!!!!
@gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMTMASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTPK
+gi|138898|sp|P21431.1|NS1_I60A0 RecName: Full=Non-structural protein 1; Short=NS1; AltName: Full=NS1A
!!!!!!!!!-+--+++++++-+----+-----+-++++----++-----+-+-+--+-++-+-++++++++-----+++-++-----+++-+--+-++---+++++-+--++++++++++--------++++-+-+++-----+++++++-----++--+--++-++-+-++-+---+-+++-++++----+++-+-+-+++++---!!!!!!!!!!
```

The SVC model runs just as fast as the LR model but has slightly different sensitivity, specificity, and MCC metrics (see below).

In fact, you can generate your own `liblinear` model and pass it on the command-line. 

#### NetChop 3.0 formatted output

You can use chipper as a drop-in replacement for NetChop in your pipeline. Chipper will write results in both the NetChop 3.0 long and short format, e.g.

```
$ chipper -i /usr/local/share/chipper/test.fa -n
chipper 0.2.0 predictions using version C-term. Threshold 0.350000

--------------------------------------
 pos  AA  C      score      Ident
--------------------------------------
   1   M  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
   2   S  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
   3   E  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
   4   T  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
   5   A  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
   6   P  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
   7   A  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
   8   A  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
   9   P  ?    nan gi|121919|sp|P10412.2|H14_HUMAN
  10   A  S    0.573356 gi|121919|sp|P10412.2|H14_HUMAN
  11   A  S    0.925881 gi|121919|sp|P10412.2|H14_HUMAN
  12   P  .    0.181648 gi|121919|sp|P10412.2|H14_HUMAN
  13   A  S    0.534707 gi|121919|sp|P10412.2|H14_HUMAN
  14   P  S    0.873750 gi|121919|sp|P10412.2|H14_HUMAN
  15   A  .    0.129173 gi|121919|sp|P10412.2|H14_HUMAN
  16   E  S    0.550686 gi|121919|sp|P10412.2|H14_HUMAN
  17   K  .    0.218508 gi|121919|sp|P10412.2|H14_HUMAN
  18   T  .    0.075059 gi|121919|sp|P10412.2|H14_HUMAN
  19   P  S    0.419142 gi|121919|sp|P10412.2|H14_HUMAN
  20   V  S    0.789228 gi|121919|sp|P10412.2|H14_HUMAN
  21   K  S    0.969208 gi|121919|sp|P10412.2|H14_HUMAN
  22   K  S    0.889635 gi|121919|sp|P10412.2|H14_HUMAN
  23   K  S    0.981874 gi|121919|sp|P10412.2|H14_HUMAN
  24   A  S    0.962042 gi|121919|sp|P10412.2|H14_HUMAN
  25   R  S    0.963641 gi|121919|sp|P10412.2|H14_HUMAN
  26   K  S    0.969662 gi|121919|sp|P10412.2|H14_HUMAN
  27   S  .    0.135979 gi|121919|sp|P10412.2|H14_HUMAN
  28   A  .    0.037417 gi|121919|sp|P10412.2|H14_HUMAN
  29   G  .    0.023903 gi|121919|sp|P10412.2|H14_HUMAN
  30   A  .    0.188232 gi|121919|sp|P10412.2|H14_HUMAN
  31   A  S    0.749295 gi|121919|sp|P10412.2|H14_HUMAN
...
...
...
 204   R  S    0.868376 gi|138898|sp|P21431.1|NS1_I60A0
 205   S  S    0.439963 gi|138898|sp|P21431.1|NS1_I60A0
 206   S  .    0.121788 gi|138898|sp|P21431.1|NS1_I60A0
 207   D  .    0.119843 gi|138898|sp|P21431.1|NS1_I60A0
 208   E  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 209   N  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 210   G  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 211   R  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 212   P  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 213   P  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 214   L  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 215   T  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 216   P  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
 217   K  ?    nan gi|138898|sp|P21431.1|NS1_I60A0
--------------------------------------

Number of cleavage sites 115. Number of amino acids 217. Protein name gi|138898|sp|P21431.1|NS1_I60A0

--------------------------------------
```

You can also output the NetChop short format, e.g.

```
$ chipper -i /usr/local/share/chipper/test.fa -s
chipper 0.2.0 predictions using version C-term. Threshold 0.350000

219 gi|121919|sp|P10412.2|H14_HUMAN
MSETAPAAPAAPAPAEKTPVKKKARKSAGAAKRKASGPPVSELITKAVAASKERSGVSLAALKKALAAAGYDVEKNNSRI
?????????SS.SS.S..SSSSSSSS....SSSSS....S..SSSSSSSS.S......SSSSSS.SSSS.S.S.S...SS
KLGLKSLVSKGTLVQTKGTGASGSFKLNKKAASGEAKPKAKKAGAAKAKKPAGAAKKPKKATGAATPKKSAKKTPKKAKK
SS.SSSSS.S..SSSSS.......SSSSSSSS....S.SSSSS.SSSSSS....SSSSSSSSSSS..SS.SSS..SSSSS
PAAAAGAKKAKSPKKAKAAKPKKAPKSPAKAKAVKPKAAKPKTAKPKAAKPKKAAAKK
SSS.S.SSSSS..SSSSSSS..S.SS.SSSSSSSSSSS.SSSS.SSSS.?????????
--------------------------------------

Number of cleavage sites 135. Number of amino acids 219. Protein name gi|121919|sp|P10412.2|H14_HUMAN

--------------------------------------
217 gi|138898|sp|P21431.1|NS1_I60A0
MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQKSLRGRGSTLGLNIETATRVGKQIVERILKEESDEALKMT
?????????.S.SSSSSSSS.S.S..S.....S.SSSS....SS.....S.S.S.SS.SS.S.SSSSSSSS.....SSSS
MASAPASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCIRMDQAIMDKNIILKANFSVIFDRLETLILLRAFTEAGAIVGEI
SS....SSSS.S..S.SS...SSSSS.SSSSSSSSSSSSS........SSSS.S.SSS..S..SSSSSSS..S..SS..S
SPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSKTLQRFAWRSSDENGRPPLTP
..SS.SS.S.SS.SS..S.SSS.SSSS....SSS.S.S.SSSSSS..?????????
--------------------------------------

Number of cleavage sites 115. Number of amino acids 217. Protein name gi|138898|sp|P21431.1|NS1_I60A0

--------------------------------------
```

Since chipper uses properties of the flanking amino acids, it does not make predictions at the ends of the peptides in the positions marked with '?'. This will change in future versions of chipper.

## Performance Comparison

Chipper compares favorably against other proteasomal cleavage detection systems. Metrics reported in "[Predicting proteasomal cleavage sites: a comparison of available methods](http://intimm.oxfordjournals.org/content/15/7/781.full)" compared with chipper. Here's the performance using the same dataset using the same test dataset.


| Method | Sensitivity | Specificity | MCC    |
|--------|-------------|-------------|--------|
|[PAProC](http://www.paproc.de/)| 45.6 | 30.0 | -0.25 |
|[FragPredict](http://www.mpiib-berlin.mpg.de/MAPPP/fragpredict.html) | 83.5 | 16.5 | 0.00 |
|[NetChop](http://www.cbs.dtu.dk/services/NetChop/) 1.0 | 39.8 | 46.3 | -0.14 |
|[NetChop](http://www.cbs.dtu.dk/services/NetChop/) 2.0 | 73.6 | 42.4 | 0.16 |
|[NetChop](http://www.cbs.dtu.dk/services/NetChop/) 3.0 | 81.0 | 48.0 | 0.31 |
|chipper (LR) | **87.0** | **74.5** | **0.620** |
|chipper (SVC) | **79.3** | **77.9** | **0.572** |

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



