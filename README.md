# sigfish

**This is under construction. Interface and parameters are not stable.**

sigfish is an experiment toolkit that attempts to directly map nanopre raw signal data to a reference.

## Building

```
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone https://github.com/hasindu2008/sigfish
cd sigfish
make
```

The commands to install zlib __development libraries__ on some popular distributions:
```sh
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```

## Usage

Current, there are two subtools: dtw and eval.

## dtw

Performs subsequence Dynamic Time Warping (sDTW) of raw signals in [S/BLOW5 format](https://www.nature.com/articles/s41587-021-01147-4) to a reference in FASTA format. This is an all-to-all alignment and is not intended for large references.

```
Usage: sigfish dtw [OPTIONS] genome.fa reads.blow5
```

Output is in a PAF-like format with the following columns:

|Col|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|1  |string|Read identifier name                       |
|2  |int   |Raw signal length (number of samples)                    |
|3  |int   |Raw signal start index  (0-based; BED-like; closed)   |
|4  |int   |Raw signal end index (0-based; BED-like; open)       |
|5  |char  |Relative strand: "+" or "-"               |
|6  |string|Reference name                     |
|7  |int   |Reference sequence length                    |
|8  |int   |start on the reference sequence (0-based; BED-like; closed)  |
|9  |int   |end on reference sequence (0-based; BED-like; open)   |
|10 |int   |Approximation of number of matching bases in the alignment                |
|11 |int   |Alignment block length on the reference in terms of bases                 |
|12 |int   |Mapping quality (0-255; 255 for missing)  |

Following optional tags are present:

|Tag|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|tp  |A   |Type of alignment: P/primary                      |
|d1  |f   |DTW score (lower the better)                     |
|d2  |f   |DTW score of the next best alignment                       |


## eval

Evaluates/compare mappings in PAF format (testset) by comparing to a truthset which is also in PAF format.  As an example, testset is the output from *sigfish dtw*, whereas, truthset can be the mapping output from [Minimap2](https://github.com/lh3/minimap2). The mapping coordinate (and the strand) for read ID in the testset is compared agianst those in the truthset. Output includes statistics for mapping accuracy considering the testset as a whole, as well as teh accuracy based on individual mapping qualities scores.

Usage:

```
sigfish eval truth.paf test.paf
```


## Acknowledgement

- The output PAF-like format output by *sigfish dtw* was inspired by [UNCALLED](https://github.com/skovaka/UNCALLED). *sigfish eval* was implemented by learning from [UNCALLED](https://github.com/skovaka/UNCALLED) pafstats.
- The methodology in [ReadUntil](https://github.com/mattloose/RUscripts) was referred to when implementing alignment component in *sigfish dtw*
- The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
- The DTW code is from [mlpy](http://mlpy.sourceforge.net/).
- The pore-models are from [Nanopolish](https://github.com/jts/nanopolish).
- Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2), [Samtools](http://samtools.sourceforge.net/).

