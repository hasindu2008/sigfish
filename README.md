# sigfish


sigfish is an experiment toolkit that attempts to directly map nanopore raw signal data to a reference. Supports R9 DNA and RNA, R10 DNA and the latest RNA004 kits from ONT.

**This is under construction. Interface and parameters are not stable. Documentation is currently minimal**

This repository of sigfish supports CPU only. For the FPGA accelerated version, please refer to the [sigfish-haru fork](https://github.com/beebdev/sigfish-haru) for which the documentation is located [HARU](https://github.com/beebdev/HARU) and the publication at [GIGAScience](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giad046/7217084).


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

Currently, there are two subtools: *dtw* and *eval*.

## dtw

Performs subsequence Dynamic Time Warping (sDTW) of raw signals in [S/BLOW5 format](https://www.nature.com/articles/s41587-021-01147-4) to a reference in FASTA format. This is an all-to-all alignment and is not intended for large references.

```
Usage: sigfish dtw [OPTIONS] genome.fa/transcriptome.fa reads.blow5
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


If you specify `--sam`, the output will be in SAM format containing ss and si tags similar to the format described [here](https://hasindu2008.github.io/f5c/docs/output#eventalign-sam-output). This output can be used to visualise the alignment using [squigualiser](https://github.com/hasindu2008/squigualiser).

Options:
```
basic options:
   -t INT                     number of processing threads [8]
   -K INT                     batch size (max number of reads loaded at once) [512]
   -B FLOAT[K/M/G]            max number of bytes loaded at once [20.0M]
   -h                         help
   -o FILE                    output to file [stdout]
   --verbose INT              verbosity level [4]
   --version                  print version
   --pore STR                 set the pore chemistry (r9, r10 or rna004) [auto]
advanced options:
   --kmer-model FILE          custom nucleotide k-mer model file (format similar to test/r9-models/r9.4_450bps.nucleotide.6mer.template.model)
   --rna                      the dataset is direct RNA
   -q INT                     the number of events in query signal to align [250]
   -p INT                     the number of events to trim at query signal start [50]
   --debug-break INT          break after processing the specified no. of batches
   --profile-cpu=yes|no       process section by section (used for profiling on CPU)
   --dtw-std                  use DTW standard instead of DTW subsequence
   --invert                   reverse the reference events instead of query
   --full-ref                 map to the full reference
   --from-end                 Map the end portion of the query instead of the beginning
   --sam                      Output in SAM format
```


## eval

Evaluates/compare mappings in PAF format (testset) by comparing to a truthset which is also in PAF format.  As an example, testset is the output from *sigfish dtw*, whereas, truthset can be the mapping output from [Minimap2](https://github.com/lh3/minimap2). The mapping position (includes contig name, coordinates and strandness) for read ID in the testset is compared agianst those in the truthset. A mapping is considered correct if the contig name and the strandness in the testset exactly matches the truthset and the mapping coodinates match the criteria `min(|diff_st|,|diff_end|) < THRESHOLD` where,

```
truth:   -----------------------------
test:                    -----------------------------------
         |<---diff_st--->|            |<-----diff_end----->|
```

THRESHOLD is 100 by default.

Output includes statistics for mapping accuracy considering the testset as a whole, as well as based on individual mapping qualities scores.

Usage:

```
sigfish eval truth.paf test.paf
```

Options:
```
basic options:
   -h                         help
   --version                  print version
   --secondary STR            consider secondary mappings. yes or no.
   --tid-only                 consider regerence name and strand only
```

## Acknowledgement

- The output PAF-like format output by *sigfish dtw* was inspired by [UNCALLED](https://github.com/skovaka/UNCALLED). *sigfish eval* was implemented by learning from [UNCALLED](https://github.com/skovaka/UNCALLED) pafstats.
- The methodology in [ReadUntil](https://github.com/mattloose/RUscripts) was referred to when implementing alignment component in *sigfish dtw*
- The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
- The DTW code is from [mlpy](http://mlpy.sourceforge.net/).
- The pore-models are from [Nanopolish](https://github.com/jts/nanopolish).
- Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2), [Samtools](http://samtools.sourceforge.net/).

# Citation

You may cite the following:

> Po Jui Shih, Hassaan Saadat, Sri Parameswaran, Hasindu Gamaarachchi, Efficient real-time selective genome sequencing on resource-constrained devices, GigaScience, Volume 12, 2023, giad046, https://doi.org/10.1093/gigascience/giad046

```
@article{shih2023efficient,
  title={Efficient real-time selective genome sequencing on resource-constrained devices},
  author={Shih, Po Jui and Saadat, Hassaan and Parameswaran, Sri and Gamaarachchi, Hasindu},
  journal={GigaScience},
  volume={12},
  pages={giad046},
  year={2023},
  publisher={Oxford University Press}
}
```