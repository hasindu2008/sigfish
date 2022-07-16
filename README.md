# sigfish

## Building

Dependencies: zlib development libraries

```
git clone https://github.com/hasindu2008/sigfish
cd sigfish
make
```

## Main tools

### dtw

Performs dynamic time warping-based alignment of raw signals in S/BLOW5 format to a reference in FASTA format. This is an all-to-all alignment and is not intended for large references.

```
Usage: sigfish dtw [OPTIONS] genome.fa reads.blow5
```

### eval

Similar to UNCALLED pafstats.


## Hepler tools

Following subtools have two possible use cases.

```
# To perform subtool operation on all reads in a BLOW5 file
sigfish <subtool> reads.blow5
# To perform subtool operation on specified read IDs in a BLOW5 file
sigfish <subtool> reads.blow5 read_id1 read_id2 ..
```

By default, a tab-delimited text file with the first row being the header is printed. You can suppress the header using `-n` flag.


## event

Event segmentation based on the method in Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).


### long output

|Col|Type  |Name            |Description                                                        |
|--:|:----:|:------:        |:-----------------------------------------                         |
|1  |string|read_id         |Read identifier name                                               |
|2  |int   |event_idx       |Event index (0-based)                                              |
|3  |int   |raw_start       |Raw signal start index for the event (0-based; BED-like; closed)   |
|4  |int   |raw_end         |Raw signal end index for the event (0-based; BED-like; open)       |
|5  |float |event_mean      |Mean level of pico-ampere scaled signal for the event              |
|6  |float |event_std       |Standard deviations of pico-ampere scaled signal for the event     |

### compact output

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |int   |raw_start       |Raw signal start index of the first event (0-based; BED-like; closed)  |
|4  |int   |raw_end         |Raw signal end index of the last event (0-based; BED-like; open)       |
|5  |int*  |events          |Command separated event lengths (based on no. raw signals samples)     |

The event 0 starts at raw signal index `raw_start` (0-based; BED-like; closed) and ends at `raw_start+event[0]` (0-based; BED-like; open).
The event 1 starts at raw signal index `raw_start+event[0]` (0-based; BED-like; closed) and ends at `raw_start+event[0]+event[1]` (0-based; BED-like; open).
Likewise, the events can be reconstructed by using the cumulative sum of `events`.


## stat

Prints signal statistics.

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |float |raw_mean        |Mean of raw signal values                                              |
|4  |float |pa_mean         |Mean of pico-amperes scaled signal                                     |
|5  |float |raw_std         |Standard deviation of raw signal values                                |
|6  |float |pa_std          |Standard deviation of pico-amperes scaled signal                       |
|7  |int   |raw_median      |Median of raw signal values                                            |
|8  |float |pa_median       |Mean of pico-amperes scaled signal                                     |

## sref

Prints synthetic reference signal

## seg

## jnn


## pa

Prints the raw signal in pico-amperes.

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |float*|pa              |Comma separated Raw signal in pico amperes                             |


## Acknowledgement
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
The DTW code is from [mlpy](http://mlpy.sourceforge.net/).
The pore-models are from [Nanopolish](https://github.com/jts/nanopolish).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2), [Samtools](http://samtools.sourceforge.net/).
