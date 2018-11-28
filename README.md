# OvpNMI - Overlapping NMI and Omega Index evaluation

An implementation of a Normalized Mutual Information (NMI) measure for sets of overlapping clusters and Omega Index.  

*ATTENTION:*
1. OvpNMI (onmi) does not suitable to evaluate multi-resolution or highly overlapping clusterings because of the used "best match" approximation (see formula 3 in the paper). To evaluate such cases [GenConvNMI](https://github.com/eXascaleInfolab/GenConvNMI) can be applied, which works slower, but suitable for all cases.
2. Omega Index evaluation is pretty slow here, use [xmeasures](https://github.com/eXascaleInfolab/xmeasures) for the order of magnitude faster evaluation.

*Internal note:* the counter example for this measure (onmi) is `rotshift_single.cnl`.

The paper: *["Normalized Mutual Information to evaluate overlapping community finding algorithms"](http://arxiv.org/abs/1110.2515) by Aaron F. McDaid, Derek Greene, Neil Hurley*  
This method is based on the method described in Appendix B at the end of:
*["Detecting the overlapping and hierarchical community structure in complex networks"](http://iopscience.iop.org/1367-2630/11/3/033015/)
by Andrea Lancichinetti, Santo Fortunato and János Kertész*

Author: Aaron F. McDaid <aaronmcdaid@gmail.com>  

This is a fork of the original [onmi](https://github.com/aaronmcdaid/Overlapping-NMI)
with the extension purposes to be used in the [PyCaBeM](https://github.com/eXascaleInfolab/PyCABeM) clustering benchmark. The extendsion and bugfixes include: modification of the I/O, node base synchronization, NMI_sqrt added, normalization border cases fixed (the case of fully overlapping multiple clusters), etc.  
Changes made by Artem Lutov <artem@exascale.info>

## Content
- [Deployment](#deployment)
	- [Dependencies](#dependencies)
	- [Compilation](#compilation)
- [Usage](#usage)
- [Related Projects](#related-projects)

# Deployment

## Dependencies
There no any dependencies for the execution or compilation.  
However, to extend the input options and automatically regenerate the input parsing,
[*gengetopt*](https://www.gnu.org/software/gengetopt) application should be installed: `$ sudo apt-get install gengetopt`.

## Compilation
Just execute `$ make`.  
To update/extend the input parameters just modify `args.ggo` and run `GenerateArgparser.sh` (calls `gengetopt`).

# Usage
```
$ onmi clsfile1 clsfile2
```
_Applicability Note:_ OvpNMI is extremely fast, but does not suitable to evaluate multi-resolution clusterings, see [GenConvNMI](https://github.com/eXascaleInfolab/GenConvNMI) instead.

Execution Options:
```
$ ./onmi -h
onmi 0.3

Compare sets of clusters by their members (nodes) using various measures (NMI,
Omega) and considering overlaps

Usage: onmi [OPTIONS] clsfile1 clsfile2

  -h, --help              Print help and exit
  -V, --version           Print version and exit
  -s, --sync=filename     synchronize the node base omitting the non-matching
                            nodes.
                            NOTE: The node base is either the first input file
                            or '-' (automatic selection of the input file
                            having the least number of nodes).
  -a, --allnmis           output all NMIs (sqrt and sum-denominators, LFK
                            besides the max-denominator)  (default=off)
  -m, --membership=FLOAT  average expected membership of nodes in the clusters,
                            > 0, typically >= 1  (default=`1')
  -o, --omega             print the Omega measure (can be slow)  (default=off)
  -t, --textid            use text ids of nodes instead of .cnl format
                            (default=off)
  -v, --verbose           detailed debugging  (default=off)
```
The input files contain list of clusters (communities, modules). A typical use case is to have
the "true" communities in one file and and those found by your algorithm in the other file.

The default input file format is CNL (cluster nodes list), where each cluster is
represented by one line. The nodes are separated by whitespace, and any non-whitespace
characters may be used in the node names. Empty lines and comments (lines starting with #) are skipped.
Example of the CNL format:
```
# The comments start with '#' like this line
# Each non-commented line is a module(cluster, community) consisting of the the member nodes separated by space / tab
1
1 2
2
```
A node id is unsigned integer by default, and it can be any word not starting with the comment symbol `#` if `-t` option is specified to use text ids.
>
- Any line starting with `#` is omitted as a comment, also as any remained part of the line starting with `#` in the *textid* mode
- Ids can't contain `:` symbol, because it is used to specify the membership share in the CNL format, which is not supported by onmi. The id part starting from the `:` symbol is omitted (trimmed).

**Note:** Please, [star this project](https://github.com/eXascaleInfolab/OvpNMI) if you use it.

# Related Projects
- [xmeasures](https://github.com/eXascaleInfolab/xmeasures)  - Extrinsic quality (accuracy) measures evaluation for the overlapping clustering on large datasets: family of mean F1-Score (including clusters labeling), Omega Index (fuzzy version of the Adjusted Rand Index) and standard NMI (for non-overlapping clusters).
- [GenConvNMI](https://github.com/eXascaleInfolab/GenConvNMI) - Overlapping NMI evaluation that is (unlike `onmi`) compatible with the original NMI and suitable for both overlapping and multi resolution (hierarchical) clusterings.
- [Clubmark](https://github.com/eXascaleInfolab/clubmark) - A parallel isolation framework for benchmarking and profiling clustering (community detection) algorithms considering overlaps (covers).
- [resmerge](https://github.com/eXascaleInfolab/resmerge)  - Resolution levels clustering merger with filtering. Flattens hierarchy/list of multiple resolutions levels (clusterings) into the single flat clustering with clusters on various resolution levels synchronizing the node base.
- [ExecTime](https://bitbucket.org/lumais/exectime/)  - A lightweight resource consumption profiler.
- [TInfES](https://github.com/eXascaleInfolab/TInfES)  - Type inference evaluation scripts and accessory apps used for the benchmarking.
