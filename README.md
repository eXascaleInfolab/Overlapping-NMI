# OvpNMI - Overlapping NMI

An implementation of a Normalized Mutual Information (NMI) measure for sets of overlapping clusters.

The paper: *["Normalized Mutual Information to evaluate overlapping community finding algorithms"](http://arxiv.org/abs/1110.2515) by Aaron F. McDaid, Derek Greene, Neil Hurley*  
This method is based on the method described in Appendix B at the end of:
*["Detecting the overlapping and hierarchical community structure in complex networks"](http://iopscience.iop.org/1367-2630/11/3/033015/)
by Andrea Lancichinetti, Santo Fortunato and János Kertész*

Author: Aaron F. McDaid <aaronmcdaid@gmail.com>  

This is a fork of the original [onmi](https://github.com/aaronmcdaid/Overlapping-NMI)
with the extension purposes (mainly modification of the I/O) to be used in the [PyCaBeM](https://github.com/eXascaleInfolab/PyCABeM) clustering benchmark.  
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
To update/extend the input parameters just modify `onmi.ggo` and run `GenerateArgparser.sh` (calls `gengetopt`).

# Usage

```
$ onmi clsfile1 clsfile2
```
Execution Options:
```
$ ./onmi -h
onmi 0.2

Compare sets of clusters by their members (nodes) using various measures (NMI,
Omega) and considering overlaps

Usage: onmi [OPTIONS] clsfile1 clsfile2

-h, --help     Print help and exit
-V, --version  Print version and exit
-v, --verbose  detailed debugging  (default=off)
-a, --allnmis  outmut all NMIs (sum-denominator and LFK besides the
                 max-denominator one)  (default=off)
-o, --omega    print the Omega measure (can be slow)  (default=off)
-t, --textid   use text ids of nodes instead of .cnl format  (default=off)
```
The input files contain list of clusters (communities, modules). A typical use case is to have
the "true" communities in one file and and those found by your algorithm in the other file.

The default input file format is CNL (cluster nodes list), where each cluster is
represented by one line. The nodes are separated by whitespace, and any non-whitespace
characters may be used in the node names. Line comments are allowed using `#`.  
Example of the CNL format:
```
# The comments start with '#' like this line
# Each non-commented line is a module(cluster, community) consisting of the the member nodes separated by space / tab
1
1 2
2
```
A node id is unsigned integer by default, and it can be any word not starting with the comment symbol `#` if `-t` option is specified to use text ids.

> Any line starting with `#` is omitted as a comment, also as any remained part of the line starting with `#` in the *textid* mode

# Related Projects
- [GenConvNMI](https://github.com/eXascaleInfolab/GenConvNMI) - Overlapping NMI evaluation that is compatible with the original NMI (unlike the `onmi`).
- [ExecTime](https://bitbucket.org/lumais/exectime/)  - A lightweight resource consumption profiler.
- [PyCABeM](https://github.com/eXascaleInfolab/PyCABeM) - Python Benchmarking Framework for the Clustering Algorithms Evaluation. Uses extrinsic (NMIs) and intrinsic (Q) measures for the clusters quality evaluation considering overlaps (nodes membership by multiple clusters).
