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

# Usage

```
$ onmi FILE1 FILE2
```
Execution Options:
```
-h, --help     Print help and exit
-V, --version  Print version and exit
-v, --verbose  detailed debugging  (default=off)
-o, --omega    print the Omega measure (can be slow)  (default=off)
```
The filenames record the sets of communities. A typical use case is to have
the "true" communities in one file and and those found by your algorithm
in the other file.

Input file format is CNL (cluster nodes list). One line per community. The nodes are
separated by whitespace, and any non-whitespace characters may be used in the
node names. Line comments are allowed using `#`. Example of the CNL format:
```
# The comments start with '#' like this line
# Each non-commented line is a module(cluster, community) consisting of the the member nodes separated by space / tab
1
1 2
2
```

# Related Projects
- [GenConvNMI](https://github.com/eXascaleInfolab/GenConvNMI) - Overlapping NMI evaluation that is compatible with the original NMI (unlike the `onmi`).
- [ExecTime](https://bitbucket.org/lumais/exectime/)  - A lightweight resource consumption profiler.
- [PyCABeM](https://github.com/eXascaleInfolab/PyCABeM) - Python Benchmarking Framework for the Clustering Algorithms Evaluation. Uses extrinsic (NMIs) and intrinsic (Q) measures for the clusters quality evaluation considering overlaps (nodes membership by multiple clusters).
