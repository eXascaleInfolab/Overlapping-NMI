#!/bin/sh
# Generate Arguments parser from the onmi.ggo

mkdir cmdline
gengetopt  --unamed-opts --output-dir cmdline < onmi.ggo
