#!/bin/sh
# Generate Arguments parser from the onmi.ggo

sh cmdline &> /dev/null
gengetopt  --unamed-opts --output-dir cmdline < onmi.ggo

echo  The arguments parser is generated
