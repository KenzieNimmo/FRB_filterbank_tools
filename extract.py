"""
To chop out part of a filterbank file. From begin_sample to begin_sample+duration.
Kenzie Nimmo 2020
"""

import sys
import numpy as np
import filterbank
import sigproc
import optparse


def write_header(outfbfile, header_params, header):
    """Given list of input filterbank objects and an output
        file object write the new filterbank header.
    """
    for hp in header_params:
        hdrval = sigproc.addto_hdr(hp, header[hp])
        outfbfile.write(hdrval)

def write_data(infbfile,begin_N,dur_N,outfbfile):
    """Given input filterbank, a start sample and number of samples to extract
        outputs the chopped filterbank file.
    """
    infbfile.seek_to_sample(begin_N)
    for i in range(begin_N,(begin_N+dur_N)):
        data = infbfile.read_sample()
        data.tofile(outfbfile)


def main():
    infbfile = filterbank.filterbank(options.infile)
    #open output fb file
    outfbfile = open(options.outname, 'wb')

    begin_N = options.begin_sample
    dur_N = options.duration

    #write header
    header = infbfile.header
    header_params = infbfile.header_params
    write_header(outfbfile,header_params,header)

    #write data
    write_data(infbfile, begin_N, dur_N, outfbfile)
    outfbfile.close()
    infbfile.close()


if __name__=='__main__':
    parser = optparse.OptionParser(usage='%prog [options] infile', \
                description="To extract N samples from a filterbank file")
    parser.add_option('-o', '--outname', dest='outname', type='string', \
                help="Output filename.", default=None)
    parser.add_option('-b', '--begin_sample', dest='begin_sample', type='int', \
                help="Begin sample number to start extract", default=0)
    parser.add_option('-d', '--duration', dest='duration', type='int', \
                help="number of bins to extract", default=None)
    (options, args) = parser.parse_args()

    if len(args)==0:
        parser.print_help()
        sys.exit(1)
    elif len(args)!=1:
        sys.stderr.write("Only one input file must be provided!\n")
    else:
        options.infile = args[-1]

    if options.outname is None:
        sys.stderr.write("An outname must be provided. " \
                            "(Use -o/--outname on command line).\n")
        sys.exit(1)

    if options.duration is None:
        sys.stderr.write("Number of bins to extract must be provided. " \
                            "(Use -d/--duration on command line).\n")
        sys.exit(1)

    main()
