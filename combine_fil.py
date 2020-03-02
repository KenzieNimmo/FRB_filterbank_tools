"""
To add two filterbank files together (such as lcp and rcp). Filterbank files must have the same number of time samples and frequency number_of_samples
Kenzie Nimmo 2020
"""
import filterbank
import numpy as np


def write_header(outfbfile, header_params, header):
    """Given list of input filterbank objects and an output
        file object write the new filterbank header.
    """
    for hp in header_params:
        hdrval = sigproc.addto_hdr(hp, header[hp])
        outfbfile.write(hdrval)

def write_data(infbfiles,outfbfile):
    """Given input filterbank, a start sample and number of samples to extract
        outputs the chopped filterbank file.
    """
    for i in range(0,testfil[0].number_of_samples):
        data=np.zeros(testfil[0].header['nchans'])
        for j in range(len(testfil)):
                data+=testfil[j].read_sample()
                data.tofile(outfbfile)

def main():
    infbfiles = filterbank.filterbank(options.infile)
    #open output fb file
    outfbfile = open(options.outname, 'wb')

    #write header
    header = infbfiles[0].header
    header_params = infbfiles[0].header_params
    write_header(outfbfile,header_params,header)

    #write data
    write_data(infbfiles, outfbfile)
    outfbfile.close()
    infbfiles.close()

if __name__=='__main__':
    parser = optparse.OptionParser(usage='%prog [options] infile', \
                description="To extract N samples from a filterbank file")
    parser.add_option('-o', '--outname', dest='outname', type='string', \
                help="Output filename.", default=None)
    (options, args) = parser.parse_args()

    if len(args)==0:
        parser.print_help()
        sys.exit(1)
    elif len(args)==1:
        sys.stderr.write("More than one input file must be provided!\n")
    else:
        options.infile = args

    if options.outname is None:
        sys.stderr.write("An outname must be provided. " \
                            "(Use -o/--outname on command line).\n")
        sys.exit(1)

    main()
