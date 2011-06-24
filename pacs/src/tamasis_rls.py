#!/usr/bin/env python
"""
Command-line interface to perform PACS map-making with a simple PACS
instrument model.

Author: Nicolas Barbey
"""


options = "hvf:"
long_options = ["help", "verbose", "filenames="]

def main():
    import os, sys, getopt, ConfigParser, time
    from tamasis.pacs.processing import tamasis_mapper_rls
    # parse command line arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], options, long_options)
    except getopt.GetoptError, err:
        print(str(err))
        usage()
        sys.exit(2)

    # default
    verbose = False

    # parse options
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-v", "--verbose"):
            verbose = True

    # read config file
    if len(args) == 0:
        print("Error: config filename argument is mandatory.\n")
        usage()
        sys.exit(2)
    config_file = args[0]
    config = ConfigParser.RawConfigParser()
    config.read(config_file)
    keywords = dict()
    sections = config.sections()
    sections.remove("main")
    # parse config
    for section in sections:
        keywords[section] = dict()
        for option in config.options(section):
            get = config.get
            # recast to bool, int or float if needed
            if section == "PacsObservation":
                if option == "reject_bad_line":
                    get = config.getboolean
                if option == "fine_sampling_factor":
                    get = config.getint
                if option in ("active_fraction",
                              "delay",
                              "calblock_extension_time"):
                    get = config.getfloat
            if section == "get_tod":
                if option in ("flatfielding",
                              "substraction_mean",
                              "raw"):
                    get = config.getboolean
            if section == "Projection":
                if option in ("oversampling",
                              "packed"):
                    get = config.getboolean
                if option == "npixels_per_sample":
                    get = config.getint
                if option == "resolution":
                    get = config.getfloat
            if section == "mapper_rls":
                if option in ("verbose",
                              "criterion"):
                    get = config.getboolean
                if option == "maxiter":
                    get = config.getint
                if option in ("hyper",
                              "tol"):
                    get = config.getfloat
            # store option using the appropriate get to recast option.
            keywords[section][option] = get(section, option)
    # special case for the main section
    data_file_list = config.get("main", "filenames").split(",")
    # if filenames argument is passed, override config file value.
    for o, a in opts:
        if o in ("-f", "--filenames"):
            data_file_list = a.split(", ")
    data_file_list = [w.lstrip().rstrip() for w in data_file_list]
    # append date string to the output file to distinguish results.
    date = time.strftime("%y%m%d_%H%M%S", time.gmtime())
    filename = data_file_list[0].split(os.sep)[-1]
    # remove extension
    fname = ".".join(filename.split(".")[:-1])
    # store results into the Data subdirectory as expected by sumatra
    output_file = "Data/map" + fname + '_' + date + '.fits'
    # run tamasis mapper
    tamasis_mapper_rls(data_file_list, output_file, keywords, verbose=verbose)

def usage():
    print(__usage__)

__usage__ = """Usage: tamasis_rls [options] [config_file]

Use tamasis regularized least-square map-making routine
using the configuration [filename]

[config_file] is the name of the configuration file. This file
contains arguments to the various processing steps of the
map-maker. For an exemple, see the file tamasis_rls.cfg in the
tamasis/pacs/src/ directory.

Options:
  -h, --help        Show this help message and exit
  -v, --verbose     Print status messages to std output.
  -f, --filenames   Overrides filenames config file value.
"""

# to call from command line
if __name__ == "__main__":
    main()
