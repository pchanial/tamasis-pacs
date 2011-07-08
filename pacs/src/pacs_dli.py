#!/usr/bin/env python

from __future__ import division

options = "hvf:o:"
long_options = ["help", "verbose", "filenames=", "output="]

def main():
    import os, sys, getopt, ConfigParser, time
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
        print("Error: config filename is mandatory.\n")
        usage()
        sys.exit(2)
    config_file = args[0]
    config = ConfigParser.RawConfigParser()
    config.read(config_file)
    keywords = dict()
    # parse config
    for section in config.sections():
        keywords[section] = dict()
        for option in config.options(section):
            get = config.get
            # recast to bool, int or float if needed
            # if option not handled here, it defaults to a string
            if section == "PacsObservation":
                if option == "reject_bad_line":
                    get = config.getboolean
                if option == "fine_sampling_factor":
                    get = config.getint
                if option in ("active_fraction",
                              "delay",
                              "calblock_extension_time"):
                    get = config.getfloat
            if section == "scanline_masking":
                if option in ("n_repetition", "n_scanline"):
                    get = config.getint
            if section == "get_tod":
                if option in ("flatfielding",
                              "substraction_mean",
                              "raw"):
                    get = config.getboolean
            if section == "deglitching":
                if option in ("length", "nsigma"):
                    get = config.getfloat
            if section == "filter_median":
                if option == "length":
                    get = config.getfloat
            if section == "Projection":
                if option in ("oversampling",
                              "packed"):
                    get = config.getboolean
                if option == "npixels_per_sample":
                    get = config.getint
                if option == "resolution":
                    get = config.getfloat
            if section == "map_mask":
                if option == "threshold":
                    get = config.getfloat
            if section == "dli":
                if option == "tau":
                    get = config.getfloat
            if section == "fmin_args":
                if option == "maxiter":
                    get = config.getfloat
                if option == "disp":
                    get = config.getboolean
            if section == "lanczos":
                if option == "maxiter":
                    get = config.getfloat
            # store option using the appropriate get to recast option.
            keywords[section][option] = get(section, option)
    # special case for the main section
    data_file_list = config.get("main", "filenames").split(", ")
    # if filenames argument is passed, override config file value.
    for o, a in opts:
        if o in ("-f", "--filenames"):
            data_file_list = a.split(", ")
    data_file_list = [w.rstrip().lstrip() for w in data_file_list]
    # append date string to the output file to distinguish results.
    date = time.strftime("%y%m%d_%H%M%S", time.gmtime())
    # extract filename from data_file
    filename = data_file_list[0].split(os.sep)[-1]
    # remove extension
    fname = ".".join(filename.split(".")[:-1])
    # store results into the Data subdirectory as expected by sumatra
    output_file = "Data/map" + fname + '_' + date + '.fits'
    # if output argument is passed, override config file value.
    for o, a in opts:
        if o in ("-o", "--output"):
            output_file = a
    # check existence of output path
    outdir = os.path.dirname(output_file)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    # run tamasis mapper
    pipeline_dli(data_file_list, output_file, keywords, verbose=verbose)

def pipeline_dli(filename, output_file, keywords, verbose=False):
    import numpy as np
    import linear_operators as lo
    import tamasis as tm
    # verbosity
    tm.var.verbose = verbose
    # define observation
    obs = tm.PacsObservation(filename, **keywords["PacsObservation"])
    # extra masking
    tm.step_scanline_masking(obs, **keywords["scanline_masking"])
    # get data
    tod = obs.get_tod(**keywords["get_tod"])
    # deglitching
    tm.step_deglitching(obs, tod, **keywords["deglitching"])
    # median filtering
    tod = tm.filter_median(tod, **keywords["filter_median"])
    # define projector
    projection = tm.Projection(obs, **keywords["Projection"])
    # build instrument model
    response = tm.ResponseTruncatedExponential(obs.pack(
            obs.instrument.detector.time_constant) / obs.SAMPLING_PERIOD)
    compression = tm.CompressionAverage(obs.slice.compression_factor)
    masking = tm.Masking(tod.mask)
    model = masking * compression * response * projection
    # set tod masked values to zero
    tod = masking(tod)
    # N^-1 operator
    invntt = tm.InvNtt(obs)

    # for the dli algorithm
    M = map_mask(tod, model, **keywords["map_mask"])

    # recast model as a new-style LinearOperator from linear_operators package
    H = lo.aslinearoperator(model) * M.T
    N = lo.aslinearoperator(invntt)
    # vectorize data so it is accepted by LinearOperators
    y = tod.ravel()
    Ds = [tm.DiscreteDifference(axis=i, shapein=projection.shapein) for i in (0, 1)]
    Ds = [lo.aslinearoperator(D) for D in Ds]
    Ds = [D * M.T for D in Ds]
    D = lo.concatenate(Ds)
    # handle tau which needs to be an ndarray
    keywords["dli"]["tau"] *= np.ones(D.shape[0])
    algo = lo.DoubleLoopAlgorithm(H, y, D, noise_covariance=N,
                                  fmin_args=keywords["fmin_args"],
                                  lanczos=keywords["lanczos"],
                                  **keywords["dli"])
    # optimize
    xe = algo()
    # reshape
    xe = (M.T * xe).reshape(projection.shapein)
    # recast as tamasis map
    xe = tm.Map(xe)
    # save
    xe.save(output_file)

def map_mask(tod, model, threshold=10):
    import linear_operators as lo
    import numpy as np

    return lo.decimate(model.T(np.ones(tod.shape)) < threshold)

def usage():
    print(__usage__)

__usage__ = """Usage: pacs_dli [options] [config_file]

Use tamasis with the double loop inference algorithm (Seeger et al.).

[config_file] is the name of the configuration file.

Options:
  -h, --help        Show this help message and exit.
  -v, --verbose     Print status messages to standard output.
  -f, --filenames   Overrides filenames configuration file value.
  -o, --output      Overrides output default value.
"""

# to call from command line
if __name__ == "__main__":
    main()
