# mappers : standard steps to go from data to a map estimate

def tamasis_mapper_rls(filenames, output_file, keywords, verbose=False):
    """
    Perform regularized least-square inversion of a simple model (tod
    mask and projection). Processing steps are as follows :
    
        - define PacsObservation instance
        - get Time Ordered Data (get_tod)
        - define projection
        - perform inversion on model
        - save file

    Arguments
    ---------
    filenames: list of strings
        List of data filenames.
    output_file : string
        Name of the output fits file.
    keywords: dict
        Dictionary containing options for all steps as dictionary.
    verbose: boolean (default False)
        Set verbosity.

    Returns
    -------
    Returns nothing. Save result as a fits file.
    """
    from scipy.sparse.linalg import cgs
    import tamasis as tm
    # verbosity
    keywords["mapper_rls"]["verbose"] = verbose
    # define observation
    obs_keys = keywords.get("PacsObservation", {})
    obs = tm.PacsObservation(filenames, **obs_keys)
    # get data
    tod_keys = keywords.get("get_tod", {})
    tod = obs.get_tod(**tod_keys)
    # define projector
    proj_keys = keywords.get("Projection", {})
    projection = tm.Projection(obs, **proj_keys)
    # define mask
    masking_tod  = tm.Masking(tod.mask)
    # define full model
    model = masking_tod * projection
    # perform map-making inversion
    mapper_keys = keywords.get("mapper_rls", {})
    map_rls = tm.mapper_rls(tod, model, solver=cgs, **mapper_keys)
    # save
    map_rls.save(output_file)

def tamasis_mapper_invntt(filenames, output_file, keywords, verbose=False):
    """
    Perform regularized least-square inversion of state of the art
    PACS model.  The PACS model includes tod mask, compression,
    response, projection and noise covariance (invntt).
    Processing steps are as follows :
    
        - define PacsObservation instance
        - mask first scanlines to avoid drift
        - get Time Ordered Data (get_tod)
        - 2nd level deglitching (with standard projector and tod median
            filtering with a  narrow window)
        - median filtering
        - define projection
        - perform inversion on model
        - save file

    Arguments
    ---------
    filenames: list of strings
        List of data filenames.
    output_file : string
        Name of the output fits file.
    keywords: dict
        Dictionary containing options for all steps as dictionary.
    verbose: boolean (default False)
        Set verbosity.

    Returns
    -------
    Returns nothing. Save result as a fits file.
    """
    from scipy.sparse.linalg import cgs
    import tamasis as tm
    # verbosity
    keywords["mapper_rls"]["verbose"] = verbose
    # define observation
    obs = tm.PacsObservation(filenames, **keywords["PacsObservation"])
    # extra masking
    scanline_masking(obs, **keywords["scanline_masking"])
    # get data
    tod = obs.get_tod(**keywords["get_tod"])
    # degltiching
    deglitching(obs, tod, **keywords["degltiching"])
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
    # perform map-making inversion
    map_rls = tm.mapper_rls(tod, model, invntt=invntt, solver=cgs,
                            **keywords["mapper_rls"])
    # save
    map_rls.save(output_file)

# processing steps

def scanline_masking(obs, n_repetition=0, n_scanline=None):
    """
    Mask everything up to the n_scanline scan line of the n_repetition repetition.
    If n_scanline is None, mask the whole repetition.

    Arguments
    ----------

    obs: PacsObservation
        The considered PACS observation instance.

    n_repetition: int
        Mask everything up to n_repetition.

    n_scanline: int or None (default: None)
        Mask everything up to n_scanline of n_repetition.
        If None mask n_repetition entirely.

    Returns
    -------

    Returns nothing, obs is altered inplace.

    Note
    ----

    scan line index starts at 1 !
    """
    if n_scanline == None:
        for i in xrange(n_repetition):
            obs.pointing.masked[obs.status.Repetition == i] = True
    else:
        for i in xrange(n_repetition - 1):
            obs.pointing.masked[obs.status.Repetition == i] = True
        for j in xrange(1, n_scanline):
            obs.pointing.masked[(obs.status.Repetition == n_repetition) *
                                (obs.status.ScanLineNumber == j)] = True

def deglitching(obs, tod, length=100, nsigma=25., method="mad"):
    """
    Include all the deglitching steps : 

       - define a sharp projector without oversampling and with
         default resolution.

       - perform tod filtering

       - second level deglitching with mad or std

    Arguments
    ---------
    obs : PacsObservation instance
    tod : Time Ordered Data
    length : int (default: 100)
        Median filtering window size.
    nsigma : float (default: 25.)
        Median filtering threshold
    method : "mad" or "std" (default: "mad")
        Filtering method (median absolute deviation or standard deviation).

    Returns
    -------
    Returns nothing. tod mask is updated in-place.
    """
    import tamasis as tm

    if method == "mad":
        method = tm.deglitch_l2mad
    elif method == "std":
        method = tm.deglitch_l2std
    else:
        raise ValueError("Unrecognized deglitching method.")

    # special deglitching projector
    proj_glitch = tm.Projection(obs,
                                method='sharp',
                                oversampling=False,
                                npixels_per_sample=6)
    # filter tod with narrow window
    tod_glitch = tm.filter_median(tod, length=length)
    # actual degltiching according to selected method (mad or std)
    tod.mask = method(tod_glitch, proj_glitch, nsigma=nsigma)

