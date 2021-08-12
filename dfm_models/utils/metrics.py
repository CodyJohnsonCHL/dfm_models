from numpy import abs, nanmax, nanmean, nanmin, nanstd, nanvar, sqrt


def rmse(model, ref):
    """
    root mean square error. Ideal value is zero

    model = numerical solution of shape M by N
    ref = analytical solution or observations of shape M by N

    returns rmse
    """
    e = model - ref
    e2 = e ** 2
    mse = nanmean(e2)
    rmse = sqrt(mse)
    return rmse


def nrmse(model, ref):
    """
    normalized (by data range) root mean square error. Ideal value is zero

    model = numerical solution of shape M by N
    ref = analytical solution or observations of shape M by N

    returns rmse
    """
    _rmse = rmse(model, ref)
    _range = nanmax(model) - nanmin(model)
    _nrmse = _rmse / _range
    return _nrmse


def r2(model, ref):
    """
    coefficient of determination. Ideal value is 1

    model = nuemrical solution of shape M by N
    ref = analytical solution or observation of shape M by N

    returns r2
    """
    e = model - ref
    e2 = e ** 2
    mse = nanmean(e2)
    var_ref = nanvar(ref)
    r2 = 1 - mse / var_ref
    return r2


def r(model, ref):
    """
    Pearson correlation coefficient.

    model = nuemrical solution of shape M by N
    ref = analytical solution or observation of shape M by N

    returns r
    """
    mod_res = model - nanmean(model)
    ref_res = ref - nanmean(ref)
    mod_sqr_res = mod_res ** 2
    ref_sqr_res = ref_res ** 2
    numerator = nanmean(mod_res * ref_res)
    denominator = sqrt(nanmean(mod_sqr_res) * nanmean(ref_sqr_res))
    r = numerator / denominator
    return r


def SI(model, ref):
    """
    Scatter index.

    model = nuemrical solution of shape M by N
    ref = analytical solution or observation of shape M by N

    returns SI
    """
    e = model - ref
    std_e = nanstd(e)
    mean_abs_ref = nanmean(abs(ref))
    SI = std_e / mean_abs_ref
    return SI


def bias(model, ref):
    """
    bias

    model = numerical solution of shape M by N
    ref = analytical solution of shape M by N

    returns bias
    """
    e = model - ref
    mean_e = nanmean(e)
    return mean_e


def nb(model, ref):
    """
    Normalzied bias.

    model = nuemrical solution of shape M by N
    ref = analytical solution or observation of shape M by N

    returns nb
    """
    e = model - ref
    mean_e = nanmean(e)
    mean_abs_ref = nanmean(abs(ref))
    nb = mean_e / mean_abs_ref

    return nb

