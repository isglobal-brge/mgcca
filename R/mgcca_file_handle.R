# ---- one HDF5 file held open for the whole of a reliability analysis -------
#
# `mgcca_sensitivity()` and `mgcca_stability()` touch their file through several
# C++ entry points in a row, and each of them opens the file, does its job and
# closes it again. Two consecutive calls are therefore a close followed by a
# reopen, and BigDataStatMeth runs a read-only pre-flight probe on every open
# it does not already know about -- the probe Windows refuses. The probe is
# skipped whenever the HDF5 instance doing the asking already holds the file
# open, so both functions take a handle of their own as soon as they have a
# context and give it up when they return: every open inside that window is
# probe-free, and no entry point had to change.
#
# Exactly one window owns the handle at a time. The registry below records the
# path that window holds, and `.mgcca_rel_release()` gives up the registry
# entry and the C++ handle in the same step, so the two cannot disagree about
# what is open.
.mgcca_held_file <- new.env(parent = emptyenv())
.mgcca_held_file$path   <- NULL
.mgcca_held_file$handle <- NULL

# Take the handle. Called AFTER the caller has armed its release, never before:
# no execution point may own an open file whose release is not yet scheduled.
.mgcca_rel_hold <- function(file) {
    if (!length(file) || is.na(file[1]) || !nzchar(file[1]))
        stop("no HDF5 file to hold open", call. = FALSE)
    # A read-write open CREATES a file it does not find, which a handle meant to
    # hold an existing analysis must never do: an absent file is an error here,
    # not an empty new one.
    if (!file.exists(file[1]))
        stop("the HDF5 file to hold open does not exist: ", file[1],
             call. = FALSE)
    if (!is.null(.mgcca_held_file$path))
        stop(.mgcca_nested_file_error(file[1]))
    handle <- mgcca_hold_file_rcpp(file[1])
    .mgcca_held_file$path   <- file[1]
    .mgcca_held_file$handle <- handle
    # A garbage-collection backstop only. The contract is the caller's
    # on.exit(); this catches a handle whose owner never got to run it.
    reg.finalizer(handle, .mgcca_rel_release)
    handle
}

# The single place a held file is given up: the caller's on.exit() and the
# garbage collector's finalizer both come here. Idempotent, and a no-op on a
# context that holds nothing, so it is safe to arm before the handle exists.
# It never raises: a finalizer that failed would be reported out of any context
# a user could act on.
.mgcca_rel_release <- function(ctx = NULL) {
    handle <- if (inherits(ctx, "externalptr")) ctx
              else if (is.environment(ctx) || is.list(ctx)) ctx$handle
              else NULL
    # Nothing of its own to give up: a context whose handle was refused, or one
    # that never asked for one, must not release somebody else's.
    if (is.null(handle)) return(invisible(NULL))
    # The registry entry and the C++ handle go together, and only THIS
    # handle's entry goes: a finalizer firing late, after a later analysis has
    # taken a handle of its own, must not clear that one's entry.
    if (identical(.mgcca_held_file$handle, handle)) {
        .mgcca_held_file$path   <- NULL
        .mgcca_held_file$handle <- NULL
    }
    if (is.environment(ctx)) ctx$handle <- NULL
    try(mgcca_release_file_rcpp(handle), silent = TRUE)
    invisible(NULL)
}

# ---- the reverse collision -------------------------------------------------
# A file mgcca holds open through its own statically linked HDF5 cannot be
# opened through BigDataStatMeth's, and the other way round: the instance that
# has to see the open is the one doing the asking. No reliability analysis
# reaches a BigDataStatMeth open at R level, but a caller that nested one
# inside a held window would get "(checkHDF5File) Cannot open file: H5Fopen
# failed" and nothing to act on. So every R-level BigDataStatMeth open in this
# package asks this first, immediately before it opens.
.mgcca_refuse_held_file <- function(file, who) {
    held <- .mgcca_held_file$path
    if (is.null(held) || !length(file) || is.na(file[1]) || !nzchar(file[1]))
        return(invisible(NULL))
    if (!identical(normalizePath(file[1], mustWork = FALSE), held))
        return(invisible(NULL))
    text <- paste0(who, ": an mgcca HDF5 handle is still open on this file, ",
                   "so it cannot be opened again here ('", held, "'). Let the ",
                   "analysis holding it finish first.")
    stop(errorCondition(text, class = "mgcca_file_held", call = NULL))
}

# The error handler of a reader that reports a failure as "this file does not
# carry it". A refusal because mgcca itself holds the file open is not an
# absence, so it goes through instead of becoming a silent NULL.
.mgcca_rethrow_held_file <- function(e, value = NULL) {
    if (inherits(e, "mgcca_file_held")) stop(e)
    value
}

.mgcca_nested_file_error <- function(file) {
    text <- paste0("an mgcca HDF5 handle is already open on '",
                   .mgcca_held_file$path, "', and only one is held at a ",
                   "time. The analysis asking for '", file, "' has to run ",
                   "on its own, not inside another one.")
    errorCondition(text, class = "mgcca_nested_file_handle", call = NULL)
}
