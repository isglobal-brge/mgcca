## Selectable outputs: which components a collection actually reads.
##
## `mgcca_results()` materialises nine components out of the HDF5 file. Most
## callers want all of them and pay nothing for it -- they are small. Some do
## not: a plotting script wants `Y`, a projection wants `weights` and `scaling`,
## and a fit with many wide tables makes the per-table blocks the expensive part
## of the read. `outputs =` lets those callers say so.
##
## The rule the rest of this file exists to enforce is that the choice is purely
## SUBTRACTIVE. It decides what is READ. It cannot decide what was COMPUTED --
## the estimator finished before any of this runs, and the file is the same file
## either way -- so any component that comes back must be bit-for-bit the one a
## full collection would have produced. And it cannot RESURRECT anything: asking
## for `scores` from a fit that was run with `scores = FALSE` still gets NULL,
## because those datasets do not exist.

## The nine components, in the order the returned object carries them, and the
## element of that object each one fills.
.mgcca_output_names <- c("Y", "corsY", "scores", "pval", "weights",
                         "scaling", "AVE", "eigen", "overlap")

## Synonyms, so the object's own element names also work as requests.
.mgcca_output_synonym <- c(pval.cor = "pval", eig = "eigen")

## Canonicalise a user request: resolve synonyms, drop duplicates, refuse
## anything unknown. Order is the user's -- it is a record of what was asked.
.mgcca_output_canonical <- function(outputs) {
    if (!is.character(outputs))
        stop("'outputs' must be a character vector of component names, one or ",
             "more of: ", paste(.mgcca_output_names, collapse = ", "))
    if (!length(outputs))
        stop("'outputs' must name at least one component; use outputs = NULL ",
             "to collect everything")
    o   <- unique(outputs)
    hit <- match(o, names(.mgcca_output_synonym))
    o[!is.na(hit)] <- unname(.mgcca_output_synonym[hit[!is.na(hit)]])
    o   <- unique(o)
    bad <- setdiff(o, .mgcca_output_names)
    if (length(bad))
        stop("unknown output(s): ", paste(bad, collapse = ", "),
             ". Valid names are: ", paste(.mgcca_output_names, collapse = ", "),
             " (plus the synonyms pval.cor and eig)")
    o
}

## The read plan: a named logical over all nine components. NULL means all of
## them, which is what every version before this one did.
.mgcca_output_plan <- function(outputs, scores, pval) {
    want <- if (is.null(outputs)) .mgcca_output_names
            else .mgcca_output_canonical(outputs)
    p <- stats::setNames(.mgcca_output_names %in% want, .mgcca_output_names)
    # A request never overrides how the fit was run.
    if (!isTRUE(scores)) p[["scores"]] <- FALSE
    if (!isTRUE(pval))   p[["pval"]]   <- FALSE
    p
}
