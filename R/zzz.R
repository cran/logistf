
.onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE))
        emmeans::.emm_register(c("logistf"), pkgname)
}
