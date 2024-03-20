#' Functions called when package is attached
#' 
#' Calls to functions to be loaded on package attachment.
#' 
#' @param libname Defaults is find.package( "MpoxModelPack")
#' @param pkgname Default is "MpoxModelPack"
#' @return None
#' @export

.onAttach <- function( libname = find.package( "MpoxModelPack"), 
                       pkgname = "MpoxModelPack" ) {
  init.pop.fn(c("mtl", "trt", "van"), 1)
  load.params.fn(VE = 0.5149,
                 contact_prop = 0.20)
}
