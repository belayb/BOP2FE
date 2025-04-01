.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is version ", packageVersion(pkgname), 
                        " of ", pkgname, ". For the latest version please check https://github.com/belayb/BOP2FE/. If you experiance an issue, please contact belaybirlie.yimer@astellas.com")
}