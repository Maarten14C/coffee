# temporary file until the migration of most 'rintcal' functions to 'rice' is completed


# internal functions to speed up reading and writing files, using the data.table R package if present
fastread <- function(fl, ...)
  if("data.frame" %in% (.packages())) # some Macs have problems with this package
    as.data.frame(data.table::fread(fl), ...) else
      read.table(fl, ...)



fastwrite <- function(fl, ...)
  if("data.frame" %in% (.packages())) # some Macs have problems with this package
    data.table::fwrite(as.data.frame(fl), ...) else
      write.table(fl, ...)



#' @name glue.ccurves
#' @title Glue calibration curves
#' @description Produce a custom curve by merging two calibration curves, e.g. a prebomb and a postbomb one for dates which straddle both curves.
#' @return The custom-made curve (invisibly)
#' @param prebomb The prebomb curve. Defaults to "IntCal20"
#' @param postbomb The postbomb curve. Defaults to "NH1" (Hua et al. 2013)
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="ccurves"}.
#' @examples
#' my.cc <- glue.ccurves()
#' @export
glue.ccurves <- function(prebomb="IntCal20", postbomb="NH1", cc.dir=c()) {
  glued <- rbind(ccurve(prebomb, FALSE, cc.dir=cc.dir), ccurve(postbomb, TRUE, cc.dir=cc.dir))
  glued <- glued[order(glued[,1]),]
  repeated <- which(diff(glued[,1]) == 0)
  if(length(repeated) > 0)
    invisible(glued[-repeated,]) else # remove any repeated years
      invisible(glued[order(glued[,1]),])
}



#' @name ccurve
#' @title Copy a calibration curve
#' @description Copy one of the calibration curves into memory.
#' @details Copy the radiocarbon calibration curve defined by cc into memory.
#' @return The calibration curve (invisible).
#' @param cc Calibration curve for 14C dates: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine),
#' \code{cc=3} for SHCal20 (southern hemisphere terrestrial). Alternatively, one can also write, e.g., "IntCal20", "Marine13". One can also make a custom-built calibration curve, e.g. using \code{mix.ccurves()}, and load this using \code{cc=4}. In this case, it is recommended to place the custom calibration curve in its own directory, using \code{cc.dir} (see below).
#' @param postbomb Use \code{postbomb=TRUE} to get a postbomb calibration curve (default \code{postbomb=FALSE}). For monthly data, type e.g. \code{ccurve("sh1-2_monthly")}
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="ccurves"}.
#' @param resample The IntCal curves come at a range of 'bin sizes'; every year from 0 to 5 kcal BP, then every 5 yr until 15 kcal BP, then every 10 yr until 25 kcal BP, and every 20 year thereafter. The curves can be resampled to constant bin sizes, e.g. \code{resample=5}. Defaults to FALSE.
#' @param glue If a postbomb curve is requested, it can be 'glued' to the pre-bomb curve. This feature is currently disabled - please use \code{glue.ccurves} instead
#' @examples
#' intcal20 <- ccurve(1)
#' marine20 <- ccurve(2)
#' shcal20 <- ccurve(3)
#' marine98 <- ccurve("Marine98")
#' pb.sh3 <- ccurve("sh3")
#' @references
#' Hammer and Levin 2017, "Monthly mean atmospheric D14CO2 at Jungfraujoch and Schauinsland from 1986 to 2016", heiDATA: Heidelberg Research Data Repository V2 \doi{10.11588/data/10100}
#'
#' Heaton et al. 2020 Marine20-the marine radiocarbon age calibration curve (0-55,000 cal BP). Radiocarbon 62, 779-820, \doi{10.1017/RDC.2020.68}
#'
#' Hogg et al. 2013 SHCal13 Southern Hemisphere Calibration, 0-50,000 Years cal BP. Radiocarbon 55, 1889-1903, \doi{10.2458/azu_js_rc.55.16783}
#'
#' Hogg et al. 2020 SHCal20 Southern Hemisphere calibration, 0-55,000 years cal BP. Radiocarbon 62, 759-778, \doi{10.1017/RDC.2020.59}
#'
#' Hua et al. 2013 Atmospheric radiocarbon for the period 1950-2010. Radiocarbon 55(4), \doi{10.2458/azu_js_rc.v55i2.16177}
#'
#' Hua et al. 2022 Atmospheric radiocarbon for the period 1950-2019. Radiocarbon 64(4), 723-745, \doi{10.1017/RDC.2021.95}
#'
#' Levin and Kromer 2004 The tropospheric 14CO2 level in mid latitudes of the Northern Hemisphere. Radiocarbon 46, 1261-1272
#'
#' Reimer et al. 2004 IntCal04 terrestrial radiocarbon age calibration, 0-26 cal kyr BP. Radiocarbon 46, 1029-1058, \doi{10.1017/S0033822200032999}
#'
#' Reimer et al. 2009 IntCal09 and Marine09 radiocarbon age calibration curves, 0-50,000 years cal BP. Radiocarbon 51, 1111-1150, \doi{10.1017/S0033822200034202}
#'
#' Reimer et al. 2013 IntCal13 and Marine13 radiocarbon age calibration curves 0-50,000 years cal BP. Radiocarbon 55, 1869-1887, \doi{10.2458/azu_js_rc.55.16947}
#'
#' Reimer et al. 2020 The IntCal20 Northern Hemisphere radiocarbon age calibration curve (0-55 cal kBP). Radiocarbon 62, 725-757, \doi{10.1017/RDC.2020.41}
#'
#' Stuiver et al. 1998 INTCAL98 radiocarbon age calibration, 24,000-0 cal BP. Radiocarbon 40, 1041-1083, \doi{10.1017/S0033822200019123}
#' @export
ccurve <- function(cc=1, postbomb=FALSE, cc.dir=NULL, resample=0, glue=FALSE) {
  if(postbomb) {
    if(cc==1 || tolower(cc) == "nh1")
      fl <- "postbomb_NH1.14C" else
      if(cc==2 || tolower(cc) == "nh2")
        fl <- "postbomb_NH2.14C" else
        if(cc==3 || tolower(cc) == "nh3")
          fl <- "postbomb_NH3.14C" else
          if(cc==4 || tolower(cc) == "sh1-2")
            fl <- "postbomb_SH1-2.14C" else
            if(cc==5 || tolower(cc) == "sh3")
              fl <- "postbomb_SH3.14C" else
              if(tolower(cc) == "nh1_monthly")
                fl <- "postbomb_NH1_monthly.14C" else
                if(tolower(cc) == "nh2_monthly")
                  fl <- "postbomb_NH2_monthly.14C" else
                  if(tolower(cc) == "nh3_monthly")
                    fl <- "postbomb_NH3_monthly.14C" else
                    if(tolower(cc) == "sh1-2_monthly")
                      fl <- "postbomb_SH1-2_monthly.14C" else
                      if(tolower(cc) == "sh3_monthly")
                        fl <- "postbomb_SH3_monthly.14C" else
                        if(tolower(cc) == "kure")
                          fl <- "Kure.14C" else
                          if(tolower(cc) == "levinkromer")
                            fl <- "LevinKromer.14C" else
                            if(tolower(cc) == "santos")
                            fl <- "Santos.14C" else
                              stop("cannot find this postbomb curve\n", call.=FALSE)
  } else
    if(cc==1 || tolower(cc) == "intcal20")
      fl <- "3Col_intcal20.14C" else
      if(cc==2 || tolower(cc) == "marine20")
        fl <- "3Col_marine20.14C" else
        if(cc==3 || tolower(cc) == "shcal20")
          fl <- "3Col_shcal20.14C" else
          if(cc==4 || tolower(cc) == "mixed")
            fl <- "mixed.14C" else
            if(tolower(cc) == "intcal13")
              fl <- "3Col_intcal13.14C" else
              if(tolower(cc) == "marine13")
                fl <- "3Col_marine13.14C" else
                if(tolower(cc) == "shcal13")
                  fl <- "3Col_shcal13.14C" else
                  if(tolower(cc) == "intcal09")
                    fl <- "3Col_intcal09.14C" else
                    if(tolower(cc) == "marine09")
                      fl <- "3Col_marine09.14C" else
                      if(tolower(cc) == "intcal04")
                        fl <- "3Col_intcal04.14C" else
                        if(tolower(cc) == "marine04")
                          fl <- "3Col_marine04.14C" else
                          if(tolower(cc) == "intcal98")
                            fl <- "3Col_intcal98.14C" else
                            if(tolower(cc) == "marine98")
                              fl <- "3Col_marine98.14C" else
                              if(tolower(cc) == "nh1")
                                fl <- "postbomb_NH1.14C" else
                                if(tolower(cc) == "nh2")
                                  fl <- "postbomb_NH2.14C" else
                                  if(tolower(cc) == "nh3")
                                    fl <- "postbomb_NH3.14C" else
                                    if(tolower(cc) == "sh1-2")
                                      fl <- "postbomb_SH1-2.14C" else
                                      if(tolower(cc) == "sh3")
                                        fl <- "postbomb_SH3.14C" else
                                         if(tolower(cc) == "nh1_monthly")
                                          fl <- "postbomb_NH1_monthly.14C" else
                                          if(tolower(cc) == "nh2_monthly")
                                            fl <- "postbomb_NH2_monthly.14C" else
                                            if(tolower(cc) == "nh3_monthly")
                                              fl <- "postbomb_NH3_monthly.14C" else
                                              if(tolower(cc) == "sh1-2_monthly")
                                                fl <- "postbomb_SH1-2_monthly.14C" else
                                                if(tolower(cc) == "sh3_monthly")
                                                  fl <- "postbomb_SH3_monthly.14C" else
                                                  if(tolower(cc) == "kure")
                                                    fl <- "kure.14C" else
                                                    if(tolower(cc) == "levinkromer")
                                                      fl <- "LevinKromer.14C" else
                                                      if(tolower(cc) == "santos")
                                                        fl <- "Santos.14C" else
                                                        if(tolower(cc) == "mixed")
                                                          fl <- "mixed.14C" else
                                                          stop("cannot find this curve", call.=FALSE)

  if(length(cc.dir) == 0)
    cc <- system.file("extdata/", fl, package='rintcal') else
      cc <- file.path(cc.dir, fl)
  cc <- fastread(cc)

  if(resample > 0) {
    yr <- seq(min(cc[,1]), max(cc[,1]), by=resample)
    mu <- approx(cc[,1], cc[,2], yr)$y
    er <- approx(cc[,1], cc[,3], yr)$y
    cc <- cbind(yr, mu, er)
  }
  invisible(cc)
}
