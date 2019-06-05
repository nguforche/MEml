#' Aortic valve replacement surgery data from the \code{joineR} package
#' 
#' @description This is longitudinal data on an observational study on detecting
#'   effects of different heart valves, differing on type of tissue, implanted
#'   in the aortic position.  The data consists of longitudinal measurements (three cardiac functions)
#'   from patients who underwent aortic valve replacement from 1991 to 2001 at the
#'   Royal Brompton Hospital, London, United Kingdom. The data was first reported in [1]
#'   where the authors used all patients during the 10 years period with at least a year of follow  
#'   up with serial echocardiographic measurements and applied a linear mixed-effect model 
#'   to predict left ventricular mass index (LVMI). Similarly, the data was used in [2] to predict longitudinal 
#'   profile of LVMI categorized as high or normal using several patient baseline characteristics and 
#'   laboratory variables. LVMI is considered increased if LVMI >134 g/m 2 in male patients and 
#'   LVMI >110 g/m 2 in female patients, thus values in this range for both sex was considered 
#'   as the positive class in MEml.

#' @usage data(heart.valve)
#' @format This is a data frame in the unbalanced format, that is, with one row 
#'   per observation. The data consists in columns for patient identification, 
#'   time of measurements, longitudinal multiple longitudinal measurements, 
#'   baseline covariates, and survival data. The column names are identified as 
#'   follows:
#'   
#'   \describe{
#'   
#'   \item{\code{num}}{number for patient identification.}
#'   
#'   \item{\code{sex}}{gender of patient (\code{0 = }Male and \code{1 = 
#'   }Female).}
#'   
#'   \item{\code{age}}{age of patient at day of surgery (years).}
#'   
#'   \item{\code{time}}{observed time point, with surgery date as the time 
#'   origin (years).}
#'   
#'   \item{\code{fuyrs}}{maximum follow up time, with surgery date as the time 
#'   origin (years).}
#'   
#'   \item{\code{status}}{censoring indicator (\code{1 = }died and \code{0 = 
#'   }lost at follow up).}
#'   
#'   \item{\code{grad}}{valve gradient at follow-up visit.}
#'   
#'   \item{\code{log.grad}}{natural log transformation of \code{grad}.}
#'   
#'   \item{\code{lvmi}}{left ventricular mass index (standardised) at follow-up 
#'   visit.}
#'   
#'   \item{\code{log.lvmi}}{natural log transformation of \code{lvmi}.}
#'   
#'   \item{\code{ef}}{ejection fraction at follow-up visit.}
#'   
#'   \item{\code{bsa}}{preoperative body surface area.}
#'   
#'   \item{\code{lvh}}{preoperative left ventricular hypertrophy.}
#'   
#'   \item{\code{prenyha}}{preoperative New York Heart Association (NYHA) 
#'   classification (\code{1 = }I/II and \code{3 = }III/IV).}
#'   
#'   \item{\code{redo}}{previous cardiac surgery.}
#'   
#'   \item{\code{size}}{size of the valve (millimeters).}
#'   
#'   \item{\code{con.cabg}}{concomitant coronary artery bypass graft.}
#'   
#'   \item{\code{creat}}{preoperative serum creatinine (\eqn{\mu}mol/mL).}
#'   
#'   \item{\code{dm}}{preoperative diabetes.}
#'   
#'   \item{\code{acei}}{preoperative use of ace inhibitor.}
#'   
#'   \item{\code{lv}}{preoperative left ventricular ejection fraction (LVEF) 
#'   (\code{1 = }good, \code{2 = }moderate, and \code{3 = }poor).}
#'   
#'   \item{\code{emergenc}}{operative urgency (\code{0 = }elective, \code{1 = 
#'   }urgent, and \code{3 = }emergency).}
#'   
#'   \item{\code{hc}}{preoperative high cholesterol (\code{0 = }absent, \code{1 
#'   = }present treated, and \code{2 = }present untreated).}
#'   
#'   \item{\code{sten.reg.mix}}{aortic valve haemodynamics (\code{1 = }stenosis,
#'   \code{2 = }regurgitation, \code{3 = }mixed).}
#'   
#'   \item{\code{hs}}{implanted aortic prosthesis type (\code{1 = }homograft 
#'   and \code{0 = }stentless porcine tissue).}
#'   
#'   }
#' @keywords datasets
#' @seealso \code{\link{mental}}, \code{\link{liver}}, \code{\link{epileptic}},
#'   \code{\link{aids}}.
#' @source Mr Eric Lim (\url{http://www.drericlim.com})
#' @docType data
#' @references
#' 
#' Lim E, Ali A, Theodorou P, Sousa I, Ashrafian H, Chamageorgakis T, Duncan M, 
#' Diggle P, Pepper J. A longitudinal study of the profile and predictors of 
#' left ventricular mass regression after stentless aortic valve replacement. 
#' \emph{Ann Thorac Surg.} 2008; \strong{85(6)}: 2026-2029.
#
#' Che Ngufor,  Holly Van Houten, Brian S. Caffo , Nilay D. Shah, Rozalina G. McCoy 
#' Mixed Effect Machine Learning: a framework for predicting longitudinal change in hemoglobin A1c, 
#' in Journal of Biomedical Informatics, 2018 
"heart.valve"