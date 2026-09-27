# Documentation for the CRC and IBD package datasets.
# One source file generates two linked help topics.

#' Age-stratified IBD microbiome data
#'
#' Processed genus-level microbiome counts and sample information for
#' inflammatory bowel disease (IBD) meta-analysis. The data contain 3,472
#' biopsy or stool samples from seven source cohorts, organized into ten
#' analysis studies. Cases include Crohn's disease and ulcerative colitis.
#' Repeated samples from the same participant are identified by
#' \code{IBD_data$subject_ids}.
#'
#' @format
#' A named list with seven elements. Shared fields use the same names as
#' \code{CRC_data}; IBD additionally includes \code{design}:
#' \describe{
#'   \item{counts}{A named list of ten numeric (double) count matrices. Rows are
#'     samples and columns are genus-level taxon identifiers. The number of
#'     retained taxa differs between studies (84--205). Counts are not
#'     converted to relative abundances.}
#'   \item{metadata}{A named list of ten sample metadata data frames, each
#'     with 60 columns. These include cohort and sample identifiers,
#'     \code{subject_accession}, \code{sample_type}, \code{disease}
#'     (\code{"CD"}, \code{"UC"}, or \code{"control"}), \code{gender},
#'     \code{race}, and \code{age}. Derived columns are \code{age.numeric}
#'     (age in years), \code{analysis_study} (cohort with HMP2 age strata),
#'     and \code{study_key} (analysis study and sample type). Additional
#'     clinical and sequencing fields may contain missing values.}
#'   \item{design}{A named list of ten design matrices, without intercept
#'     columns. The first column is \code{disease}; remaining columns are
#'     the study-specific adjustment covariates \code{gender}, \code{race},
#'     and/or \code{age_adult}.}
#'   \item{subject_ids}{A named list of ten factors derived from
#'     \code{metadata[[i]]$subject_accession}. Each factor has one entry per
#'     sample; equal values identify repeated samples within that study.
#'     These are cluster identifiers, not covariate design matrices.}
#'   \item{taxonomy}{A 249-by-6 character matrix. Rows are taxon identifiers;
#'     columns \code{Rank1} through \code{Rank6} give kingdom, phylum,
#'     class, order, family, and genus, respectively. Values retain prefixes
#'     such as \code{f__} and \code{g__}; a prefix without a name denotes an
#'     unclassified taxon at that rank. This is the full taxonomy table, so
#'     not all rows are retained in every study.}
#'   \item{disease}{A named list of ten one-column matrices containing the
#'     disease indicator: 1 for IBD (CD or UC) and 0 for controls.
#'     Each matrix equals \code{design[[i]][, "disease", drop = FALSE]}.}
#'   \item{adjust}{A named list of ten matrices containing the adjustment
#'     columns of \code{design}, excluding \code{disease}. Covariate sets
#'     are allowed to differ between studies.}
#' }
#'
#' @details
#' All six study-level lists have identical names and ordering. Within each
#' study, rows of \code{counts}, \code{metadata}, \code{design},
#' \code{disease}, and \code{adjust}, and entries of \code{subject_ids},
#' refer to the same samples in the same order. Match count-matrix column
#' names to \code{rownames(IBD_data$taxonomy)}; do not align taxa by column position
#' across studies.
#'
#' HMP2 is divided into pediatric (age below 18 years) and adult (age at
#' least 18 years) studies. HMP2 samples with missing age are excluded.
#' MucosalIBD is retained as one study: ages range from 4 to 18 years, and
#' the 18-year-old subgroup has only one control sample.
#'
#' The analysis studies and their sample sizes are:
#' \tabular{lr}{
#'   Study \tab Samples \cr
#'   CS-PRISM.biopsy \tab 175 \cr
#'   CS-PRISM.stool \tab 319 \cr
#'   Herfarth_CCFA_Microbiome_3B_combined.stool \tab 756 \cr
#'   HMP2.adult.biopsy \tab 78 \cr
#'   HMP2.pediatric.biopsy \tab 76 \cr
#'   Jansson_Lamendella_Crohns.stool \tab 614 \cr
#'   MucosalIBD.biopsy \tab 132 \cr
#'   Pouchitis.biopsy \tab 451 \cr
#'   RISK.biopsy \tab 620 \cr
#'   RISK.stool \tab 251
#' }
#'
#' Preparation excludes the seven samples with the lowest total read depths
#' in the source count table and the PROTECT, LSS-PRISM, and BIDMC-FMT
#' cohorts. Within each analysis study, taxa with zero counts in more than
#' 90 percent of samples are removed.
#'
#' In the saved design matrices, \code{gender} is coded 0 for female and
#' 1 otherwise; \code{race} is coded 0 for observed non-white values and
#' 1 otherwise; \code{age_adult} is coded 0 for age below 18 and 1 otherwise.
#' Gender and age are forward-filled before coding. Residual missing
#' gender, race, or age values are assigned to the category coded 1 by the
#' preprocessing procedure; the original metadata are retained. A covariate
#' is included only when fewer than half its source values are missing and
#' neither binary category accounts for 95 percent or more of samples.
#'
#' The dataset is stored as a single object, \code{IBD_data}, in
#' \code{IBD_data.rda}. Its list packaging and serialization conventions match
#' \code{CRC_data.rda}. Both objects have the fields \code{counts},
#' \code{metadata}, \code{subject_ids}, \code{taxonomy}, \code{disease},
#' and \code{adjust}. Only IBD stores \code{design}, which combines disease
#' and study-specific adjustment covariates. Source-specific metadata and taxonomy
#' labels are retained; the shared field names do not imply that the
#' cohorts have identical metadata columns or taxon identifiers.
#'
#' @seealso \code{\link{CRC_data}}
#'
#' @source
#' Prepared from the genus-level IBD count, metadata, and taxonomy files
#' using \code{IBD_process_age_split.R}, then packaged using
#' \code{package_data_export.R}.
#'
#' @examples
#' data("IBD_data")
#' names(IBD_data)
#' vapply(IBD_data$counts, nrow, integer(1))
#'
#' # Inspect the adult HMP2 study.
#' study <- "HMP2.adult.biopsy"
#' dim(IBD_data$counts[[study]])
#' table(IBD_data$disease[[study]][, "disease"])
#' nlevels(IBD_data$subject_ids[[study]])
#'
#' # Extract family and genus annotations for a composition analysis.
#' tax <- IBD_data$taxonomy[, c("Rank5", "Rank6"), drop = FALSE]
#' colnames(tax) <- c("Rank1", "Rank2")
#'
#' @docType data
#' @name IBD_data
#' @usage IBD_data
#' @keywords datasets
"IBD_data"


#' Processed colorectal cancer microbiome data
#'
#' Genus-level microbiome counts and sample information for colorectal
#' cancer (CRC) meta-analysis. The data contain 573 samples from five
#' country-based study groups, including 283 CRC cases and 290 controls.
#' Preprocessing follows \code{CRC_process.R}. Field names are aligned
#' with the processed IBD dataset where applicable.
#'
#' @format
#' A named list with six elements. Shared fields use the same names as
#' \code{IBD_data}; CRC does not store a redundant disease-only \code{design}:
#' \describe{
#'   \item{counts}{A named list of five numeric (double) count matrices. Rows are
#'     samples and columns are genus-level taxon identifiers. Each study
#'     retains 104--122 taxa after prevalence filtering. Counts are not
#'     converted to relative abundances or rescaled to a common depth.}
#'   \item{metadata}{A named list of five sample metadata data frames,
#'     each with 18 columns. These include \code{Sample_ID},
#'     \code{External_ID}, \code{Age}, \code{Gender}, \code{BMI},
#'     \code{Country}, \code{Study}, \code{Group}, \code{Library_Size},
#'     and \code{block}, together with tumor stage, localization,
#'     colonoscopy timing, fecal occult blood test (FOBT), diabetes,
#'     vegetarian diet, and smoking information. \code{Group} is
#'     \code{"CRC"} for cases and \code{"CTR"} for controls.
#'     Original metadata names and missing values are retained.}
#'   \item{subject_ids}{\code{NULL}. The original CRC analysis treats
#'     samples as independent and does not specify subject-level clusters.
#'     No cluster identifiers are inferred from \code{External_ID}, which
#'     includes repeated placeholder values.}
#'   \item{taxonomy}{A 133-by-6 character matrix. Rows are genus-level
#'     taxon identifiers; columns \code{Rank1} through \code{Rank6} give
#'     kingdom, phylum, class, order, family, and genus, respectively.
#'     The original annotation strings are preserved, including numeric
#'     prefixes and unresolved groups such as \code{"NA Enterobacteriaceae
#'     gen. [C Escherichia/Shigella]"}. Unlike the IBD taxonomy, these
#'     strings do not use \code{k__}, \code{f__}, or \code{g__} prefixes.
#'     This is the full genus taxonomy table; not all rows are retained
#'     in every study.}
#'   \item{disease}{A named list of five one-column matrices containing
#'     the disease indicator: 1 for CRC and 0 for controls, derived directly
#'     from \code{metadata[[i]]$Group}. The original CRC analysis does not
#'     adjust for age, gender, BMI, or other covariates.}
#'   \item{adjust}{A named list of five zero-column matrices, each with
#'     one row per sample. These retain the IBD-compatible field structure
#'     while indicating that no adjustment covariates were specified.}
#' }
#'
#' @details
#' The lists \code{counts}, \code{metadata}, \code{disease},
#' and \code{adjust} have identical study names and ordering. Within each
#' study, their rows refer to the same samples in the same order. Row names
#' use the original \code{metadata[[i]]$Sample_ID}; the syntactic versions
#' of sample IDs in the source count table are matched before processing.
#' Match count-matrix column names to \code{rownames(CRC_data$taxonomy)};
#' do not align taxa by column position across studies.
#'
#' Studies are grouped by the original \code{Country} field, as in
#' \code{CRC_process.R}, rather than by \code{block}. In particular, the
#' before- and after-colonoscopy blocks of the Chinese cohort remain in
#' one study. The retained sample sizes are:
#' \tabular{llrrr}{
#'   Country key \tab Cohort \tab Samples \tab CRC \tab Controls \cr
#'   AUS \tab AT-CRC \tab 109 \tab 46 \tab 63 \cr
#'   CHI \tab CN-CRC \tab 127 \tab 73 \tab 54 \cr
#'   FRA \tab FR-CRC \tab 113 \tab 52 \tab 61 \cr
#'   GER \tab DE-CRC \tab 120 \tab 60 \tab 60 \cr
#'   USA \tab US-CRC \tab 104 \tab 52 \tab 52
#' }
#'
#' The low-depth sample \code{CCIS12370844ST-4-0} (original row 19,
#' total genus count 611) is excluded. Within each study, taxa with zero
#' counts in more than 90 percent of samples are removed. The
#' \code{Library_Size} metadata field records the original sequencing
#' library size and need not equal the sum of retained genus counts.
#'
#' Taxonomy follows the source annotations. The first annotation for each
#' genus is retained, matching the original processing script. For
#' \code{375288 Parabacteroides}, this retains the source family annotation
#' \code{171551 Porphyromonadaceae} rather than the alternative annotation
#' found later in the source taxonomy table.
#'
#' The package dataset is stored as a single object, \code{CRC_data}, to
#' avoid overwriting identically named components of other datasets. It is
#' saved in \code{CRC_data.rda}, with the same list packaging and serialization
#' conventions as \code{IBD_data.rda}. The fields \code{counts},
#' \code{metadata}, \code{subject_ids}, \code{taxonomy}, \code{disease},
#' and \code{adjust} have matching names in both objects. CRC does not store
#' \code{design}, because its only predictor is already stored in \code{disease}.
#' Unlike IBD, CRC has no subject-level clustering or adjustment covariates:
#' \code{subject_ids} is \code{NULL} and \code{adjust} contains zero-column
#' matrices. Source-specific metadata and taxonomy labels are retained.
#'
#' @seealso \code{\link{IBD_data}}
#'
#' @source
#' Prepared from the CRC genus count, metadata, and taxonomy files using
#' \code{CRC_process_package.R}, which reproduces the data-preparation
#' rules in \code{CRC_process.R} without running association tests, then
#' packaged using \code{package_data_export.R}.
#'
#' @examples
#' data("CRC_data")
#' names(CRC_data)
#' vapply(CRC_data$counts, nrow, integer(1))
#'
#' # Inspect the French cohort.
#' study <- "FRA"
#' dim(CRC_data$counts[[study]])
#' table(CRC_data$disease[[study]][, "disease"])
#' identical(rownames(CRC_data$counts[[study]]),
#'           rownames(CRC_data$metadata[[study]]))
#'
#' # Extract family and genus annotations for a composition analysis.
#' tax <- CRC_data$taxonomy[, c("Rank5", "Rank6"), drop = FALSE]
#' colnames(tax) <- c("Rank1", "Rank2")
#'
#' @docType data
#' @name CRC_data
#' @usage CRC_data
#' @keywords datasets
"CRC_data"
