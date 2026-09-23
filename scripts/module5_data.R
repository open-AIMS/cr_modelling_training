# ---------------------------------------------------------------------------
# The datasets module 5 fits, prepared exactly as the module shows them.
#
# Sourced by scripts/generate_module5_fits.R, which makes the fits the module
# loads, and by scripts/screen_module5_models.R, which fits every admissible
# equation to the same data to choose which equations the module uses. Keeping
# the preparation in one file means the screen and the taught fits cannot
# drift onto different data.
# ---------------------------------------------------------------------------

data(coral_pam, package = "bayesnec")
data(alga, package = "bayesnec")
data(lum31, package = "bayesnec")
data(nassarius, package = "bayesnec")

# A log predictor needs the control somewhere other than log(0). One decade
# below the lowest tested level is the convention the module uses throughout.
pam_data <- subset(coral_pam, climate == "2018")
pam_data$log.x <- log(pam_data$diuron + min(pam_data$diuron[pam_data$diuron > 0]) / 10)

# Cladocopium proliferum against contaminant B, the series module 5 fits. It
# replaced contaminant A on 2026-09-23, after scripts/screen_module5_models.R
# showed that no equation described contaminant A's counts or growth rate
# within check_fit(): its response falls steeply over the top doses and then
# rises from 420 to 958 cells, and the single dispersion parameter inflated to
# absorb that left every equation simulating five to ten times the observed
# spread at the control. Contaminant B declines monotonically, with a
# coefficient of variation between 0.02 and 0.09.
alga_data <- subset(alga, species == "c_proliferum" & contaminant == "B")
alga_data$log.x <- log(alga_data$dose + min(alga_data$dose[alga_data$dose > 0]) / 10)

# Contaminant A, kept so the screen that rejected it can be rerun.
alga_data_a <- subset(alga, species == "c_proliferum" & contaminant == "A")
alga_data_a$log.x <- log(alga_data_a$dose + min(alga_data_a$dose[alga_data_a$dose > 0]) / 10)

# One plate rather than the whole series, because module 8 is where the
# plate-to-plate structure is modelled. This plate is chosen because every
# reading on it is positive, which a gamma requires; the copper plates include
# readings floored at zero. Only the 15-minute reading is used: each well is
# read at 15 and at 30 minutes, and pooling the two treats a repeat reading of
# one well as an independent replicate while confounding the spread at high
# concentrations with exposure time (RF, 2026-09-23).
lum_data <- subset(lum31, plate == "Zn 28Mar24 Rep1B" & minutes == 15)
lum_data$log.x <- log(lum_data$conc)

# The binomial example is a tracked CSV rather than a packaged dataset, because
# it predates the others and is already real: counts of survivors out of a
# total, from Gerard Ricardo's repository.
binom_data <- read.csv("vignettes/example_binomial.csv")
binom_data$log.x <- log(binom_data$raw_x)

# The hurdle example. Growth is referenced to a baseline mean rather than to
# each snail's own starting size, so a few survivors are recorded at or below
# zero; those are measurement noise around a true value near zero rather than
# deaths, and are nudged off the boundary. The control is placed one decade
# below the lowest tested dose so that a log predictor is defined.
snail <- subset(nassarius, contaminant == "A")
snail_pos_min <- min(snail$growth[snail$growth > 0])
snail$growth <- ifelse(snail$alive == 1 & snail$growth <= 0,
                       snail_pos_min / 2, snail$growth)
snail$log_dose <- log(snail$dose + min(snail$dose[snail$dose > 0]) / 10)
