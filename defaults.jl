# RESK's default parameters
# ------------------------------------------------
# Edit these to your preference.
# ================================================


# Space bounds
# ------------------------------------------------
DEF_X_MAX = 100
DEF_Y_MAX = 10
DEF_Z_MAX = 10


# Population parameters
# ------------------------------------------------
DEF_N_DEMES_STARTFILL = 5 # Number of demes filled at start
DEF_CAPACITY = 20
DEF_PROLIF_RATE = 1.8 # Proliferation rate


# Gene parameters
# ------------------------------------------------
DEF_MUT_RATE = 0.7567 # Genome-wide mutation rate
DEF_MIGR_RATE = 0.1 # Migration rate
DEF_SEL_COEF = 0.002 # Selection coefficient
DEF_PROP_OF_DEL_MUTS = 0.9 # Proportion of deleterious mutations in nature

# Finite-sites
DEF_DOMIN_COEF = 0 # Dominance coefficient
DEF_N_SEL_LOCI = 500 # Number of selected loci
DEF_N_LOCI = 1000 # Number of loci

# Infinite-sites
DEF_N_SEGR_REGIONS = 20 # Number of segregating regions


# Expansion parameters
# ------------------------------------------------
DEF_X_MAX_BURNIN = 5
DEF_R_MAX_BURNIN = 3 # Radius that bounds the burn-in area
DEF_N_GENS_BURNIN = 10 # Number of burn-in generations
DEF_X_MAX_EXP = DEF_X_MAX
DEF_Y_MAX_EXP = DEF_Y_MAX
DEF_R_MAX_EXP = 20 # Radius that bounds the expansion area
DEF_N_GENS_EXP = 40 # Number of expansion generations
DEF_MIGR_MODE = "ort" # Migration mode
DEF_DATA_TO_GENERATE = "FP"


# Other
# ------------------------------------------------
DEF_AUTO_GRAPHS = false # Scale graphs automatically?