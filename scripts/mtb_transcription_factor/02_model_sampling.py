"""
Generate flux samples from iMAT models for all the TFOE strains
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib

# External Imports
import cobra  # type:ignore
from cobra.sampling import OptGPSampler  # type:ignore
import metworkpy  # type:ignore

# Local Imports

# Setup Path
try:
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
except NameError:
    BASE_PATH = pathlib.Path(".").absolute()
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
TFOE_MODELS_PATH = CACHE_PATH / "tf_models"
FLUX_SAMPLES_PATH = CACHE_PATH / "tf_model_flux_samples"

# Make directories if required
CACHE_PATH.mkdir(parents=True, exist_ok=True)
TFOE_MODELS_PATH.mkdir(parents=True, exist_ok=True)
FLUX_SAMPLES_PATH.mkdir(parents=True, exist_ok=True)

# Setup Logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=BASE_PATH
    / "logs"
    / "transcription_factors"
    / "02_model_sampling.log",
    filemode="w",
    level=logging.INFO,
)


# Run Parameters
PROCESSES = 12
THINNING = 500
NUM_SAMPLES = 1_000

# Set cobra to use the hybrid solver
cobra.Configuration().solver = "hybrid"


# Generate samples for all the models
logger.info("Generating samples for all the tf models")
for model_path in (TFOE_MODELS_PATH).glob("*.json"):
    model_name = model_path.stem
    logger.info(f"Sampling from {model_name}")
    flux_sample_out_path = FLUX_SAMPLES_PATH / f"{model_name}.parquet"
    if flux_sample_out_path.exists():
        continue  # Already sampled from this
    model = metworkpy.read_model(model_path)
    logger.info("Starting OptGP Sampler")
    sampler = OptGPSampler(model=model, thinning=THINNING, processes=PROCESSES)
    samples = sampler.sample(NUM_SAMPLES)
    valid_samples = samples[sampler.validate(samples) == "v"]  # type:ignore
    logger.info(
        f"Validated samples, {len(valid_samples) / NUM_SAMPLES:.2%} valid"
    )
    logger.info("Saving samples")
    valid_samples.to_parquet(
        flux_sample_out_path,
        index=False,
    )

print("Finished Sampling :-)")
