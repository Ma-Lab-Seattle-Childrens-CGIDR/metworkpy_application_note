"""
Generate flux samples from iMAT models for all the TFOE strains
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type:ignore
from cobra.sampling import OptGPSampler  # type:ignore
import metworkpy  # type:ignore

# Local Imports

# Setup Path
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
TFOE_MODELS_PATH = CACHE_PATH / "tf_models"
FLUX_SAMPLES_PATH = CACHE_PATH / "tf_model_flux_samples"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

# Make directories if required
CACHE_PATH.mkdir(parents=True, exist_ok=True)
FLUX_SAMPLES_PATH.mkdir(parents=True, exist_ok=True)

# Setup Logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "02_model_sampling.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Set cobra to use the hybrid solver
cobra.Configuration().solver = CONFIG["cobra"]["solver"]


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
    sampler = OptGPSampler(
        model=model,
        thinning=CONFIG["mtb_tf"]["sampling"]["thinning"],
        processes=CONFIG["processes"],
    )
    samples = sampler.sample(CONFIG["mtb_tf"]["sampling"]["num-samples"])
    valid_samples = samples[sampler.validate(samples) == "v"]  # type:ignore
    logger.info(
        f"Validated samples, {
            len(valid_samples)
            / CONFIG['mtb_tf']['sampling']['num-samples']:.2%} valid"
    )
    logger.info("Saving samples")
    valid_samples.to_parquet(
        flux_sample_out_path,
        index=False,
    )
logger.info("Finished Sampling! ;)")
