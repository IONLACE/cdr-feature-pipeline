from pathlib import Path
import os
import yaml

ROOT = Path(__file__).resolve().parent

with open(ROOT / "config.yaml") as f:
    CONFIG = yaml.safe_load(f)

# Allow container/runtime overrides without changing tracked config.yaml.
if os.environ.get("CDR_OUT_DIR"):
    CONFIG["out_dir"] = os.environ["CDR_OUT_DIR"]
