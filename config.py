from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parent

with open(ROOT / "config.yaml") as f:
    CONFIG = yaml.safe_load(f)
