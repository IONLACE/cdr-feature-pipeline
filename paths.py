from pathlib import Path
import re
from config import CONFIG, ROOT

class Paths:
    def __init__(self, root: Path, config: dict):
        self.root = root

        self.inp_pdb_dir = self._as_path(self._expand(config['paths']['inp_pdb_dir'], config), root)
        self.meta_dir = self._as_path(self._expand(config['paths']['meta_dir'], config), root)
        self.features_dir = self._as_path(self._expand(config['paths']['features_dir'], config), root)
        if config['paths'].get('logs_dir', config['paths'].get('log_dir')) is None:
            raise KeyError("Missing required config key: paths.log_dir (or paths.logs_dir)")
        self.logs_dir = self._as_path(
            self._expand(config['paths'].get('logs_dir', config['paths'].get('log_dir')), config),
            root,
        )
        self.log_dir = self.logs_dir


    def _expand(self, value: str, config: dict) -> str:
        pattern = re.compile(r"\$\{([^}]+)\}")

        def repl(match):
            key = match.group(1)
            if key not in config:
                raise KeyError(f"Missing config key for placeholder: {key}")
            return str(config[key])

        return pattern.sub(repl, value)

    def _as_path(self, value: str, root: Path) -> Path:
        path_value = Path(value)
        if path_value.is_absolute():
            return path_value
        return root / path_value
        


PATHS = Paths(ROOT, CONFIG) 
# print(f"Root directory: {PATHS.root}")
# print(f"Input PDB directory: {PATHS.inp_pdb_dir}")
# print(f"Meta directory: {PATHS.meta_dir}")
# print(f"Features directory: {PATHS.features_dir}")
# print(f"Logs directory: {PATHS.logs_dir}")
