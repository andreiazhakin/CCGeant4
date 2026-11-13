"""Load the CCPACK reader that ships beside the runtime CSV tools."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def load_read_pack(converter_dir=None):
    directory = Path(converter_dir) if converter_dir else Path(__file__).parent
    source = directory / "cc_binary.py"
    if not source.is_file():
        raise FileNotFoundError(source)
    spec = spec_from_file_location("_ccgeant4_pack_reader", source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load the CCPACK reader from {source}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.read_pack
