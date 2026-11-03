"""MassCubeX: streaming MS1 feature extraction. Licensed under Apache-2.0."""

__version__ = "0.5.1"
__all__ = ["Config", "ProcessingResult", "process_file", "__version__"]


def __getattr__(name):
    if name == "Config":
        from .config import Config

        return Config
    if name in {"ProcessingResult", "process_file"}:
        from . import pipeline

        return getattr(pipeline, name)
    raise AttributeError(name)
