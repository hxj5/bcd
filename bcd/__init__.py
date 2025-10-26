# Benchmarking of CNA Detection.


from . import app
from . import call_accuracy

from .app import VERSION
from .call_accuracy.main import bcd_main as call_accuracy_main


# what `__all__` does:
# https://stackoverflow.com/questions/44834/what-does-all-mean-in-python
# https://stackoverflow.com/questions/42950256/how-to-import-private-functions-with-import
__all__ = ["__version__"]
__version__ = VERSION
