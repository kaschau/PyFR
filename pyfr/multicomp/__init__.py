import functools
import re

# Physical constants
RU = 8314.46261815324         # Universal gas constant [J/(kmol·K)]
AVOGADRO = 6.02214076e26     # Avogadro's number [1/kmol]
KB = 1.380649e-23            # Boltzmann constant [J/K]
EPS0 = 8.8541878128e-12      # Vacuum permittivity [F/m]


def clean_csigns(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        s = fn(*args, **kwargs)
        s = re.sub(r'\+\s*\-', '- ', s)
        s = re.sub(r'\-\s*\-', '+ ', s)
        s = re.sub(r'\-\s*\+', '- ', s)
        s = re.sub(r'\+\s*\+', '+ ', s)
        s = re.sub(r'=\s*\+\s*', '= ', s)
        return s
    return wrapper
