"""Explicit PSI-MS spectrum types that do not represent mass spectra."""

# The EMR parent and its emission/absorption children describe intensities over
# wavelength, not m/z. Keep these explicit: names, units and missing MS levels
# alone must never cause an otherwise unlabelled mass spectrum to be skipped.
NON_MS_SPECTRUM_TERMS = frozenset({"MS:1000804", "MS:1000805", "MS:1000806"})
