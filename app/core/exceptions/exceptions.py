class NonIntegerItemException(Exception):
    """The list contains non-integers."""
    pass


class NonDNAException(Exception):
    """The sequence contains non-DNA characters."""
    pass

  
class MissingDataException(Exception):
    """Data is missing from the file."""
    pass
