
class Iterator:
    """
    Protocol for iterating over an object.
    """

    def __init__(self, _iterable):
        """
        Store iterable and initialize count.

        Args:

            _iterable (iterable object) - must have a 'length' and be indexable

        """
        self._iterable = _iterable
        self.count = 0

    def __iter__(self):
        """ Iterate over contents. """
        return self

    def __next__(self):
        """ Advance iterator. """
        if self.count >= len(self._iterable):
            raise StopIteration
        else:
            self.count += 1
            return self._iterable[self.count - 1]
