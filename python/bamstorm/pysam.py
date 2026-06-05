"""
Drop-in pysam compatibility shim.

Replace:
    import pysam

With:
    import bamstorm.pysam as pysam

AlignmentFile and AlignedSegment have the same constructor signatures and
field names as their pysam counterparts.  Unsupported methods (pileup, mate,
find_introns, etc.) raise AttributeError as usual.
"""

from bamstorm._core import AlignmentFile, BamRecord as AlignedSegment, count

__all__ = ["AlignmentFile", "AlignedSegment", "count"]
