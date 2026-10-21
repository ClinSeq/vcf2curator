"""Annotation resource loaders and lookups (OncoKB, CGC, hotspot, curation).

Each loader turns an on-disk reference file into an in-memory lookup used to
annotate variants during parsing. Paths are supplied by the caller (CLI
arguments) rather than hardcoded.
"""
