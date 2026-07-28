"""
citations.py
=============

Utilities for discovering the DOIs attached to the classes/functions used
(directly or indirectly) by a given callable, and fetching bibtex citations
for them.

The primary entry point is :func:`gather_bibtex`, which statically walks the
call graph rooted at a function, method, or class -- following calls made in
its body, the methods of any classes involved, and their base classes -- and
collects the ``DOI`` class attribute wherever one is defined. It then queries
doi.org's content-negotiation endpoint for a bibtex entry for each DOI found.
"""

import inspect

import requests

DOI_BIBTEX_URL = "https://doi.org/{doi}"


def _resolve_name(name, obj):
    """
    Try to resolve `name` (a dotted attribute chain's root name) to a live
    object, using whatever namespaces are reachable from `obj` (a function,
    method, or class): globals, closure variables, and -- if `obj` is itself
    a method -- the enclosing class's attributes.
    """
    func = inspect.unwrap(obj.__func__) if inspect.ismethod(obj) else obj

    if inspect.isfunction(func):
        if name in func.__globals__:
            return func.__globals__[name]
        closure_names = func.__code__.co_freevars
        if name in closure_names and func.__closure__:
            idx = closure_names.index(name)
            try:
                return func.__closure__[idx].cell_contents
            except ValueError:
                return None
    return None


def _class_methods(cls):
    """Yield every function object defined directly in `cls.__dict__`."""
    for value in vars(cls).values():
        if isinstance(value, (staticmethod, classmethod)):
            value = value.__func__
        if inspect.isfunction(value):
            yield value


def doi_to_bibtex(doi, timeout=10):
    """
    Fetch a bibtex-formatted citation for a single DOI using doi.org's
    content-negotiation service.

    :param doi: The DOI to look up (e.g. "10.1086/145971").
    :type doi: str
    :param timeout: Request timeout in seconds.
    :type timeout: float
    :return: The bibtex entry text.
    :rtype: str
    :raises RuntimeError: If the DOI cannot be resolved to a bibtex entry.
    """
    url = DOI_BIBTEX_URL.format(doi=doi)
    headers = {"Accept": "application/x-bibtex"}
    try:
        response = requests.get(url, headers=headers, timeout=timeout)
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Failed to fetch bibtex for DOI {doi}: {e}") from e

    if not response.ok:
        raise RuntimeError(
            f"Failed to fetch bibtex for DOI {doi}: HTTP {response.status_code}"
        )
    return response.text.strip()
