"""
citations.py
=============

Utilities for discovering the DOIs attached to the classes/functions used
(directly or indirectly) by a given callable, and fetching bibtex citations
for them.

The primary entry point is :func:`gather_bibtex`, which statically walks the
call graph rooted at a function, method, or class -- following calls made in
its body, the methods of any classes involved, and their base classes -- and
collects the ``DOI`` class attribute wherever one is defined. That attribute
may be a single DOI string or a list/tuple of DOI strings, for work that
should be credited to more than one paper. It then queries doi.org's
content-negotiation endpoint for a bibtex entry for each DOI found.
"""

import ast
import inspect
import textwrap

import requests

__all__ = ["find_dois", "doi_to_bibtex", "gather_bibtex"]

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


def _iter_referenced_names(func):
    """
    Parse the source of `func` and yield every plain name it references
    (loads). This covers direct calls (`Qux()`, `self.foo()`), as well as
    classes/functions that are only referenced indirectly, e.g. stored in a
    dict and dispatched through a lookup (`registry[key](**kwargs)`).
    """
    try:
        source = inspect.getsource(func)
    except (OSError, TypeError):
        return
    source = textwrap.dedent(source)
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return

    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load):
            yield node.id


def _iter_doi_strings(doi):
    """
    Normalize a ``DOI`` class attribute into individual DOI strings. The
    attribute may be a single string, or a list/tuple of strings for work
    that needs several papers cited. Anything else (or an empty entry) is
    ignored.
    """
    if isinstance(doi, str):
        candidates = [doi]
    elif isinstance(doi, (list, tuple, set, frozenset)):
        candidates = doi
    else:
        return

    for candidate in candidates:
        if isinstance(candidate, str) and candidate.strip():
            yield candidate.strip()


def _class_methods(cls):
    """Yield every function object defined directly in `cls.__dict__`."""
    for value in vars(cls).values():
        if isinstance(value, (staticmethod, classmethod)):
            value = value.__func__
        if inspect.isfunction(value):
            yield value


def find_dois(target, _seen=None):
    """
    Recursively walk the call graph rooted at `target` (a function, bound/
    unbound method, or class) and collect every ``DOI`` class attribute
    encountered along the way. A ``DOI`` attribute holding a list/tuple of
    strings contributes each of its elements individually.

    A class is "visited" by inspecting its own methods (including
    ``__init__``) and its base classes. A function/method is "visited" by
    statically finding the names it calls and resolving+recursing into
    whichever of those names refer to functions or classes.

    :param target: The function, method, or class to start the search from.
    :param _seen: Internal set of already-visited objects, used to avoid
        infinite recursion on circular/self-referential call graphs.
    :return: Set of DOI strings found.
    :rtype: set[str]
    """
    if _seen is None:
        _seen = set()

    key = id(target)
    if key in _seen:
        return set()
    _seen.add(key)

    dois = set()

    if inspect.isclass(target):
        # DOIs are collected from anywhere in the MRO, since an inherited
        # DOI still applies. Method bodies, however, are only inspected for
        # methods the class defines itself -- recursing into *every* method
        # of every base class would also pull in unrelated sibling classes
        # (e.g. a base class factory method that references other
        # subclasses), which is not part of what `target` actually uses.
        for klass in target.__mro__:
            if klass is object:
                continue
            dois |= set(_iter_doi_strings(klass.__dict__.get("DOI")))
        for method in _class_methods(target):
            dois |= find_dois(method, _seen)
        return dois

    if inspect.isfunction(target) or inspect.ismethod(target):
        for name in _iter_referenced_names(target):
            resolved = _resolve_name(name, target)
            if resolved is None:
                continue
            if inspect.isclass(resolved) or inspect.isfunction(resolved):
                dois |= find_dois(resolved, _seen)
        return dois

    return dois


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


def gather_bibtex(target, timeout=10):
    """
    Find every DOI attached to classes/functions used (directly or
    indirectly) by `target`, and return a single string containing bibtex
    entries for all of them.

    :param target: The function, method, or class to start the search from.
    :param timeout: Per-request timeout in seconds, passed to
        :func:`doi_to_bibtex`.
    :return: Newline-separated bibtex entries, one per DOI found.
    :rtype: str
    """
    dois = find_dois(target)
    entries = []
    for doi in sorted(dois):
        try:
            entries.append(doi_to_bibtex(doi, timeout=timeout))
        except RuntimeError as e:
            entries.append(f"% Could not resolve DOI {doi}: {e}")
    return "\n\n".join(entries)
