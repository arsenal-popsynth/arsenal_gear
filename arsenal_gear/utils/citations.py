"""
citations.py
============

Attach literature citations to classes and functions, and collect every
citation reachable from a given entry point as bibtex.

:func:`cite` decorates a class or a function with one or more identifiers --
a DOI (``"10.1086/145971"``) or, for work predating DOIs, a NASA ADS bibcode
(``"1996A&A...315..105R"``, written literally, no URL escaping)::

    @cite("10.1086/145971")
    class Salpeter(IMF):
        ...

:func:`gather_bibtex` then walks the call graph rooted at a function, method,
or class -- its body, the methods of any classes involved, and their base
classes -- and renders bibtex for every citation it finds along the way.

Bibtex text is only ever read, from the ``citations.bib`` shipped alongside
this module, so citations work offline and without credentials. Identifiers
are taken on trust here: checking their shape, looking them up on NASA ADS
and writing that file all belong to
:mod:`arsenal_gear.utils.refresh_citations`, which a maintainer runs after
adding or changing a ``@cite``::

    python -m arsenal_gear.utils.refresh_citations
"""

from __future__ import annotations

import ast
import inspect
import textwrap
from pathlib import Path

__all__ = [
    "REGISTRY",
    "CitationError",
    "cite",
    "find_citations",
    "gather_bibtex",
    "load_bibliography",
]

#: Attribute the :func:`cite` decorator writes its parsed identifiers to.
CITATIONS_ATTR = "__citations__"

#: Bibliography shipped with the package; the only source of bibtex text.
BIBLIOGRAPHY_PATH = Path(__file__).with_name("citations.bib")


class CitationError(Exception):
    """Raised when a citation cannot be resolved to bibtex."""


#: Every identifier seen by :func:`cite` since import. Populated as a side
#: effect of decoration, so once all of arsenal_gear has been imported this
#: is the complete set of citations used anywhere in the package -- which is
#: what the refresher walks. Updated in place, never rebound.
REGISTRY: set[str] = set()


def cite(*identifiers: str):
    """
    Decorator attaching citations to a class or function.

    The decorated object is returned unwrapped, with the identifiers
    recorded on it, which keeps it fully introspectable --
    :func:`find_citations` relies on being able to read the source of
    decorated functions. Stacking the decorator accumulates citations rather
    than replacing them.

    Identifiers are recorded as given, bar surrounding whitespace; whether
    they are well formed is the refresher's business.

    :param identifiers: One or more DOIs or ADS bibcodes.
    :return: A decorator that annotates and returns its argument.
    """
    parsed = tuple(i.strip() for i in identifiers)
    REGISTRY.update(parsed)

    def decorator(obj):
        # staticmethod/classmethod objects cannot carry attributes, so
        # annotate the function they wrap; `_class_methods` unwraps to match.
        target = obj.__func__ if isinstance(obj, (staticmethod, classmethod)) else obj
        # Read through __dict__ rather than getattr so that a class picks up
        # only its own citations here, never an inherited copy of a base
        # class's. Inheritance is handled deliberately in `find_citations`.
        existing = tuple(vars(target).get(CITATIONS_ATTR, ()))
        merged = existing + tuple(c for c in parsed if c not in existing)
        try:
            setattr(target, CITATIONS_ATTR, merged)
        except (AttributeError, TypeError) as e:
            raise CitationError(
                f"Cannot attach citations to {target!r}: it does not support "
                f"attribute assignment."
            ) from e
        return obj

    return decorator


def _own_citations(obj) -> tuple[str, ...]:
    """
    Identifiers attached directly to `obj`, ignoring anything it merely
    inherits.
    """
    return tuple(vars(obj).get(CITATIONS_ATTR, ()))


def _resolve_name(name, obj):
    """
    Try to resolve `name` to a live object, using the namespaces reachable
    from `obj` (a function or method): its globals, then its closure.
    """
    func = obj.__func__ if inspect.ismethod(obj) else obj

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


def _class_methods(cls):
    """Yield every function object defined directly in `cls.__dict__`."""
    for value in vars(cls).values():
        if isinstance(value, (staticmethod, classmethod)):
            value = value.__func__
        if inspect.isfunction(value):
            yield value


def find_citations(target, _seen=None) -> set[str]:
    """
    Recursively walk the call graph rooted at `target` (a function, bound/
    unbound method, or class) and collect every citation attached with
    :func:`cite` along the way.

    A class is "visited" by inspecting its own methods (including
    ``__init__``) and its base classes. A function/method is "visited" by
    statically finding the names it calls and resolving+recursing into
    whichever of those names refer to functions or classes.

    :param target: The function, method, or class to start the search from.
    :param _seen: Internal set of already-visited objects, used to avoid
        infinite recursion on circular/self-referential call graphs.
    :return: Set of citation identifiers found.
    :rtype: set[str]
    """
    if _seen is None:
        _seen = set()

    key = id(target)
    if key in _seen:
        return set()
    _seen.add(key)

    if inspect.isclass(target):
        citations = set()
        # Citations are collected from anywhere in the MRO, on the class
        # itself and on its methods, since an inherited citation still
        # applies. Method *bodies*, however, are only walked for methods the
        # class defines itself -- recursing into every method of every base
        # class would also pull in unrelated sibling classes (e.g. a base
        # class factory method that references other subclasses), which is
        # not part of what `target` actually uses.
        for klass in target.__mro__:
            if klass is object:
                continue
            citations |= set(_own_citations(klass))
            for method in _class_methods(klass):
                citations |= set(_own_citations(method))
        for method in _class_methods(target):
            citations |= find_citations(method, _seen)
        return citations

    if inspect.isfunction(target) or inspect.ismethod(target):
        citations = set(_own_citations(target))
        for name in _iter_referenced_names(target):
            resolved = _resolve_name(name, target)
            if resolved is None:
                continue
            if inspect.isclass(resolved) or inspect.isfunction(resolved):
                citations |= find_citations(resolved, _seen)
        return citations

    return set()


def load_bibliography(path: Path | None = None) -> dict[str, str]:
    """
    Read the bibliography, which is a plain .bib file in which each entry is
    preceded by a ``%%ID <identifier>`` comment recording the citation it
    answers to::

        %%ID 10.1086/145971
        @ARTICLE{1955ApJ...121..161S,
          ...
        }

    Text outside an entry is ignored by bibtex itself, so the file stays
    usable as an ordinary bibliography. Anything before the first ``%%ID``
    is a header and is skipped. The matching writer lives in
    :mod:`arsenal_gear.utils.refresh_citations`, which generates the file.

    The file is taken as read: it ships with the package and is regenerated
    wholesale, so there is nothing to invalidate and an identifier it does
    not mention is simply missing.

    :param path: Bibliography file to read; defaults to the shipped one.
    :return: Mapping of citation identifier to bibtex entry text.
    :rtype: dict[str, str]
    """
    path = BIBLIOGRAPHY_PATH if path is None else path
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return {}

    entries: dict[str, str] = {}
    current_key = None
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith("%%ID "):
            if current_key is not None:
                entries[current_key] = "\n".join(lines).strip()
            current_key = line[len("%%ID ") :].strip()
            lines = []
        elif current_key is not None:
            lines.append(line)
    if current_key is not None:
        entries[current_key] = "\n".join(lines).strip()

    return {key: entry for key, entry in entries.items() if entry}


def gather_bibtex(target) -> str:
    """
    Find every citation attached to classes/functions used (directly or
    indirectly) by `target`, and return a single string containing bibtex
    entries for all of them.

    Everything is read from the shipped bibliography, in a single pass, so
    this never touches the network. Citations the bibliography is missing
    become bibtex comments rather than raising, so one stale reference does
    not cost you the rest of the bibliography.

    :param target: The function, method, or class to start the search from.
    :return: Newline-separated bibtex entries, one per citation found.
    :rtype: str
    """
    bibliography = load_bibliography()
    entries = []
    for citation in sorted(find_citations(target)):
        entry = bibliography.get(citation)
        entries.append(entry or f"% Could not resolve {citation}")
    return "\n\n".join(entries)
