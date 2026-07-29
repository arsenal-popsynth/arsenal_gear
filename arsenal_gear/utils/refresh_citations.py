"""
refresh_citations.py
====================

Maintainer tool for regenerating the bibliography that ships with
arsenal_gear as static data. This is the only part of the package that
performs network access; :mod:`arsenal_gear.utils.citations` only reads the
bibliography.

Run it after adding or changing a :func:`~arsenal_gear.utils.citations.cite`
decorator anywhere in the package, and commit the result::

    python -m arsenal_gear.utils.refresh_citations

It imports every arsenal_gear module (which is what registers the citations),
looks up bibtex for any identifier the bibliography is missing, and rewrites
``arsenal_gear/utils/citations.bib``. End users then get every citation with
no network connection and no credentials.

Everything is fetched from NASA ADS, whatever kind of identifier it started
as, in two requests total: one search to turn DOIs into bibcodes, then one
export. ADS has no anonymous API, so this needs a (free) ADS API token in
``$ADS_DEV_KEY``. Without one nothing can be fetched, though entries already
in the bibliography are kept.
"""

from __future__ import annotations

import argparse
import importlib
import os
import pkgutil
import re
import sys
from pathlib import Path

import requests

import arsenal_gear

from .citations import (
    BIBLIOGRAPHY_PATH,
    REGISTRY,
    CitationError,
    load_bibliography,
)

ADS_SEARCH_URL = "https://api.adsabs.harvard.edu/v1/search/query"
ADS_EXPORT_URL = "https://api.adsabs.harvard.edu/v1/export/bibtex"
ADS_TOKEN_VARS = ("ADS_DEV_KEY", "ADS_API_TOKEN")

#: Identifiers per request. ADS accepts far more, but chunking keeps search
#: query strings to a sane length and the failure blast radius small.
ADS_CHUNK_SIZE = 50

# The start of a bibtex entry, capturing its key: the "Foo_1955" of
# "@article{Foo_1955, ...}". ADS keys its entries by bibcode, which is what
# lets a batched export be split back up per citation.
BIBTEX_ENTRY_RE = re.compile(r"^@\w+\s*\{\s*([^,\s]+)\s*,", re.MULTILINE)

# A DOI is the prefix "10.", a registrant code of at least four digits, a
# slash, and a non-empty suffix.
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")

# An ADS bibcode is exactly 19 characters: a four digit year, a five
# character journal abbreviation, a four character volume, a one character
# qualifier, a four character page, and the first author's initial. Unused
# fields are padded with dots.
BIBCODE_RE = re.compile(
    r"^\d{4}[A-Za-z0-9&.]{5}[A-Za-z0-9.]{4}[A-Za-z.][A-Za-z0-9.]{4}[A-Za-z.]$"
)


def _ads_headers() -> dict[str, str]:
    """
    Authorization header for the ADS API.

    :raises CitationError: If no token is set in the environment.
    """
    token = next(
        (os.environ[var] for var in ADS_TOKEN_VARS if os.environ.get(var)), None
    )
    if token is None:
        raise CitationError(
            f"NASA ADS has no anonymous API. Set ${ADS_TOKEN_VARS[0]} to an ADS "
            f"API token to refresh {BIBLIOGRAPHY_PATH.name}."
        )
    return {"Authorization": f"Bearer {token}"}


def _chunks(items: list, size: int = ADS_CHUNK_SIZE):
    """Split `items` into consecutive lists of at most `size` elements."""
    for start in range(0, len(items), size):
        yield items[start : start + size]


def _split_entries(text: str) -> dict[str, str]:
    """
    Split a run of concatenated bibtex entries into a mapping of entry key to
    entry text. ADS keys its entries by bibcode, and does not necessarily
    return them in the order they were requested, so keying by the entry
    itself is what makes a batched export usable.

    :param text: One or more bibtex entries.
    :return: Mapping of bibtex entry key to entry text.
    :rtype: dict[str, str]
    """
    matches = list(BIBTEX_ENTRY_RE.finditer(text))
    entries = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        entries[match.group(1)] = text[match.start() : end].strip()
    return entries


def resolve_bibcodes(dois, timeout: float = 10) -> dict[str, str]:
    """
    Look up the ADS bibcode for each of `dois`, batching the whole lot into
    as few search queries as possible.

    DOIs with no ADS record are simply absent from the result rather than
    raising, so one unindexed paper does not sink the rest of the batch.

    :param dois: Iterable of DOI strings.
    :param timeout: Per-request timeout in seconds.
    :return: Mapping of DOI to bibcode, for those that resolved.
    :rtype: dict[str, str]
    :raises CitationError: If no ADS token is set, or a query fails.
    """
    wanted = list(dict.fromkeys(dois))
    resolved: dict[str, str] = {}

    for chunk in _chunks(wanted):
        query = " OR ".join(f'doi:"{doi}"' for doi in chunk)
        try:
            response = requests.get(
                ADS_SEARCH_URL,
                params={"q": query, "fl": "bibcode,doi", "rows": len(chunk) + 10},
                headers=_ads_headers(),
                timeout=timeout,
            )
        except requests.exceptions.RequestException as e:
            raise CitationError(f"ADS search failed: {e}") from e
        if not response.ok:
            raise CitationError(f"ADS search failed: HTTP {response.status_code}")

        try:
            docs = response.json()["response"]["docs"]
        except (ValueError, KeyError) as e:
            raise CitationError(f"Unexpected ADS search response: {e}") from e

        requested = set(chunk)
        for doc in docs:
            # A record may list several DOIs (e.g. a preprint alongside the
            # published version); only the ones actually asked for count.
            for doi in doc.get("doi", []):
                if doi in requested:
                    resolved.setdefault(doi, doc["bibcode"])

    return resolved


def export_bibtex(bibcodes, timeout: float = 10) -> dict[str, str]:
    """
    Export bibtex for each of `bibcodes`, batching the whole lot into as few
    export requests as possible.

    :param bibcodes: Iterable of ADS bibcode strings.
    :param timeout: Per-request timeout in seconds.
    :return: Mapping of bibcode to bibtex entry text, for those that
        resolved.
    :rtype: dict[str, str]
    :raises CitationError: If no ADS token is set, or an export fails.
    """
    wanted = list(dict.fromkeys(bibcodes))
    entries: dict[str, str] = {}

    for chunk in _chunks(wanted):
        try:
            response = requests.post(
                ADS_EXPORT_URL,
                json={"bibcode": chunk},
                headers=_ads_headers(),
                timeout=timeout,
            )
        except requests.exceptions.RequestException as e:
            raise CitationError(f"ADS export failed: {e}") from e
        if not response.ok:
            raise CitationError(f"ADS export failed: HTTP {response.status_code}")

        try:
            export = response.json()["export"]
        except (ValueError, KeyError) as e:
            raise CitationError(f"Unexpected ADS export response: {e}") from e
        entries.update(_split_entries(export))

    return entries


def is_doi(identifier: str) -> bool:
    """
    Whether `identifier` is a DOI rather than an ADS bibcode. This is the one
    place the difference matters: the export API speaks bibcodes, so a DOI
    needs a search lookup first. Everything in
    :mod:`arsenal_gear.utils.citations` treats identifiers as opaque.

    :param identifier: A well formed citation identifier.
    :return: True for a DOI, False for a bibcode.
    :rtype: bool
    """
    return bool(DOI_RE.match(identifier))


def malformed(identifiers) -> list[str]:
    """
    Pick out identifiers that are neither a DOI nor an ADS bibcode.

    Shape is checked here rather than in
    :mod:`arsenal_gear.utils.citations` because this is the only code that
    can do anything about it: a citation has to be looked up before it
    resolves, so a typo cannot escape a refresh.

    :param identifiers: Iterable of citation identifiers.
    :return: Those matching neither form, sorted.
    :rtype: list[str]
    """
    return sorted(
        i for i in identifiers if not (DOI_RE.match(i) or BIBCODE_RE.match(i))
    )


def fetch_entries(citations, timeout: float = 10) -> dict[str, str]:
    """
    Fetch bibtex for a collection of citation identifiers, of either kind, in
    as few requests as possible: one search pass to turn DOIs into bibcodes,
    then one export pass for every bibcode involved.

    Citations that cannot be resolved are absent from the result rather than
    raising, so callers can report them individually.

    :param citations: Iterable of citation identifiers.
    :param timeout: Per-request timeout in seconds.
    :return: Mapping of citation identifier to bibtex entry text.
    :rtype: dict[str, str]
    :raises CitationError: If no ADS token is set, or a request fails.
    """
    citations = list(citations)
    resolved = resolve_bibcodes([c for c in citations if is_doi(c)], timeout=timeout)

    # Several citations can map onto one bibcode (a DOI and the bibcode for
    # the same paper), so track every identifier each bibcode answers for.
    ids_by_bibcode: dict[str, list[str]] = {}
    for citation in citations:
        bibcode = resolved.get(citation) if is_doi(citation) else citation
        if bibcode is not None:
            ids_by_bibcode.setdefault(bibcode, []).append(citation)

    exported = export_bibtex(ids_by_bibcode, timeout=timeout)
    return {
        identifier: exported[bibcode]
        for bibcode, identifiers in ids_by_bibcode.items()
        if bibcode in exported
        for identifier in identifiers
    }


def save_bibliography(entries: dict[str, str], path: Path | None = None) -> None:
    """
    Write the bibliography: every entry, sorted by identifier for a stable
    diff, each preceded by the ``%%ID <identifier>`` comment that
    :func:`~arsenal_gear.utils.citations.load_bibliography` reads it back by.

    :param entries: Mapping of citation identifier to bibtex entry text.
    :param path: Bibliography file to write; defaults to the shipped bibliography.
    """
    path = BIBLIOGRAPHY_PATH if path is None else path
    chunks = [
        "%% arsenal_gear bibliography -- generated, do not edit by hand.",
        "%% Refresh with: python -m arsenal_gear.utils.refresh_citations",
        "",
    ]
    for key in sorted(entries):
        chunks.append(f"%%ID {key}")
        chunks.append(entries[key].strip())
        chunks.append("")
    path.write_text("\n".join(chunks), encoding="utf-8")


def import_all_modules() -> list[str]:
    """
    Import every arsenal_gear submodule, so that the decorators throughout it
    run and register their citations.

    :return: Names of modules that could not be imported.
    :rtype: list[str]
    """
    failed = []
    for module in pkgutil.walk_packages(
        arsenal_gear.__path__, f"{arsenal_gear.__name__}."
    ):
        try:
            importlib.import_module(module.name)
        # A module that will not import must not stop the sweep: the rest of
        # the package still has citations to register.
        except Exception as e:  # pylint: disable=broad-exception-caught
            failed.append(f"{module.name}: {e}")
    return failed


def refresh(force: bool = False, timeout: float = 10) -> int:
    """
    Regenerate the bibliography from the citations registered across the
    package.

    :param force: Re-fetch even identifiers already in the bibliography.
    :param timeout: Per-request timeout in seconds.
    :return: Process exit status: 0 if every citation resolved, 1 otherwise.
    :rtype: int
    """
    for failure in import_all_modules():
        print(f"warning: could not import {failure}", file=sys.stderr)

    citations = sorted(REGISTRY)
    if not citations:
        print("No citations registered; nothing to do.", file=sys.stderr)
        return 1

    # citations.py takes identifiers on trust, so this is where a typo is
    # caught. Bad ones are reported and set aside rather than being sent to
    # ADS, where they would come back as an unhelpful "no record found".
    bad = malformed(citations)
    rejected = set(bad)
    citations = [c for c in citations if c not in rejected]

    bibliography = {} if force else load_bibliography()
    missing = [c for c in citations if c not in bibliography]
    kept = len(citations) - len(missing)

    # One batch for the whole package: a request-level failure (no token, ADS
    # unreachable) takes out every citation at once and is worth reporting
    # only once, whereas individual absences are reported per citation.
    fetched: dict[str, str] = {}
    batch_error = None
    if missing:
        try:
            fetched = fetch_entries(missing, timeout=timeout)
        except CitationError as e:
            batch_error = str(e)

    failures = [
        f"{c!r} is neither a DOI (10.xxxx/...) nor a 19 character ADS bibcode "
        f"(e.g. 1996A&A...315..105R)"
        for c in bad
    ]
    if batch_error:
        failures.append(batch_error)
    for citation in missing:
        if citation in fetched:
            bibliography[citation] = fetched[citation]
            print(f"fetched {citation}")
        elif batch_error is None:
            failures.append(f"No NASA ADS record found for {citation}")

    # Identifiers that no longer appear anywhere in the package are dropped,
    # so the bibliography cannot accumulate entries for deleted citations.
    dropped = sorted(set(bibliography) - set(citations))
    for key in dropped:
        del bibliography[key]
        print(f"dropped stale {key}")

    save_bibliography(bibliography)
    print(
        f"\n{BIBLIOGRAPHY_PATH}: {len(bibliography)} entries "
        f"({len(fetched)} fetched, {kept} already present, {len(dropped)} dropped)"
    )

    for failure in failures:
        print(f"FAILED: {failure}", file=sys.stderr)
    return 1 if failures else 0


def main(argv=None) -> int:
    """Command line entry point."""
    parser = argparse.ArgumentParser(
        description="Regenerate the bibliography shipped with arsenal_gear."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="re-fetch every citation instead of only the missing ones",
    )
    parser.add_argument(
        "--timeout", type=float, default=10, help="per-request timeout in seconds"
    )
    args = parser.parse_args(argv)
    return refresh(force=args.force, timeout=args.timeout)


if __name__ == "__main__":
    sys.exit(main())
