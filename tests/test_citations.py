"""
test_citations
==========

This file contains unit tests for the DOI-discovery and bibtex-fetching
utilities in arsenal_gear.utils.citations.
"""

import pytest

from arsenal_gear.formation.dist_funcs.imf import (
    IMF,
    Chabrier,
    Kroupa1993,
    Kroupa2001,
    MillerScalo,
    Salpeter,
)
from arsenal_gear.utils.citations import doi_to_bibtex, find_dois, gather_bibtex

SALPETER_DOI = "10.1086/145971"
KROUPA2001_DOI = "10.1046/j.1365-8711.2001.04022.x"
KROUPA1993_DOI = "10.1093/mnras/262.3.545"


def build_default_population(name):
    """A standalone function that only reaches an IMF class indirectly,
    through the `IMF.get_imf` registry dispatch."""
    imf = IMF.get_imf(name)
    return imf.sample(10)


class FakeResponse:
    def __init__(self, text, status_code=200):
        self.text = text
        self.status_code = status_code
        self.ok = status_code < 400


@pytest.mark.parametrize(
    "imf_class,expected_doi",
    [
        (Salpeter, SALPETER_DOI),
        (Kroupa2001, KROUPA2001_DOI),
        (Kroupa1993, KROUPA1993_DOI),
    ],
)
def test_find_dois_on_class_with_doi(imf_class, expected_doi):
    assert find_dois(imf_class) == {expected_doi}


@pytest.mark.parametrize("imf_class", [MillerScalo, Chabrier])
def test_find_dois_on_class_without_doi(imf_class):
    assert find_dois(imf_class) == set()


def test_find_dois_does_not_leak_sibling_classes():
    """Salpeter and Kroupa2001 are unrelated siblings that both inherit from
    IMF. Finding DOIs for one must not pick up the other's, even though they
    share a base class that also happens to define a factory method
    (`IMF.get_imf`) referencing every registered IMF."""
    assert KROUPA2001_DOI not in find_dois(Salpeter)
    assert SALPETER_DOI not in find_dois(Kroupa2001)


def test_find_dois_follows_indirect_dispatch():
    """`build_default_population` never calls Salpeter/Kroupa* by name --
    it only calls `IMF.get_imf`, which dispatches through a dict. All
    registered IMFs with a DOI should still be found."""
    dois = find_dois(build_default_population)
    assert dois == {SALPETER_DOI, KROUPA2001_DOI, KROUPA1993_DOI}


def test_doi_to_bibtex_requests_correct_url_and_headers(monkeypatch):
    captured = {}

    def fake_get(url, headers=None, **_kwargs):
        captured["url"] = url
        captured["headers"] = headers
        return FakeResponse("@article{Salpeter_1955, ...}")

    monkeypatch.setattr("arsenal_gear.utils.citations.requests.get", fake_get)

    result = doi_to_bibtex(SALPETER_DOI)

    assert captured["url"] == f"https://doi.org/{SALPETER_DOI}"
    assert captured["headers"] == {"Accept": "application/x-bibtex"}
    assert result == "@article{Salpeter_1955, ...}"


def test_doi_to_bibtex_raises_on_http_error(monkeypatch):
    def fake_get(*_args, **_kwargs):
        return FakeResponse("not found", status_code=404)

    monkeypatch.setattr("arsenal_gear.utils.citations.requests.get", fake_get)

    with pytest.raises(RuntimeError):
        doi_to_bibtex(SALPETER_DOI)


def test_gather_bibtex_combines_all_found_dois(monkeypatch):
    def fake_get(url, **_kwargs):
        doi = url.removeprefix("https://doi.org/")
        return FakeResponse(f"@article{{{doi}}}")

    monkeypatch.setattr("arsenal_gear.utils.citations.requests.get", fake_get)

    result = gather_bibtex(Salpeter)

    assert f"@article{{{SALPETER_DOI}}}" in result
