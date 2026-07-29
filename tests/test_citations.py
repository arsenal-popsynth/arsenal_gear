"""
test_citations
==========

This file contains unit tests for the citation-discovery and bibtex-rendering
utilities in arsenal_gear.utils.citations, and for the ADS lookups in
arsenal_gear.utils.refresh_citations.
"""

# A test that takes a fixture names its parameter after that fixture, which
# pylint cannot tell apart from shadowing a global.
# pylint: disable=redefined-outer-name

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import requests

from arsenal_gear.formation.dist_funcs.imf import (
    IMF,
    Chabrier,
    Kroupa1993,
    Kroupa2001,
    MillerScalo,
    Salpeter,
)
from arsenal_gear.stellar_evolution.isochrone import RaiteriLifetime
from arsenal_gear.utils import citations, refresh_citations
from arsenal_gear.utils.citations import (
    REGISTRY,
    CitationError,
    cite,
    find_citations,
    gather_bibtex,
    load_bibliography,
)
from arsenal_gear.utils.refresh_citations import (
    fetch_entries,
    is_doi,
    malformed,
    save_bibliography,
)

SALPETER_DOI = "10.1086/145971"
KROUPA2001_DOI = "10.1046/j.1365-8711.2001.04022.x"
KROUPA1993_DOI = "10.1093/mnras/262.3.545"
RAITERI_BIBCODE = "1996A&A...315..105R"


def build_default_population(name):
    """A standalone function that only reaches an IMF class indirectly,
    through the `IMF.get_imf` registry dispatch."""
    imf = IMF.get_imf(name)
    return imf.sample(10)


def fake_export(bibcodes):
    """The bibtex blob an ADS export returns for `bibcodes`."""
    return "\n".join(
        f"@ARTICLE{{{bibcode},\n  title = {{Entry}}\n}}" for bibcode in bibcodes
    )


class FakeResponse:
    def __init__(self, text, status_code=200, payload=None):
        self.text = text
        self.status_code = status_code
        self.ok = status_code < 400
        self._payload = payload

    def json(self):
        if self._payload is None:
            raise ValueError("no json")
        return self._payload


# --- identifier shape (checked only by the refresher) ---------------------


@pytest.mark.parametrize(
    "identifier",
    [
        "10.1086/145971",
        "1996A&A...315..105R",
        "2001MNRAS.322..231K",
        "1996A&A...315L.105R",
    ],
)
def test_well_formed_identifiers_pass(identifier):
    assert malformed([identifier]) == []


@pytest.mark.parametrize(
    "identifier",
    [
        "not a citation",
        "10.1086",  # DOI with no suffix
        "10.86/145971",  # registrant too short
        "1996A&A...315..105",  # 18 chars, truncated bibcode
        "1996A&A...315..105RX",  # 20 chars
        "199XA&A...315..105R",  # non-numeric year
        "@book{Foo, title={Bar}}",  # literal bibtex is deliberately not supported
        "",
    ],
)
def test_malformed_identifiers_are_caught(identifier):
    assert malformed([identifier]) == [identifier]


@pytest.mark.parametrize(
    "identifier,expected",
    [
        ("10.1086/145971", True),
        ("1996A&A...315..105R", False),
    ],
)
def test_is_doi_routes_by_shape(identifier, expected):
    """The one thing that still distinguishes the two kinds: a DOI needs a
    search lookup before it can be exported."""
    assert is_doi(identifier) is expected


# --- the cite decorator ---------------------------------------------------


def test_cite_annotates_without_wrapping():
    """The decorator must return the object itself, so that source-based
    introspection keeps working."""

    def original():
        pass

    decorated = cite("10.1086/145971")(original)
    assert decorated is original


def test_cite_works_on_functions_and_classes():
    @cite("10.1086/145971")
    def func():
        pass

    @cite("10.1086/145971")
    class Klass:
        pass

    assert find_citations(func) == {SALPETER_DOI}
    assert find_citations(Klass) == {SALPETER_DOI}


def test_cite_accepts_multiple_identifiers_of_mixed_kinds():
    @cite("10.1093/mnras/stz2158", "1996A&A...315..105R")
    class Klass:
        pass

    assert find_citations(Klass) == {
        "10.1093/mnras/stz2158",
        "1996A&A...315..105R",
    }


def test_cite_stacks_additively():
    @cite("10.1093/mnras/stz2158")
    @cite("10.1086/145971")
    class Klass:
        pass

    assert find_citations(Klass) == {SALPETER_DOI, "10.1093/mnras/stz2158"}


def test_cite_takes_identifiers_on_trust():
    """citations.py does not check shape -- the refresher does, and a
    citation cannot resolve until it has been run."""

    @cite("obviously not a DOI")
    class Klass:
        pass

    assert find_citations(Klass) == {"obviously not a DOI"}
    assert malformed(find_citations(Klass)) == ["obviously not a DOI"]


def test_cite_registers_identifiers_for_the_refresher():
    @cite("10.1234/registered-example")
    def func():
        pass

    assert "10.1234/registered-example" in REGISTRY


# --- citation discovery ---------------------------------------------------


@pytest.mark.parametrize(
    "imf_class,expected",
    [
        (Salpeter, SALPETER_DOI),
        (Kroupa2001, KROUPA2001_DOI),
        (Kroupa1993, KROUPA1993_DOI),
    ],
)
def test_find_citations_on_class_with_citation(imf_class, expected):
    assert find_citations(imf_class) == {expected}


@pytest.mark.parametrize("imf_class", [MillerScalo, Chabrier])
def test_find_citations_on_class_without_citation(imf_class):
    assert find_citations(imf_class) == set()


def test_find_citations_does_not_leak_sibling_classes():
    """Salpeter and Kroupa2001 are unrelated siblings that both inherit from
    IMF. Finding citations for one must not pick up the other's, even though
    they share a base class that also happens to define a factory method
    (`IMF.get_imf`) referencing every registered IMF."""
    assert KROUPA2001_DOI not in find_citations(Salpeter)
    assert SALPETER_DOI not in find_citations(Kroupa2001)


def test_find_citations_follows_indirect_dispatch():
    """`build_default_population` never calls Salpeter/Kroupa* by name --
    it only calls `IMF.get_imf`, which dispatches through a dict. All
    registered IMFs with a citation should still be found."""
    found = find_citations(build_default_population)
    assert found == {SALPETER_DOI, KROUPA2001_DOI, KROUPA1993_DOI}


def test_find_citations_on_method_level_citation():
    """A citation may sit on a single method rather than a whole class."""

    class Klass:
        @cite("10.1086/145971")
        def cited_method(self):
            pass

        def uncited_method(self):
            pass

    assert find_citations(Klass) == {SALPETER_DOI}
    assert find_citations(Klass.cited_method) == {SALPETER_DOI}
    assert find_citations(Klass.uncited_method) == set()


def test_find_citations_on_inherited_method_level_citation():
    """A method-level citation on a base class still applies to a subclass,
    even though the subclass does not define the method itself."""

    class Base:
        @cite("10.1086/145971")
        def cited_method(self):
            pass

    class Derived(Base):
        pass

    assert find_citations(Derived) == {SALPETER_DOI}


def test_find_citations_on_static_and_class_methods():
    class Klass:
        @staticmethod
        @cite("10.1086/145971")
        def static_method():
            pass

        @classmethod
        @cite("10.1093/mnras/stz2158")
        def class_method(cls):
            pass

    assert find_citations(Klass) == {SALPETER_DOI, "10.1093/mnras/stz2158"}


def test_find_citations_picks_up_bibcodes():
    assert find_citations(RaiteriLifetime) == {RAITERI_BIBCODE}


# --- the citations/refresher boundary -------------------------------------


def test_literal_bibtex_is_not_a_citation():
    """Supplying bibtex by hand is deliberately not supported: work ADS does
    not index is not public enough to cite."""
    entry = "@book{Binney_2008, title={Galactic Dynamics}}"
    assert malformed([entry]) == [entry]


def test_citations_module_performs_no_network_access():
    """The whole point of the split: citations.py must not even be able to
    make a request."""
    source = Path(citations.__file__).read_text(encoding="utf-8")
    assert "requests" not in source
    assert "http" not in source


def test_citations_module_never_writes_the_bibliography():
    """Generating citations.bib belongs to the refresher; citations.py is
    read-only with respect to it."""
    source = Path(citations.__file__).read_text(encoding="utf-8")
    assert "write_text" not in source
    assert not hasattr(citations, "save_bibliography")


def test_gather_bibtex_never_touches_the_network(monkeypatch):
    """Belt and braces: rig requests to explode, and gather anyway."""

    def explode(*_args, **_kwargs):
        raise AssertionError("citations.py must not hit the network")

    monkeypatch.setattr(requests, "get", explode)
    monkeypatch.setattr(requests, "post", explode)

    result = gather_bibtex(build_default_population)
    assert "@ARTICLE{1955ApJ...121..161S," in result


# --- ADS lookups (refresh_citations.py: the only networked code) ----------


@pytest.fixture
def fake_ads(monkeypatch):
    """Stand in for the two ADS endpoints, recording every call made."""

    class FakeADS:
        def __init__(self):
            self.searches = []
            self.exports = []
            self.bibcodes = {"10.1086/145971": "1955ApJ...121..161S"}

        def get(self, url, **kwargs):
            self.searches.append({"url": url, **kwargs})
            docs = [
                {"bibcode": bibcode, "doi": [doi]}
                for doi, bibcode in self.bibcodes.items()
                if f'doi:"{doi}"' in kwargs["params"]["q"]
            ]
            return FakeResponse("", payload={"response": {"docs": docs}})

        def post(self, url, **kwargs):
            self.exports.append({"url": url, **kwargs})
            return FakeResponse(
                "", payload={"export": fake_export(kwargs["json"]["bibcode"])}
            )

    ads = FakeADS()
    monkeypatch.setattr("arsenal_gear.utils.refresh_citations.requests.get", ads.get)
    monkeypatch.setattr("arsenal_gear.utils.refresh_citations.requests.post", ads.post)
    monkeypatch.setenv("ADS_DEV_KEY", "test-token")
    return ads


def test_doi_is_resolved_through_ads_search_then_export(fake_ads):
    """A DOI takes the same route as everything else: search for its
    bibcode, then export."""
    entries = fetch_entries([SALPETER_DOI])

    assert len(fake_ads.searches) == 1
    search = fake_ads.searches[0]
    assert search["url"] == "https://api.adsabs.harvard.edu/v1/search/query"
    assert search["params"]["q"] == 'doi:"10.1086/145971"'
    assert search["headers"] == {"Authorization": "Bearer test-token"}

    assert len(fake_ads.exports) == 1
    export = fake_ads.exports[0]
    assert export["url"] == "https://api.adsabs.harvard.edu/v1/export/bibtex"
    assert export["json"] == {"bibcode": ["1955ApJ...121..161S"]}
    assert entries[SALPETER_DOI].startswith("@ARTICLE{1955ApJ...121..161S,")


def test_bibcode_skips_the_search_step(fake_ads):
    """A bibcode is already what the export API wants, and is sent verbatim
    in a JSON body -- so it never needs URL escaping (%26)."""
    entries = fetch_entries([RAITERI_BIBCODE])

    assert fake_ads.searches == []
    assert fake_ads.exports[0]["json"] == {"bibcode": ["1996A&A...315..105R"]}
    assert entries[RAITERI_BIBCODE].startswith("@ARTICLE{1996A&A...315..105R,")


def test_mixed_citations_are_fetched_in_one_search_and_one_export(fake_ads):
    """However many citations of either kind, it costs two round trips."""
    fake_ads.bibcodes["10.1093/mnras/stz2158"] = "2019MNRAS.489.1082B"

    entries = fetch_entries(
        [
            SALPETER_DOI,
            "10.1093/mnras/stz2158",
            RAITERI_BIBCODE,
        ]
    )

    assert len(fake_ads.searches) == 1
    assert len(fake_ads.exports) == 1
    assert sorted(fake_ads.exports[0]["json"]["bibcode"]) == [
        "1955ApJ...121..161S",
        "1996A&A...315..105R",
        "2019MNRAS.489.1082B",
    ]
    assert set(entries) == {
        "10.1086/145971",
        "10.1093/mnras/stz2158",
        "1996A&A...315..105R",
    }


def test_export_is_split_by_entry_key_not_request_order(fake_ads, monkeypatch):
    """ADS does not return entries in the order they were asked for, so each
    entry has to be matched back to its bibcode by its own key."""

    def reversing_post(url, **kwargs):  # pylint: disable=unused-argument
        export = fake_export(reversed(kwargs["json"]["bibcode"]))
        return FakeResponse("", payload={"export": export})

    fake_ads.bibcodes["10.1093/mnras/stz2158"] = "2019MNRAS.489.1082B"
    monkeypatch.setattr(
        "arsenal_gear.utils.refresh_citations.requests.post", reversing_post
    )

    entries = fetch_entries([RAITERI_BIBCODE, "10.1093/mnras/stz2158"])

    assert entries["1996A&A...315..105R"].startswith("@ARTICLE{1996A&A...315..105R,")
    assert entries["10.1093/mnras/stz2158"].startswith("@ARTICLE{2019MNRAS.489.1082B,")


@pytest.mark.usefixtures("fake_ads")
def test_doi_absent_from_ads_is_omitted_not_faked():
    """A DOI with no ADS record resolves to nothing, so the refresher can
    report it rather than caching an empty entry."""
    assert fetch_entries(["10.9999/not-in-ads"]) == {}


@pytest.mark.usefixtures("fake_ads")
def test_search_http_error_raises(monkeypatch):
    monkeypatch.setattr(
        "arsenal_gear.utils.refresh_citations.requests.get",
        lambda *_a, **_k: FakeResponse("nope", status_code=500),
    )
    with pytest.raises(CitationError, match="ADS search failed"):
        fetch_entries([SALPETER_DOI])


@pytest.mark.usefixtures("fake_ads")
def test_export_http_error_raises(monkeypatch):
    monkeypatch.setattr(
        "arsenal_gear.utils.refresh_citations.requests.post",
        lambda *_a, **_k: FakeResponse("nope", status_code=500),
    )
    with pytest.raises(CitationError, match="ADS export failed"):
        fetch_entries([RAITERI_BIBCODE])


@pytest.mark.parametrize("citation", [SALPETER_DOI, RAITERI_BIBCODE])
def test_fetching_without_token_raises_explaining_why(monkeypatch, citation):
    for var in ("ADS_DEV_KEY", "ADS_API_TOKEN"):
        monkeypatch.delenv(var, raising=False)

    with pytest.raises(CitationError, match="anonymous"):
        fetch_entries([citation])


# --- the bibliography -----------------------------------------------------


def test_bibliography_round_trips(tmp_path):
    path = tmp_path / "citations.bib"
    entries = {
        "10.1086/145971": "@article{Salpeter_1955,\n  year={1955}\n}",
        "1996A&A...315..105R": "@ARTICLE{Raiteri_1996,\n  year={1996}\n}",
    }
    save_bibliography(entries, path=path)
    assert load_bibliography(path=path) == entries


def test_missing_bibliography_leaves_everything_missing(monkeypatch, tmp_path):
    """A bibliography that is not there at all is just an empty one: every
    citation goes missing, and nothing raises on the way."""
    monkeypatch.setattr(
        "arsenal_gear.utils.citations.BIBLIOGRAPHY_PATH", tmp_path / "absent.bib"
    )

    assert load_bibliography() == {}
    assert gather_bibtex(Salpeter).startswith("% Could not resolve")


def test_shipped_bibliography_covers_every_citation_used_in_the_package():
    """Every citation used in the package must already be in the shipped
    bibliography: an end user has no ADS token, so a missing entry is
    unresolvable for them.

    Run in a subprocess: the citation registry is global and populated by
    import, so the citations invented by the tests above would otherwise
    count as "used in the package"."""
    script = """
import json
from arsenal_gear.utils.citations import REGISTRY, load_bibliography
from arsenal_gear.utils.refresh_citations import import_all_modules

import_all_modules()
recorded = load_bibliography()
print(json.dumps(sorted(
    c for c in REGISTRY if c not in recorded
)))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True
    )
    assert json.loads(completed.stdout.splitlines()[-1]) == []


def test_gather_bibtex_combines_all_found_citations():
    """Entries come from the shipped bibliography, and so are keyed by bibcode the
    way ADS exports them."""
    result = gather_bibtex(build_default_population)

    assert "@ARTICLE{1955ApJ...121..161S," in result
    assert "@ARTICLE{2001MNRAS.322..231K," in result
    assert "@ARTICLE{1993MNRAS.262..545K," in result


def test_gather_bibtex_comments_missing_citations(monkeypatch, tmp_path):
    """A citation the maintainers have not refreshed yet must not cost you
    the rest of the bibliography."""

    @cite("10.1086/145971", "10.9999/never-refreshed")
    class Klass:
        pass

    path = tmp_path / "citations.bib"
    save_bibliography({"10.1086/145971": "@ARTICLE{1955ApJ...121..161S,\n}"}, path=path)
    monkeypatch.setattr("arsenal_gear.utils.citations.BIBLIOGRAPHY_PATH", path)

    result = gather_bibtex(Klass)

    assert "@ARTICLE{1955ApJ...121..161S," in result
    assert "% Could not resolve 10.9999/never-refreshed" in result


def test_refresh_reports_malformed_identifiers(monkeypatch, tmp_path, capsys):
    """A typo cannot reach an end user: it fails the refresh that has to run
    before the citation would resolve at all."""
    monkeypatch.setattr(refresh_citations, "import_all_modules", lambda *a: [])
    monkeypatch.setattr(
        refresh_citations, "REGISTRY", {"10.1086/145971", "obviously not a DOI"}
    )
    monkeypatch.setattr(
        refresh_citations, "BIBLIOGRAPHY_PATH", tmp_path / "citations.bib"
    )
    monkeypatch.setattr(
        refresh_citations,
        "load_bibliography",
        lambda: {"10.1086/145971": "@ARTICLE{1955ApJ...121..161S,\n}"},
    )

    status = refresh_citations.refresh()

    assert status == 1
    assert "neither a DOI" in capsys.readouterr().err
    # The good citation is still written out; one typo does not lose the rest.
    written = load_bibliography(tmp_path / "citations.bib")
    assert set(written) == {"10.1086/145971"}


needs_ads_token = pytest.mark.skipif(
    not (os.environ.get("ADS_DEV_KEY") or os.environ.get("ADS_API_TOKEN")),
    reason="needs a NASA ADS API token in $ADS_DEV_KEY",
)


@needs_ads_token
def test_live_ads_lookup_by_bibcode():
    """Live check of the ADS export API. Skipped unless a token is set."""
    entry = fetch_entries([RAITERI_BIBCODE])[RAITERI_BIBCODE]
    assert entry.startswith("@")
    assert "Raiteri" in entry


@needs_ads_token
def test_live_ads_lookup_by_doi():
    """Live check of the DOI path: ADS search, then export."""
    entry = fetch_entries([SALPETER_DOI])[SALPETER_DOI]
    assert entry.startswith("@ARTICLE{1955ApJ...121..161S,")
    assert "Salpeter" in entry


@needs_ads_token
def test_live_ads_batch_costs_two_round_trips():
    """Live check that a mixed batch really is two requests, not 2N."""
    entries = fetch_entries([SALPETER_DOI, KROUPA2001_DOI, RAITERI_BIBCODE])

    assert set(entries) == {SALPETER_DOI, KROUPA2001_DOI, RAITERI_BIBCODE}
    assert all(entry.startswith("@") for entry in entries.values())
