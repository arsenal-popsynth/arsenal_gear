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

import requests

DOI_BIBTEX_URL = "https://doi.org/{doi}"


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
