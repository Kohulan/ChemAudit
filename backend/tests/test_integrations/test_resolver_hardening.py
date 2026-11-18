"""
Regression tests for the /resolve hardening (security scan findings F2 and F23).
"""

import logging
import time
from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from app.services.integrations.wikidata import WikidataClient


def test_extract_title_basic():
    c = WikidataClient()
    assert c._extract_title("https://en.wikipedia.org/wiki/Aspirin#History") == "Aspirin"
    assert c._extract_title("https://en.wikipedia.org/wiki/Acetic_acid") == "Acetic acid"
    assert c._extract_title("en.wikipedia.org/wiki/Caf%C3%A9") == "Café"
    assert c._extract_title("https://example.org/Aspirin") is None
    assert c._extract_title("https://en.wikipedia.org/wiki/") is None


def test_extract_title_pathological_input_is_fast():
    start = time.perf_counter()
    WikidataClient()._extract_title("wikipedia.org/wiki/" + "#" * 20_000 + "\nx")
    assert time.perf_counter() - start < 0.5


async def test_resolve_rejects_long_identifier(client):
    r = await client.post("/api/v1/resolve", json={"identifier": "x" * 2049})
    assert r.status_code == 422


async def test_resolver_logs_identifier_with_newlines_escaped(caplog):
    from app.services.integrations import resolver

    async def boom(_identifier):
        raise RuntimeError("upstream failed")

    with patch.dict(resolver._RESOLVERS, {"name": boom}), caplog.at_level(logging.ERROR):
        await resolver.resolve_identifier("evil\nFAKE LOG LINE", "name")
    assert "evil\nFAKE" not in caplog.text
    assert "evil\\nFAKE" in caplog.text


def _capturing_client():
    """Patch target for httpx.AsyncClient that records requested URLs."""
    calls = []
    response = MagicMock()
    response.raise_for_status.side_effect = RuntimeError("stop")
    inner = MagicMock()
    inner.get = AsyncMock(side_effect=lambda url, **kw: calls.append(url) or response)
    inner.post = AsyncMock(side_effect=lambda url, **kw: calls.append(url) or response)
    ctx = MagicMock()
    ctx.__aenter__ = AsyncMock(return_value=inner)
    ctx.__aexit__ = AsyncMock(return_value=False)
    return MagicMock(return_value=ctx), calls


@pytest.mark.parametrize(
    "module,cls,method",
    [
        ("pubchem", "PubChemClient", "search_by_name"),
        ("pubchem", "PubChemClient", "search_by_inchikey"),
        ("chembl", "ChEMBLClient", "get_molecule"),
        ("unichem", "UniChemClient", "get_cross_references"),
        ("surechembl", "SureChEMBLClient", "_get_schembl_id_via_unichem"),
        ("surechembl", "SureChEMBLClient", "_get_compound_data"),
    ],
)
async def test_path_segments_are_url_encoded(module, cls, method):
    import importlib

    mod = importlib.import_module(f"app.services.integrations.{module}")
    factory, calls = _capturing_client()
    with patch.object(mod.httpx, "AsyncClient", factory):
        try:
            await getattr(getattr(mod, cls)(), method)("../../x y")
        except Exception:
            pass
    assert calls, "no request made"
    assert "../" not in calls[0] and "%2F" in calls[0]
