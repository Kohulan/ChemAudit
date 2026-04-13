"""Tests for MCP server configuration and safety.

Validates the MCP tag allowlist, admin endpoint exclusion,
startup safety assertion, and MCP endpoint mount.
"""


from app.main import _MCP_EXCLUDED_TAGS, MCP_INCLUDE_TAGS

# ---------------------------------------------------------------------------
# Tag count and content
# ---------------------------------------------------------------------------


def test_mcp_include_tags_count():
    """MCP_INCLUDE_TAGS must contain exactly 16 tags."""
    assert len(MCP_INCLUDE_TAGS) == 16, (
        f"Expected 16 MCP include tags, got {len(MCP_INCLUDE_TAGS)}: {MCP_INCLUDE_TAGS}"
    )


def test_mcp_include_tags_contains_required():
    """All required MCP tags must be present in the include list."""
    required_tags = [
        "health",
        "validation",
        "scoring",
        "standardization",
        "alerts",
        "profiler",
        "safety",
        "diagnostics",
        "qsar-ready",
        "structure-filter",
        "dataset-intelligence",
        "batch",
        "export",
        "integrations",
        "resolve",
        "profiles",
    ]
    for tag in required_tags:
        assert tag in MCP_INCLUDE_TAGS, f"Required tag '{tag}' missing from MCP_INCLUDE_TAGS"


# ---------------------------------------------------------------------------
# Admin tag exclusion (D-02 / D-03)
# ---------------------------------------------------------------------------


def test_mcp_excludes_admin_tags():
    """Admin tags must NOT appear in MCP_INCLUDE_TAGS."""
    admin_tags = ["api-keys", "config", "session", "bookmarks", "permalinks", "history"]
    for tag in admin_tags:
        assert tag not in MCP_INCLUDE_TAGS, (
            f"Admin tag '{tag}' must not be in MCP_INCLUDE_TAGS"
        )


def test_admin_exclusion_assertion():
    """D-03 safety check: intersection of include and excluded tags must be empty."""
    leaked = set(MCP_INCLUDE_TAGS) & _MCP_EXCLUDED_TAGS
    assert leaked == set(), f"Admin tags leaked into MCP include list: {leaked}"


# ---------------------------------------------------------------------------
# Quality checks
# ---------------------------------------------------------------------------


def test_mcp_no_duplicate_tags():
    """MCP_INCLUDE_TAGS must not contain duplicate entries."""
    assert len(MCP_INCLUDE_TAGS) == len(set(MCP_INCLUDE_TAGS)), (
        "Duplicate tags found in MCP_INCLUDE_TAGS"
    )


# ---------------------------------------------------------------------------
# MCP endpoint mount
# ---------------------------------------------------------------------------


def test_mcp_mount_endpoint_exists():
    """MCP endpoint at /mcp must be registered as a route on the app.

    The SSE-based MCP endpoint streams data, so we verify it exists by
    inspecting FastAPI's registered routes rather than making an HTTP
    request that would hang waiting for SSE events.
    """
    from app.main import app

    route_paths = [getattr(r, "path", "") for r in app.routes]
    # fastapi-mcp registers GET /mcp for SSE connections and POST /mcp/messages
    assert "/mcp" in route_paths, (
        f"MCP endpoint /mcp not found in app routes. Available: {route_paths}"
    )


# ---------------------------------------------------------------------------
# OpenAPI tag cross-check (D-03 full assertion)
# ---------------------------------------------------------------------------


def test_openapi_admin_endpoints_not_in_mcp_tags():
    """Admin endpoint paths must not carry tags that overlap with MCP include tags.

    This verifies the full D-03 assertion against the real OpenAPI schema:
    no path tagged with an MCP-included tag should be an admin-only endpoint.
    """
    from app.main import app

    openapi = app.openapi()
    mcp_tag_set = set(MCP_INCLUDE_TAGS)

    # Known admin paths that must never be exposed via MCP
    admin_paths = {"/api/v1/keys", "/api/v1/admin/config"}

    for path, methods in openapi.get("paths", {}).items():
        for method, operation in methods.items():
            if method in ("get", "post", "put", "delete", "patch"):
                operation_tags = set(operation.get("tags", []))
                if operation_tags & mcp_tag_set:
                    assert path not in admin_paths, (
                        f"Admin endpoint {path} has MCP-included tag(s): "
                        f"{operation_tags & mcp_tag_set}"
                    )
