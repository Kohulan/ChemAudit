"""Redis clients must be pooled singletons, not created per call."""

import app.core.rate_limit as rl


def test_sync_redis_client_is_singleton():
    rl._SYNC_REDIS = None  # reset module cache
    c1 = rl.get_sync_redis_client()
    c2 = rl.get_sync_redis_client()
    assert c1 is c2


async def test_async_redis_client_is_singleton():
    import app.core.security as sec

    sec._ASYNC_REDIS = None
    c1 = await sec.get_redis_client()
    c2 = await sec.get_redis_client()
    assert c1 is c2
