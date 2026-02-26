"""
Scoring Profile Service

CRUD operations for custom scoring profiles with 8 preset templates.
Presets are immutable (cannot be updated or deleted) but can be duplicated.
"""

import json
import logging
from typing import List, Optional

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.db.models.profile import ScoringProfile

logger = logging.getLogger(__name__)

# =============================================================================
# 8 Preset Profile Definitions
# =============================================================================

PRESET_PROFILES = [
    {
        "name": "Drug-like (Lipinski)",
        "description": "Lipinski Rule of Five: MW<500, LogP<5, HBD<=5, HBA<=10",
        "thresholds": {
            "mw": {"min": 0, "max": 500},
            "logp": {"min": -5, "max": 5},
            "hbd": {"min": 0, "max": 5},
            "hba": {"min": 0, "max": 10},
        },
        "weights": {"mw": 1.0, "logp": 1.0, "hbd": 1.0, "hba": 1.0},
    },
    {
        "name": "Lead-like",
        "description": "Lead-like space: MW 200-450, LogP -1..4, RotBonds<=7",
        "thresholds": {
            "mw": {"min": 200, "max": 450},
            "logp": {"min": -1, "max": 4},
            "rotatable_bonds": {"min": 0, "max": 7},
        },
        "weights": {"mw": 1.0, "logp": 1.0, "rotatable_bonds": 1.0},
    },
    {
        "name": "Fragment-like (Rule of 3)",
        "description": "Fragment screening: MW<300, LogP<=3, HBD<=3, HBA<=3, RotBonds<=3",
        "thresholds": {
            "mw": {"min": 0, "max": 300},
            "logp": {"min": -3, "max": 3},
            "hbd": {"min": 0, "max": 3},
            "hba": {"min": 0, "max": 3},
            "rotatable_bonds": {"min": 0, "max": 3},
        },
        "weights": {
            "mw": 1.0,
            "logp": 1.0,
            "hbd": 1.0,
            "hba": 1.0,
            "rotatable_bonds": 1.0,
        },
    },
    {
        "name": "CNS-penetrant",
        "description": "Blood-brain barrier permeable: MW<400, LogP 1-5, TPSA<90, HBD<=2",
        "thresholds": {
            "mw": {"min": 0, "max": 400},
            "logp": {"min": 1, "max": 5},
            "tpsa": {"min": 0, "max": 90},
            "hbd": {"min": 0, "max": 2},
        },
        "weights": {"mw": 1.0, "logp": 1.0, "tpsa": 1.0, "hbd": 1.0},
    },
    {
        "name": "Ghose (Amgen)",
        "description": "Ghose filter: MW 160-480, LogP -0.4..5.6, atoms 20-70, refractivity 40-130",
        "thresholds": {
            "mw": {"min": 160, "max": 480},
            "logp": {"min": -0.4, "max": 5.6},
            "num_atoms": {"min": 20, "max": 70},
            "refractivity": {"min": 40, "max": 130},
        },
        "weights": {"mw": 1.0, "logp": 1.0, "num_atoms": 1.0, "refractivity": 1.0},
    },
    {
        "name": "Veber (GSK)",
        "description": "Veber rules for oral bioavailability: RotBonds<=10, TPSA<=140",
        "thresholds": {
            "rotatable_bonds": {"min": 0, "max": 10},
            "tpsa": {"min": 0, "max": 140},
        },
        "weights": {"rotatable_bonds": 1.0, "tpsa": 1.0},
    },
    {
        "name": "PPI-like",
        "description": "Protein-protein interaction space: MW 400-800, LogP 2-6, RotBonds 5-15",
        "thresholds": {
            "mw": {"min": 400, "max": 800},
            "logp": {"min": 2, "max": 6},
            "rotatable_bonds": {"min": 5, "max": 15},
        },
        "weights": {"mw": 1.0, "logp": 1.0, "rotatable_bonds": 1.0},
    },
    {
        "name": "NP-like",
        "description": "Natural product space: MW 200-800, RotBonds<=10, HBD<=5, rings>=2",
        "thresholds": {
            "mw": {"min": 200, "max": 800},
            "rotatable_bonds": {"min": 0, "max": 10},
            "hbd": {"min": 0, "max": 5},
            "rings": {"min": 2, "max": 20},
        },
        "weights": {"mw": 1.0, "rotatable_bonds": 1.0, "hbd": 1.0, "rings": 1.0},
    },
]


class ProfileService:
    """Service for scoring profile CRUD and preset management."""

    async def create(
        self,
        db: AsyncSession,
        name: str,
        description: Optional[str],
        thresholds: dict,
        weights: dict,
    ) -> ScoringProfile:
        """Create a new user profile (is_preset=False)."""
        profile = ScoringProfile(
            name=name,
            description=description,
            thresholds=json.dumps(thresholds),
            weights=json.dumps(weights),
            is_preset=False,
            is_active=True,
        )
        db.add(profile)
        await db.commit()
        await db.refresh(profile)
        return profile

    async def get(self, db: AsyncSession, profile_id: int) -> Optional[ScoringProfile]:
        """Get a profile by ID."""
        result = await db.execute(
            select(ScoringProfile).where(
                ScoringProfile.id == profile_id,
                ScoringProfile.is_active.is_(True),
            )
        )
        return result.scalar_one_or_none()

    async def list_all(
        self, db: AsyncSession, include_presets: bool = True
    ) -> List[ScoringProfile]:
        """List all active profiles."""
        query = select(ScoringProfile).where(ScoringProfile.is_active.is_(True))
        if not include_presets:
            query = query.where(ScoringProfile.is_preset.is_(False))
        query = query.order_by(ScoringProfile.is_preset.desc(), ScoringProfile.name)
        result = await db.execute(query)
        return list(result.scalars().all())

    async def update(
        self,
        db: AsyncSession,
        profile_id: int,
        updates: dict,
    ) -> Optional[ScoringProfile]:
        """
        Update a user profile. Refuses to update presets.

        Args:
            db: Database session
            profile_id: Profile ID
            updates: Dict of fields to update

        Returns:
            Updated profile, or None if not found

        Raises:
            ValueError: If profile is a preset
        """
        profile = await self.get(db, profile_id)
        if profile is None:
            return None
        if profile.is_preset:
            raise ValueError("Cannot modify a preset profile. Duplicate it first.")

        if "name" in updates and updates["name"] is not None:
            profile.name = updates["name"]
        if "description" in updates and updates["description"] is not None:
            profile.description = updates["description"]
        if "thresholds" in updates and updates["thresholds"] is not None:
            profile.thresholds = json.dumps(updates["thresholds"])
        if "weights" in updates and updates["weights"] is not None:
            profile.weights = json.dumps(updates["weights"])

        await db.commit()
        await db.refresh(profile)
        return profile

    async def delete(self, db: AsyncSession, profile_id: int) -> bool:
        """
        Soft-delete a user profile. Refuses to delete presets.

        Returns:
            True if deleted, False if not found

        Raises:
            ValueError: If profile is a preset
        """
        profile = await self.get(db, profile_id)
        if profile is None:
            return False
        if profile.is_preset:
            raise ValueError("Cannot delete a preset profile.")

        profile.is_active = False
        await db.commit()
        return True

    async def duplicate(
        self,
        db: AsyncSession,
        profile_id: int,
        new_name: str,
    ) -> Optional[ScoringProfile]:
        """
        Duplicate any profile (including presets) as a new user profile.

        Returns:
            New profile, or None if source not found
        """
        source = await self.get(db, profile_id)
        if source is None:
            return None

        new_profile = ScoringProfile(
            name=new_name,
            description=f"Copy of {source.name}",
            thresholds=source.thresholds,  # Already JSON string
            weights=source.weights,
            is_preset=False,
            is_active=True,
        )
        db.add(new_profile)
        await db.commit()
        await db.refresh(new_profile)
        return new_profile

    async def export_json(
        self, db: AsyncSession, profile_id: int
    ) -> Optional[dict]:
        """Export a profile as a JSON-serializable dict."""
        profile = await self.get(db, profile_id)
        if profile is None:
            return None
        return {
            "name": profile.name,
            "description": profile.description,
            "thresholds": json.loads(profile.thresholds),
            "weights": json.loads(profile.weights),
        }

    async def import_json(self, db: AsyncSession, data: dict) -> ScoringProfile:
        """Import a profile from a JSON dict."""
        return await self.create(
            db,
            name=data["name"],
            description=data.get("description"),
            thresholds=data.get("thresholds", {}),
            weights=data.get("weights", {}),
        )

    async def seed_presets(self, db: AsyncSession) -> int:
        """
        Seed 8 preset profiles if they do not already exist (idempotent).

        Returns:
            Number of presets created
        """
        # Check how many presets already exist
        result = await db.execute(
            select(ScoringProfile).where(ScoringProfile.is_preset.is_(True))
        )
        existing = {p.name for p in result.scalars().all()}

        created = 0
        for preset in PRESET_PROFILES:
            if preset["name"] not in existing:
                profile = ScoringProfile(
                    name=preset["name"],
                    description=preset["description"],
                    thresholds=json.dumps(preset["thresholds"]),
                    weights=json.dumps(preset["weights"]),
                    is_preset=True,
                    is_active=True,
                )
                db.add(profile)
                created += 1

        if created:
            await db.commit()
            logger.info("Seeded %d preset scoring profiles", created)

        return created


def _profile_to_response_dict(profile: ScoringProfile) -> dict:
    """Convert a ScoringProfile ORM instance to a response-compatible dict."""
    return {
        "id": profile.id,
        "name": profile.name,
        "description": profile.description,
        "thresholds": json.loads(profile.thresholds) if profile.thresholds else {},
        "weights": json.loads(profile.weights) if profile.weights else {},
        "is_preset": profile.is_preset,
        "is_active": profile.is_active,
        "created_at": profile.created_at,
        "updated_at": profile.updated_at,
    }
