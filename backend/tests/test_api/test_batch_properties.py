"""Test that batch results include rotatable_bonds and aromatic_rings."""
from app.services.batch.tasks import _process_single_molecule


def test_druglikeness_includes_rotatable_bonds_and_aromatic_rings():
    """rotatable_bonds and aromatic_rings should appear in druglikeness dict."""
    mol_data = {"smiles": "c1ccc(CC(=O)O)cc1", "name": "test", "index": 0}
    result = _process_single_molecule(mol_data)
    dl = result["scoring"]["druglikeness"]
    assert "rotatable_bonds" in dl
    assert "aromatic_rings" in dl
    assert isinstance(dl["rotatable_bonds"], int)
    assert isinstance(dl["aromatic_rings"], int)
    # Phenylacetic acid: 1 aromatic ring, 2 rotatable bonds
    assert dl["aromatic_rings"] == 1
    assert dl["rotatable_bonds"] >= 1
