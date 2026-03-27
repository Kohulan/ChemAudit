"""
InChI Layer Diff Service (DIAG-02).

Provides pure string parsing of InChI strings and layer-by-layer comparison.
No RDKit dependency — InChI format is slash-delimited with single-character layer keys.
"""

# Mapping from InChI layer key character to human-readable layer name
LAYER_KEY_NAMES = {
    "c": "connections",
    "h": "hydrogens",
    "q": "charge",
    "b": "double_bond_stereo",
    "t": "tetrahedral_stereo",
    "m": "stereo_parity",
    "s": "stereo_type",
    "i": "isotope",
    "p": "protons",
}

# Display order for layer rows (unknown layers sort to end)
LAYER_DISPLAY_ORDER = [
    "formula",
    "connections",
    "hydrogens",
    "charge",
    "double_bond_stereo",
    "tetrahedral_stereo",
    "stereo_parity",
    "stereo_type",
    "isotope",
    "protons",
]

_LAYER_ORDER_INDEX = {name: i for i, name in enumerate(LAYER_DISPLAY_ORDER)}


def parse_inchi_layers(inchi: str) -> dict:
    """Parse an InChI string into its constituent layers.

    Accepts any string starting with "InChI=" (including non-standard versions).
    The formula layer is the first segment after the version prefix that starts
    with an uppercase letter.

    Args:
        inchi: InChI string, e.g. "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"

    Returns:
        Dict mapping layer name to layer value string (without the leading key char).
        Special key "version" holds the InChI version string (e.g. "1S").
        Special key "formula" holds the molecular formula.
    """
    parts = inchi.split("/")
    result: dict[str, str] = {}

    if not parts:
        return result

    # First part is version, e.g. "InChI=1S"
    version_part = parts[0]
    if "=" in version_part:
        result["version"] = version_part.split("=", 1)[1]
    else:
        result["version"] = version_part

    for segment in parts[1:]:
        if not segment:
            continue
        first_char = segment[0]
        if first_char.isupper():
            # Molecular formula layer
            result["formula"] = segment
        elif first_char in LAYER_KEY_NAMES:
            layer_name = LAYER_KEY_NAMES[first_char]
            result[layer_name] = segment[1:]  # Strip the leading key character
        else:
            # Unknown layer — store under its raw key char
            result[first_char] = segment[1:]

    return result


def diff_inchi_layers(inchi_a: str, inchi_b: str) -> dict:
    """Compare two InChI strings layer by layer.

    Args:
        inchi_a: First InChI string.
        inchi_b: Second InChI string.

    Returns:
        Dict with keys:
            identical (bool): True if all layers match.
            layer_rows (list[dict]): One row per layer with keys:
                layer (str), value_a (str | None), value_b (str | None), match (bool).
            layers_a (dict): Parsed layers from inchi_a.
            layers_b (dict): Parsed layers from inchi_b.
    """
    layers_a = parse_inchi_layers(inchi_a)
    layers_b = parse_inchi_layers(inchi_b)

    # Build union of layer keys, excluding "version"
    all_keys = set(layers_a.keys()) | set(layers_b.keys())
    all_keys.discard("version")

    def _sort_key(layer_name: str) -> int:
        return _LAYER_ORDER_INDEX.get(layer_name, 99)

    sorted_keys = sorted(all_keys, key=_sort_key)

    layer_rows = []
    all_match = True
    for key in sorted_keys:
        val_a = layers_a.get(key)
        val_b = layers_b.get(key)
        match = val_a == val_b
        if not match:
            all_match = False
        layer_rows.append({
            "layer": key,
            "value_a": val_a,
            "value_b": val_b,
            "match": match,
        })

    return {
        "identical": all_match,
        "layer_rows": layer_rows,
        "layers_a": {k: v for k, v in layers_a.items() if k != "version"},
        "layers_b": {k: v for k, v in layers_b.items() if k != "version"},
    }
