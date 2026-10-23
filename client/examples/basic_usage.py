"""
Basic usage examples for ChemAudit client.

Demonstrates single molecule validation, alert screening, standardization, and scoring.
"""

from chemaudit import ChemAuditClient

# Configuration
API_BASE_URL = "http://localhost:8000"
API_KEY = None  # Set to your API key for higher rate limits


def main():
    print("ChemAudit Client - Basic Usage Examples\n")
    print("=" * 60)

    # Create client (use context manager for automatic cleanup)
    with ChemAuditClient(base_url=API_BASE_URL, api_key=API_KEY) as client:
        # Example 1: Simple validation
        print("\n1. Simple Validation")
        print("-" * 60)
        result = client.validate("CCO")
        print(f"Molecule: {result.molecule_info.canonical_smiles}")
        print(f"Overall Score: {result.overall_score}/100")
        print(f"Number of Issues: {len(result.issues)}")
        if result.issues:
            print("\nIssues Found:")
            for issue in result.issues:
                print(f"  - {issue.check_name}: {issue.message}")

        # Example 2: Comprehensive validation
        print("\n2. Validation with an explicit input format")
        print("-" * 60)
        result = client.validate(
            molecule="c1ccccc1C(=O)O",  # Benzoic acid
            format="smiles",
        )
        print(f"Molecule: {result.molecule_info.canonical_smiles}")
        print(f"Validation Score: {result.overall_score}/100")
        print(f"Molecular Weight: {result.molecule_info.molecular_weight:.2f}")
        print(f"InChIKey: {result.molecule_info.inchikey}")

        # Example 3: Structural alert screening
        print("\n3. Structural Alert Screening")
        print("-" * 60)
        alerts = client.screen_alerts(
            molecule="c1ccccc1",
            catalogs=["PAINS", "BRENK"],
        )
        print(f"Catalogs Screened: {', '.join(alerts.screened_catalogs)}")
        print(f"Total Alerts: {alerts.total_alerts}")
        if alerts.alerts:
            print("\nAlerts Found:")
            for alert in alerts.alerts[:3]:  # Show first 3
                print(f"  - {alert.catalog_source}: {alert.description}")
                print(f"    Severity: {alert.severity}")
        if alerts.educational_note:
            print(f"\nNote: {alerts.educational_note}")

        # Example 4: ML-readiness scoring
        print("\n4. ML-Readiness Scoring")
        print("-" * 60)
        scores = client.score("CC(C)Cc1ccc(cc1)C(C)C(=O)O")  # Ibuprofen
        print("Scores:")
        if scores.ml_readiness:
            ml = scores.ml_readiness
            print(f"  - ML-Readiness: {ml['score']} ({ml['label']})")
        if scores.np_likeness:
            print(f"  - NP-Likeness: {scores.np_likeness['score']:.2f}")
        if scores.scaffold:
            print(f"  - Murcko Scaffold: {scores.scaffold['scaffold_smiles']}")

        # Example 5: Molecule standardization
        print("\n5. Molecule Standardization")
        print("-" * 60)
        result = client.standardize(
            molecule="CCN.Cl",  # Ethylamine hydrochloride
            tautomer=False,  # Preserve stereochemistry
        )
        print(f"Input: {result.result.original_smiles}")
        print(f"Standardized: {result.result.standardized_smiles}")
        print("\nSteps:")
        for step in result.result.steps_applied:
            print(f"  - {step['step_name']}: {step.get('changes', '')}")

        # Example 6: Invalid molecule handling
        print("\n6. Error Handling")
        print("-" * 60)
        try:
            result = client.validate("INVALID_SMILES_123")
        except Exception as e:
            print(f"Error Type: {type(e).__name__}")
            print(f"Error Message: {e}")

    print("\n" + "=" * 60)
    print("Examples completed successfully!")


if __name__ == "__main__":
    main()
