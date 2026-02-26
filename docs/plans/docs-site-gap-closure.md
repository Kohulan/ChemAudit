# Documentation Gap Closure Plan

**Branch:** `features_v1`
**Scope:** 7 new pages, 5 page updates, sidebar restructure, API reference overhaul
**Estimated files touched:** ~15 docs-site markdown files + `sidebars.ts`

---

## Context

The `features_v1` branch implements 80+ features across 6 GSD phases (~180 files changed).
The docs-site was largely written before phases 3–6 and covers only the original feature set.
This plan closes all documentation gaps identified in the Feb 26 audit.

---

## Phase 1: New User Guide Pages (5 new pages)

### 1.1 — Bookmarks & History (`docs-site/docs/user-guide/bookmarks-history.md`)

**sidebar_position: 8** (after exporting-results)

Sections to write:

1. **Overview** — what bookmarks and history provide (persistent result snapshots, audit trail)
2. **Bookmarks**
   - Creating a bookmark (BookmarkButton on single validation — star icon, toast)
   - What gets saved (SMILES, InChIKey, name, tags, notes, source, full result snapshot in IndexedDB)
   - Bookmarks page (`/bookmarks`) — search, tag filter chips, expandable rows
   - Managing bookmarks — edit tags/notes, delete individual, bulk delete
   - Opening a bookmark — restores full validation state from IndexedDB snapshot
   - Submit as Batch — select bookmarks → create new batch job from them
   - API reference summary: `POST/GET/PUT/DELETE /bookmarks`, `POST /bookmarks/batch-submit`, `DELETE /bookmarks/bulk`
   - Example curl commands for CRUD
3. **Validation History**
   - What gets logged (every validation event — SMILES, InChIKey, outcome, score, source, timestamp)
   - History page (`/history`) — stats grid (total/passed/warnings/failed)
   - Filtering — date range pickers, outcome dropdown, source filter, SMILES search
   - Paginated results table (25/page, newest first)
   - API reference summary: `GET /history`, `GET /history/stats`
   - Example curl with filter params
4. **Session Scope** — bookmarks and history are scoped to session cookie (anonymous) or API key (authenticated); data isolated between users
5. **Data Retention & Privacy**
   - 30-day auto-purge for anonymous sessions (Celery Beat)
   - GDPR erasure — "Purge My Data" button on Privacy page calls `DELETE /me/data`
   - What gets deleted (all bookmarks + audit entries for session; IndexedDB cleared locally)
   - API-key-scoped data is preserved (not auto-purged)
6. **Next Steps** — links to batch processing, scoring profiles, exporting results

---

### 1.2 — Scoring Profiles (`docs-site/docs/user-guide/scoring/profiles.md`)

**Add to Scoring category in sidebar**, position after aggregator-likelihood

Sections to write:

1. **Overview** — what scoring profiles are (customizable property thresholds + weights for compound filtering)
2. **8 Built-in Presets** — table with name, description, key thresholds for each:
   - Drug-like (Lipinski), Lead-like, Fragment-like (Ro3), CNS-penetrant, Ghose (Amgen), Veber (GSK), PPI-like, NP-like
   - Note: presets are immutable (can duplicate, cannot edit/delete)
3. **Custom Profile Builder** — UI walkthrough:
   - 8 threshold sliders (MW, LogP, TPSA, HBD, HBA, RotBonds, AromRings, Fsp3) with dual min/max range
   - Weight inputs (0–2, step 0.1) per property
   - Save/Cancel/Reset buttons
   - Export/Import JSON for sharing profiles
4. **How Profile Scoring Works**
   - Desirability function: `desirability(value, min, max)` — linear falloff 0–1 outside range
   - Composite score: QED-style weighted geometric mean × 100
   - Per-property desirability breakdown in results
5. **Using Profiles in Batch Processing**
   - ProfileSidebar in upload flow — select profile before upload
   - `profile_id` parameter on `POST /batch/upload`
   - Profile Score column in results table
   - Profile compliance rate in batch statistics (≥80 threshold)
   - Re-scoring subset with different profile via SubsetActionPanel
6. **API Reference** — `GET/POST/PUT/DELETE /profiles`, `POST /profiles/{id}/duplicate`, `GET /profiles/{id}/export`, `POST /profiles/import`
   - Example: create custom profile, apply to batch, query results sorted by profile_score
7. **Profiles Page** (`/profiles`) — PresetPicker tab, Custom Builder tab
8. **Next Steps** — links to batch processing, exporting results

---

### 1.3 — Batch Analytics & Visualizations (`docs-site/docs/user-guide/batch-analytics.md`)

**sidebar_position: 3** (after batch-processing, before structural-alerts)

Sections to write:

1. **Overview** — automatic + on-demand analytics computed after batch processing completes
2. **Automatic Analytics** (run immediately after batch completion)
   - **Deduplication** — 4 levels: exact (canonical SMILES), tautomer, stereo-stripped (InChI), salt-form (ChEMBL parent); total_unique counts per level
   - **Statistics** — per-property stats (mean/median/std/Q1/Q3/IQR/min/max) for validation_score, QED, SA_score, ML_readiness, Fsp3; pairwise Pearson correlations; IQR-fence outlier detection; composite quality score (0–100)
3. **On-Demand Analytics** (user-triggered, for larger datasets)
   - **Scaffold Analysis** — Murcko + generic scaffolds; Shannon entropy diversity; frequency distribution (top 50 + "Other"); R-group decomposition with user-supplied SMARTS core
   - **Chemical Space** — PCA (randomized SVD, pure numpy) or t-SNE (openTSNE, ≤2000 molecules); Morgan ECFP4 2048-bit fingerprints; 2D coordinates + variance explained
   - **MMP & Activity Cliffs** — BRICS single-cut fragmentation; matched molecular pairs (cap 1000); SALI activity cliff index; LLE (requires activity column); batch limit 5000 molecules
   - **Similarity** — find similar molecules (query by SMILES or index, top-k BulkTanimotoSimilarity); nearest neighbors with isolation score; similarity matrix (dense ≤500, sparse 500–2000, refused >2000)
4. **Interactive Visualizations** (BatchAnalyticsPanel)
   - **Distributions Tab**:
     - ScoreHistogram — 4-category colored bars (Excellent/Good/Moderate/Poor); click to select all molecules in band (linked brushing)
     - PropertyScatterPlot — configurable X/Y axes (MW, LogP, TPSA, QED, SA Score, Fsp3); color-by property dropdown; click to toggle selection
     - AlertFrequencyChart — horizontal bars sorted by count; dynamic height
     - ValidationTreemap — hierarchical issue grouping (category → issue → count)
   - **Chemical Space Tab**:
     - ScaffoldTreemap — scaffold frequency treemap; click to select scaffold group
     - ChemicalSpaceScatter — canvas-rendered (handles >800 points); PCA/t-SNE toggle; brush rectangle selection; color-by dropdown; DPR-aware; PNG export
   - **Linked Brushing** — selections sync across all charts via `useBrushSelection` reducer (SET/TOGGLE/ADD_RANGE/CLEAR)
5. **Molecule Comparison Panel**
   - Select 1–2 molecules from results table → "Compare" floating button appears
   - Side-by-side 2D structures (MoleculeViewer 180×140px)
   - Property comparison table with best-value green/red highlighting (Overall Score, QED, SA Score, Fsp3, Lipinski Violations, Alert Count)
   - MoleculePropertyRadar — 6-axis Recharts radar with overlaid profiles + dataset average
   - ECFP4 Tanimoto similarity score between the pair (via `POST /validate/similarity`)
6. **Batch Timeline** — 4-phase horizontal timeline (Upload → Validation → Analytics → Complete) with live status colors and icons
7. **API Reference Summary**
   - `GET /batch/{job_id}/analytics` — get all analytics status + completed results
   - `POST /batch/{job_id}/analytics/{type}` — trigger: `scaffold`, `chemical_space`, `mmp`, `similarity_search`, `rgroup`
   - `POST /validate/similarity` — ECFP4 Tanimoto between two SMILES
   - Example: trigger scaffold analysis, poll for completion, retrieve results
8. **Next Steps** — links to exporting results, scoring profiles, structural alerts

---

### 1.4 — Subset Actions & Sharing (`docs-site/docs/user-guide/subset-actions.md`)

**sidebar_position: 4** (after batch-analytics, before structural-alerts — or merge into batch-analytics)

**Decision: Make this a standalone page for discoverability.**

Sections to write:

1. **Overview** — work with subsets of batch results without re-uploading
2. **Selecting Molecules** — row checkboxes (individual + page select-all); cross-page selection maintained; selection count in floating "Actions" button
3. **SubsetActionPanel** (slide-over)
   - **Re-validate** — creates new batch job with selected molecules; API: `POST /batch/{job_id}/subset/revalidate`
   - **Re-score with Profile** — apply different scoring profile inline (no Celery, synchronous); profile picker dropdown; API: `POST /batch/{job_id}/subset/score-inline`
   - **Export Subset** — export selected molecules only in any format; API: `POST /batch/{job_id}/subset/export`
4. **Batch → Single Navigation**
   - "Open in Single Validation" button per row in BatchResultsTable
   - Navigates to `/` with `fromBatch` state — shows back-to-batch navigation bar with molecule name/index
   - Full single validation view (all 6 tabs) for that molecule
5. **Sharing & Permalinks**
   - Share button creates permalink via `POST /permalinks` — short_id, 30-day expiry
   - URL copied to clipboard automatically
   - Resolving: `GET /report/{short_id}` returns job_id + snapshot (410 Gone if expired)
   - Single-molecule stateless permalink: `/?smiles=<encoded>` (no DB storage, no expiry)
6. **Notifications on Completion**
   - **Email** — `notification_email` param on `POST /batch/upload` or global `NOTIFICATION_EMAIL` config; HTML template with molecule counts and job link; SMTP + STARTTLS
   - **Webhook** — `WEBHOOK_URL` + `WEBHOOK_SECRET` env vars; HMAC-SHA256 signed HTTP POST on batch completion; payload: event, job_id, status, counts, scores, report_url, timestamp; 3 retries with exponential backoff
   - Configuration: set env vars in `.env` or `docker-compose.yml`
7. **Next Steps** — links to exporting results, batch analytics

---

### 1.5 — IUPAC & Input Detection (`docs-site/docs/user-guide/iupac-conversion.md`)

**sidebar_position: 7** (after database-integrations, before exporting-results)

Sections to write:

1. **Overview** — ChemAudit can accept IUPAC names directly as input and converts them to SMILES automatically
2. **How It Works**
   - Auto-detection heuristic in input box: classifies as `smiles` / `iupac` / `ambiguous`
   - Detection rules: SMILES-special chars (`()[]=#@/\+-`) → smiles; spaces or IUPAC suffixes (ane/ene/ine/ol/al/one/acid/yl/ide) → iupac; purely alphabetic ≥4 chars → iupac
   - Input type badge shown below input (e.g., "Detected: IUPAC name")
   - Force-override buttons when detection is ambiguous
3. **Conversion Pipeline**
   - Tries OPSIN first (offline, Java-based, fast) — `NameToStructure.parseChemicalName()`
   - Falls back to PubChem PUG REST API (online, slower) — `/compound/name/{name}/property/CanonicalSMILES/txt`
   - Returns: converted SMILES + conversion source (OPSIN or PubChem)
   - `InputInterpretation` response schema: `detected_type`, `original_input`, `converted_smiles`, `conversion_source`
4. **UI Experience**
   - IUPAC conversion result badge shows original name → converted SMILES
   - Validation proceeds with converted SMILES — all 6 tabs work as normal
   - Placeholder text updated: "Enter SMILES, InChI, or IUPAC name (e.g., aspirin)"
5. **API Usage**
   - `POST /validate` with `{"molecule": "aspirin"}` — auto-detects and converts
   - Response includes `input_interpretation` field with detection/conversion details
   - Example curl + response JSON
6. **Supported Names** — works with systematic names (2-acetoxybenzoic acid), common names (aspirin), trade names (if in PubChem)
7. **Limitations** — OPSIN requires Java JRE (provisioned in Docker); PubChem fallback requires internet; very obscure names may fail; ambiguous names may need force-override
8. **Next Steps** — links to single validation, database integrations

---

## Phase 2: Existing Page Updates (5 pages)

### 2.1 — Update `single-validation.md`

Add/modify these sections:

1. **Input Formats table** — add row for IUPAC names with link to IUPAC conversion page
2. **Tabs list** — update from 5 tabs to 6 tabs:
   - Add "Scoring Profiles" tab with brief description (consensus score, lead/fragment-likeness, property breakdowns, bioavailability radar, BOILED-Egg) and link to scoring/profiles page
   - Add "Deep Validation" tab description with link (already partially documented)
3. **Scoring Profiles Tab section** (new section after existing scoring info):
   - ConsensusScoreCard — consensus drug-likeness 0–5 across 5 rule sets
   - LeadFragmentCard — lead-likeness, Rule of 3, salt inventory, ligand efficiency
   - PropertyBreakdownCard — TPSA/LogP per-atom breakdowns, Bertz complexity, Fsp3 detail
   - BioavailabilityCard — 6-axis radar chart + BOILED-Egg scatter
   - AtomContributionViewer — per-atom heatmap
   - Link to scoring/profiles for full documentation
4. **Standardization Provenance section** (new section in standardization tab info):
   - ProvenanceTimeline — vertical stage timeline
   - Per-stage changes tracked: charge, bond, ring, radical, fragment removal
   - DVAL cross-references (e.g., "DVAL-01: N undefined stereocenters")
   - `include_provenance=true` option
5. **Bookmarking** — mention BookmarkButton (star icon) and link to bookmarks-history page
6. **RDKit Version tooltip** — mention the RDKit version source tooltip on canonical SMILES in molecule info
7. **Severity Configuration Panel** — expand the brief mention into a proper sub-section:
   - Gear icon → slide-in panel
   - Per-check severity override (ERROR/WARNING/INFO or reset to default)
   - Overrides persist to localStorage
   - Dynamic verdict recomputed from effective severities
8. **IUPAC input** — mention auto-detection and link to iupac-conversion page

---

### 2.2 — Update `batch-processing.md`

Add/modify these sections:

1. **Processing Options table** — add rows for:
   - `profile_id` — apply scoring profile to batch
   - `notification_email` — email address for completion notification
   - `include_analytics` — auto-compute analytics (default: true)
2. **Upload flow** — mention ProfileSidebar (select scoring profile before upload)
3. **Sort Fields table** — add `profile_score` sort field
4. **Results section** — mention:
   - Profile Score column (conditional, shown when profile was applied)
   - Profile compliance rate in statistics
   - "Open in Single Validation" per-row button
   - Molecule comparison (select 2 → Compare button)
5. **Analytics section** (new, brief — link to batch-analytics page):
   - After batch completes, automatic analytics run (deduplication, statistics)
   - On-demand analytics available (scaffold, chemical space, MMP)
   - Interactive visualizations below results table
   - Link to full batch-analytics page
6. **Subset Actions section** (new, brief — link to subset-actions page):
   - Select molecules → Actions button → re-validate, re-score, export subset
   - Link to full subset-actions page
7. **Sharing section** (new, brief):
   - Share button creates shareable permalink
   - Webhook/email notifications on completion
   - Link to subset-actions page
8. **Next Steps** — add links to batch-analytics, subset-actions, scoring profiles

---

### 2.3 — Update `exporting-results.md`

Add/modify these sections:

1. **Available Export Formats table** — add 4 new rows:
   - Fingerprint Matrix (`.zip`) — Morgan/MACCS/RDKit FPs in CSV + numpy formats
   - Deduplicated Set (`.zip`) — summary + annotated CSVs with group IDs
   - Scaffold-Grouped (`.csv`) — molecules with Murcko scaffold + group assignment
   - Property Matrix (`.zip`) — flat CSV + 4-sheet Excel (Descriptors, Scores, Alerts, Properties)
2. **Format Details** — add 4 new sub-sections:
   - **Fingerprint Matrix Export** — 9 files in ZIP: 3 FP types (Morgan ECFP4 2048-bit, MACCS 167-bit, RDKit path 2048-bit) × 3 formats (CSV, NPY, NPZ); skips error molecules; ideal for ML pipeline input
   - **Deduplicated Set Export** — 2 files in ZIP: `dedup_summary.csv` (one row per unique group: representative, canonical SMILES, group size, member indices, dedup level) + `dedup_annotated.csv` (all molecules with group_id, is_representative, dedup_level, score, InChIKey); ideal for dataset curation
   - **Scaffold-Grouped Export** — single CSV with `scaffold_smiles` and integer `scaffold_group` columns; acyclic molecules get empty scaffold and group 0; ideal for SAR analysis
   - **Property Matrix Export** — 2 files in ZIP: `properties_flat.csv` + `properties.xlsx` with 4 sheets (Descriptors: MW/LogP/TPSA/HBD/HBA/RotBonds/AromRings/Fsp3/HAC/RingCount; Scores: overall/ML/QED/SA/NP/Lipinski; Alerts: per-catalog counts; Properties: SMILES/name/InChIKey/canonical); ideal for comprehensive analysis
3. **PDF Section Selection** — document the `sections` parameter:
   - 8 configurable sections: `validation_summary`, `score_distribution`, `alert_frequency`, `chemical_space`, `scaffold_treemap`, `statistics`, `correlation_matrix`, `mmp_pairs`
   - All sections included by default
   - Example: `?format=pdf&sections=validation_summary,score_distribution`
4. **Excel Images** — document `include_images` parameter for embedded 2D structure PNGs
5. **ExportDialog UI** — mention the 9-format dialog with icons, descriptions, "New" badges, file name preview, format-specific extensions

---

### 2.4 — Update `standardization.md`

Add a new section:

1. **Standardization Provenance** (after existing "Understanding Results"):
   - Enable with `include_provenance=true` on `/standardize`
   - Returns per-stage provenance records with typed change arrays
   - Stage types: Tautomer, ChargeNormalization, BondNormalization, FragmentRemoval, Stereo
   - Change types tracked:
     - ChargeChange — atom_idx, element, before/after charge, rule_name, SMARTS
     - BondChange — bond_idx, atom1, atom2, before/after type, rule_name
     - RingChange — ring_atoms, ring_size, before/after aromaticity type
     - RadicalChange — atom_idx, element, before/after radical electrons
     - Fragment removal — SMILES, name (from 55-entry dictionary), role, MW
   - Tautomer provenance: canonical_smiles, num_tautomers_enumerated, complexity flag (>100 tautomers)
   - Stereo detail: per-center type/before_config/after_config/reason
   - DVAL cross-references: links to deep validation results (e.g., "DVAL-03: N tautomers enumerated")
   - API example with `include_provenance=true` and sample response JSON
   - UI: ProvenanceTimeline (vertical, auto-expands changed stages) with ProvenanceStageCards

---

### 2.5 — Update `api/endpoints.md`

Add these missing endpoint sections:

1. **Bookmarks** section:
   - `GET /bookmarks` — paginated list with filters (tag, search, source)
   - `POST /bookmarks` — create (SMILES, name, tags, notes)
   - `GET /bookmarks/{id}` — single bookmark
   - `PUT /bookmarks/{id}` — update tags/notes
   - `DELETE /bookmarks/{id}` — delete single
   - `POST /bookmarks/batch-submit` — submit as batch job
   - `DELETE /bookmarks/bulk` — bulk delete by IDs
2. **History** section:
   - `GET /history` — paginated audit trail with filters (date_from, date_to, outcome, source, smiles_search)
   - `GET /history/stats` — total, outcome distribution, source distribution
3. **Scoring Profiles** section:
   - `GET /profiles` — list all (8 presets + user)
   - `POST /profiles` — create custom
   - `GET /profiles/{id}` — single
   - `PUT /profiles/{id}` — update (400 if preset)
   - `DELETE /profiles/{id}` — soft-delete (400 if preset)
   - `POST /profiles/{id}/duplicate` — duplicate any (including presets)
   - `GET /profiles/{id}/export` — export JSON
   - `POST /profiles/import` — import JSON
4. **Batch Analytics** section:
   - `GET /batch/{job_id}/analytics` — all analytics status + results
   - `POST /batch/{job_id}/analytics/{type}` — trigger expensive analytics
5. **Batch Subset Actions** section:
   - `POST /batch/{job_id}/subset/revalidate` — re-validate subset
   - `POST /batch/{job_id}/subset/rescore` — re-score with profile
   - `POST /batch/{job_id}/subset/export` — export subset
   - `POST /batch/{job_id}/subset/score-inline` — synchronous inline scoring
6. **Similarity** section:
   - `POST /validate/similarity` — ECFP4 Tanimoto between two SMILES
7. **Permalinks** section:
   - `POST /permalinks` — create batch report permalink
   - `GET /report/{short_id}` — resolve permalink (410 if expired)
8. **Session** section:
   - `DELETE /me/data` — GDPR erasure (all bookmarks + history for session)

Each endpoint includes: method, path, description, request body/params, response JSON example, rate limit note.

---

## Phase 3: Sidebar Restructure & Introduction Update

### 3.1 — Update `sidebars.ts`

New sidebar structure:

```typescript
const sidebars: SidebarsConfig = {
  docs: [
    { type: 'doc', id: 'intro', label: 'Introduction' },
    {
      type: 'category',
      label: 'Getting Started',
      collapsed: false,
      items: [
        'getting-started/installation',
        'getting-started/configuration',
        'getting-started/first-validation',
      ],
    },
    {
      type: 'category',
      label: 'User Guide',
      collapsed: false,
      items: [
        'user-guide/single-validation',
        'user-guide/batch-processing',
        'user-guide/batch-analytics',        // NEW
        'user-guide/subset-actions',          // NEW
        'user-guide/structural-alerts',
        {
          type: 'category',
          label: 'Scoring',
          items: [
            'user-guide/scoring/overview',
            'user-guide/scoring/ml-readiness',
            'user-guide/scoring/drug-likeness',
            'user-guide/scoring/safety-filters',
            'user-guide/scoring/admet',
            'user-guide/scoring/np-likeness',
            'user-guide/scoring/scaffold-analysis',
            'user-guide/scoring/aggregator-likelihood',
            'user-guide/scoring/profiles',    // NEW
          ],
        },
        'user-guide/standardization',
        'user-guide/database-integrations',
        'user-guide/iupac-conversion',        // NEW
        'user-guide/exporting-results',
        'user-guide/bookmarks-history',       // NEW
      ],
    },
    {
      type: 'category',
      label: 'API Reference',
      collapsed: true,
      items: [
        'api/overview',
        'api/authentication',
        'api/endpoints',
        'api/websocket',
        'api/error-handling',
        'api/rate-limits',
      ],
    },
    {
      type: 'category',
      label: 'Deployment',
      collapsed: true,
      items: [
        'deployment/docker',
        'deployment/production',
        'deployment/monitoring',
      ],
    },
    { type: 'doc', id: 'troubleshooting', label: 'Troubleshooting' },
  ],
};
```

### 3.2 — Update `intro.md`

1. **What ChemAudit Does** bullet list — add:
   - **Bookmark** results and track validation **history**
   - **Compare** molecules side-by-side with property radar overlays
   - **Analyze** batch datasets with interactive charts and chemical space maps
   - **Share** results via permalinks and receive **notifications** on completion
2. **Key Features** — add sections:
   - **Batch Analytics & Visualizations** — interactive charts, chemical space mapping, deduplication, scaffold analysis, MMP detection
   - **Bookmarks & History** — save result snapshots, audit trail, tag filtering, submit bookmarks as batch
   - **Scoring Profiles** — 8 built-in presets, custom profile builder, batch profile scoring
   - **IUPAC Name Support** — auto-detect IUPAC names, convert via OPSIN/PubChem
   - **Notifications** — email and webhook on batch completion
3. **Navigation** — add links to new pages

---

## Phase 4: Cross-Reference & Quality Pass

### 4.1 — Cross-link audit

Ensure all new pages have:
- Proper `sidebar_position` frontmatter matching the new sidebar order
- "Next Steps" sections with relevant cross-links
- Consistent formatting (tables, code blocks, admonitions)

### 4.2 — Update existing "Next Steps" sections

Pages that should link to new content:
- `batch-processing.md` → batch-analytics, subset-actions, scoring/profiles
- `exporting-results.md` → batch-analytics, subset-actions
- `scoring/overview.md` → scoring/profiles
- `single-validation.md` → bookmarks-history, iupac-conversion, scoring/profiles
- `standardization.md` → (provenance is internal, just ensure self-contained)
- `database-integrations.md` → iupac-conversion
- `first-validation.md` → update tabs list to 6, mention bookmarks
- `configuration.md` → add SMTP/webhook env vars to config table
- `troubleshooting.md` → add troubleshooting for IUPAC (Java/OPSIN), bookmarks, webhooks

### 4.3 — Update `configuration.md`

Add new environment variables to the config table:
- `NOTIFICATION_EMAIL` — global email for batch completion notifications
- `WEBHOOK_URL` — webhook endpoint URL
- `WEBHOOK_SECRET` — HMAC signing secret for webhooks
- `OPSIN_JAR_PATH` — path to OPSIN JAR (auto-provisioned in Docker)
- `SMTP_HOST`, `SMTP_PORT`, `SMTP_USER`, `SMTP_PASS` — email server settings
- `BASE_URL` — public URL for report links in emails/webhooks

---

## Execution Order

```
Phase 1 (New pages — do first, most content):
  1.1  bookmarks-history.md          ~250 lines
  1.2  scoring/profiles.md           ~300 lines
  1.3  batch-analytics.md            ~400 lines  (largest page)
  1.4  subset-actions.md             ~250 lines
  1.5  iupac-conversion.md           ~180 lines

Phase 2 (Existing page updates — after new pages exist for cross-links):
  2.1  single-validation.md update   ~120 lines added
  2.2  batch-processing.md update    ~100 lines added
  2.3  exporting-results.md update   ~200 lines added
  2.4  standardization.md update     ~100 lines added
  2.5  api/endpoints.md update       ~350 lines added

Phase 3 (Structural):
  3.1  sidebars.ts update            ~5 lines changed
  3.2  intro.md update               ~60 lines added

Phase 4 (Polish):
  4.1  Cross-link audit              all pages
  4.2  Next Steps updates            ~8 pages
  4.3  configuration.md update       ~20 lines added
```

**Total estimated new content: ~2,300 lines across 15 files.**

---

## Quality Checklist

For each new/updated page, verify:

- [ ] Frontmatter correct (sidebar_position, title, description)
- [ ] Feature matches actual implementation (cross-check with codebase)
- [ ] API examples use correct endpoints, params, and response shapes
- [ ] All code blocks have language annotations (bash, json, typescript, python)
- [ ] Admonitions used appropriately (:::tip, :::info, :::warning, :::caution)
- [ ] Tables render correctly in Docusaurus
- [ ] No broken internal links
- [ ] "Next Steps" section present with relevant links
- [ ] No duplicate content (link to primary page instead of repeating)
- [ ] Consistent tone and formatting with existing docs

---

## What NOT to Document

- Internal implementation details (React hook internals, Celery task signatures)
- Visual polish (count-up animations, Framer Motion transitions, ScoreChart variants)
- IndexedDB storage internals (mention as "locally cached", don't explain IDB API)
- Fragment dictionary contents (55 entries — too detailed)
- SMARTS patterns in aggregator/alerts (already in scoring pages)
- Provenance schema field names (show example JSON instead)
