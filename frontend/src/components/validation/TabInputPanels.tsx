import { type ComponentProps } from 'react';
import { AlertTriangle, ChevronDown, FlaskConical, Info, Layers, Play, Sparkles } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { InfoTooltip, DoiLink } from '../ui/Tooltip';
import { cn } from '../../lib/utils';
import { DeepValidationTab } from './DeepValidationTab';
import { ScoringProfilesTab } from '../scoring-profiles';
import { ProfilerAccordion } from '../profiler/ProfilerAccordion';
import { SafetyAccordion } from '../safety/SafetyAccordion';

/**
 * Left-column input panels for each analysis tab of the Single Validation
 * page. All are presentational — the page owns the molecule input, results,
 * loading flags, and handlers, and passes them in as props. Each panel is
 * rendered by the page under its own `activeTab === ...` guard.
 */

interface ValidateButtonProps {
  onValidate: () => void;
  disabled: boolean;
  loading: boolean;
}

/** Primary "Validate" trigger shared by every tab that gates its content on a
 * successful basic validation. */
function ValidateButton({ onValidate, disabled, loading }: ValidateButtonProps) {
  return (
    <ClayButton
      variant="primary"
      onClick={onValidate}
      disabled={disabled}
      loading={loading}
      leftIcon={<Play className="w-4 h-4" />}
    >
      Validate
    </ClayButton>
  );
}

interface ValidateTabPanelProps {
  onValidate: () => void;
  onScore: () => void;
  /** !molecule.trim() || isAnyLoading */
  actionsDisabled: boolean;
  hasResult: boolean;
  validateLoading: boolean;
  scoreLoading: boolean;
  /** !result && molecule has content — show the "score after validation" hint. */
  showScoreHint: boolean;
}

export function ValidateTabPanel({
  onValidate,
  onScore,
  actionsDisabled,
  hasResult,
  validateLoading,
  scoreLoading,
  showScoreHint,
}: ValidateTabPanelProps) {
  return (
    <div className="space-y-4">
      <div className="flex items-start gap-3">
        <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
          <Info className="w-4 h-4" />
        </div>
        <p className="text-[var(--color-text-secondary)] text-sm">
          Run Validate first. Score becomes available once the canonical SMILES is computed.
        </p>
      </div>
      <div className="flex flex-wrap gap-3">
        <ValidateButton onValidate={onValidate} disabled={actionsDisabled} loading={validateLoading} />
        <ClayButton
          variant="accent"
          onClick={onScore}
          disabled={actionsDisabled || !hasResult}
          loading={scoreLoading}
          leftIcon={<Sparkles className="w-4 h-4" />}
          title={!hasResult ? 'Run Validate first — Score uses the canonical SMILES from validation' : undefined}
        >
          Score
        </ClayButton>
      </div>
      {showScoreHint && (
        <p className="text-xs text-[var(--color-text-muted)]">
          Score becomes available after validation completes.
        </p>
      )}
    </div>
  );
}

interface DeepValidationTabPanelProps {
  onValidate: () => void;
  actionsDisabled: boolean;
  validateLoading: boolean;
  /** result?.all_checks ?? null — null shows the Validate button instead of the deep checks. */
  checks: ComponentProps<typeof DeepValidationTab>['checks'] | null;
  onHighlightAtoms: ComponentProps<typeof DeepValidationTab>['onHighlightAtoms'];
}

export function DeepValidationTabPanel({
  onValidate,
  actionsDisabled,
  validateLoading,
  checks,
  onHighlightAtoms,
}: DeepValidationTabPanelProps) {
  return (
    <div className="space-y-4">
      <div className="flex items-start gap-3">
        <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
          <Info className="w-4 h-4" />
        </div>
        <p className="text-[var(--color-text-secondary)] text-sm">
          Stereoisomer enumeration, tautomer detection, composition guards, and complexity flags.
          Requires a successful basic validation first.
        </p>
      </div>
      {!checks && (
        <ValidateButton onValidate={onValidate} disabled={actionsDisabled} loading={validateLoading} />
      )}
      {checks && (
        <DeepValidationTab
          checks={checks}
          onHighlightAtoms={onHighlightAtoms}
        />
      )}
    </div>
  );
}

interface ScoringProfilesTabPanelProps {
  onValidate: () => void;
  actionsDisabled: boolean;
  validateLoading: boolean;
  hasResult: boolean;
  /** result ? resolvedSmiles : '' */
  smiles: string;
}

export function ScoringProfilesTabPanel({
  onValidate,
  actionsDisabled,
  validateLoading,
  hasResult,
  smiles,
}: ScoringProfilesTabPanelProps) {
  return (
    <div className="space-y-4">
      {!hasResult && (
        <>
          <div className="flex items-start gap-3">
            <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
              <Info className="w-4 h-4" />
            </div>
            <p className="text-[var(--color-text-secondary)] text-sm">
              Drug-likeness, lead-likeness, property breakdowns, and bioavailability radar
              across multiple consensus profiles. Run Validate first to populate.
            </p>
          </div>
          <ValidateButton onValidate={onValidate} disabled={actionsDisabled} loading={validateLoading} />
        </>
      )}
      <ScoringProfilesTab smiles={smiles} />
    </div>
  );
}

interface CompoundProfileTabPanelProps {
  onValidate: () => void;
  actionsDisabled: boolean;
  validateLoading: boolean;
  hasResult: boolean;
  /** canonicalSmiles || '' */
  smiles: string;
  profile: ComponentProps<typeof ProfilerAccordion>['profile'];
  profileLoading: boolean;
  profileError: ComponentProps<typeof ProfilerAccordion>['error'];
}

export function CompoundProfileTabPanel({
  onValidate,
  actionsDisabled,
  validateLoading,
  hasResult,
  smiles,
  profile,
  profileLoading,
  profileError,
}: CompoundProfileTabPanelProps) {
  return (
    <div className="space-y-4">
      <div className="flex items-start gap-3">
        <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
          <FlaskConical className="w-4 h-4" />
        </div>
        <div className="flex-1 min-w-0">
          <div className="flex items-start gap-1.5">
            <p className="text-[var(--color-text-secondary)] text-sm flex-1">
              Comprehensive molecular profiling with physicochemical properties, drug-likeness assessment,
              synthesizability comparison, and 3D shape analysis.
            </p>
            <InfoTooltip
              title="Profiling Metrics"
              position="bottom"
              content={
                <div className="text-xs space-y-2">
                  <div>
                    <p className="font-semibold text-white">PFI (Property Forecast Index)</p>
                    <p className="text-white/70">cLogP + aromatic ring count. &lt;5 low, 5-7 moderate, &gt;7 high risk.</p>
                    <p className="text-white/50 italic">Young et al. Drug Discov Today (2011)</p>
                    <DoiLink doi="10.1016/j.drudis.2011.06.001" />
                  </div>
                  <div>
                    <p className="font-semibold text-white">#Stars (Outlier Count)</p>
                    <p className="text-white/70">Properties outside 95th-percentile drug-like ranges. 0 = drug-like, 3+ = outlier.</p>
                    <p className="text-white/50 italic">Jorgensen &amp; Duffy. Adv Drug Deliv Rev (2002)</p>
                    <DoiLink doi="10.1016/S0169-409X(02)00008-X" />
                  </div>
                  <div>
                    <p className="font-semibold text-white">Abbott Bioavailability Score</p>
                    <p className="text-white/70">4-class oral bioavailability probability (11%, 17%, 56%, 85%).</p>
                    <p className="text-white/50 italic">Martin. J Med Chem (2005)</p>
                    <DoiLink doi="10.1021/jm0492002" />
                  </div>
                  <div>
                    <p className="font-semibold text-white">Drug-likeness Rules</p>
                    <p className="text-white/70">Lipinski, Veber, Egan, Muegge, Ghose filter evaluation.</p>
                    <p className="text-white/50 italic">Lipinski et al. Adv Drug Deliv Rev (2001)</p>
                    <DoiLink doi="10.1016/S0169-409X(00)00129-0" />
                  </div>
                  <div>
                    <p className="font-semibold text-white">SA Comparison</p>
                    <p className="text-white/70">Synthetic accessibility: SA Score + SCScore + SYBA side-by-side.</p>
                    <p className="text-white/50 italic">Ertl &amp; Schuffenhauer. J Cheminform (2009)</p>
                    <DoiLink doi="10.1186/1758-2946-1-8" />
                  </div>
                </div>
              }
            />
          </div>
        </div>
      </div>
      {!hasResult && (
        <ValidateButton onValidate={onValidate} disabled={actionsDisabled} loading={validateLoading} />
      )}
      <ProfilerAccordion
        smiles={smiles}
        profile={profile}
        isLoading={profileLoading}
        error={profileError}
      />
    </div>
  );
}

interface AlertsTabPanelProps {
  /** Catalogs available in the collapsible ChEMBL group. */
  chemblCatalogs: ReadonlyArray<{ id: string; label: string }>;
  selectedCatalogs: string[];
  onToggleCatalog: (catalog: string) => void;
  chemblExpanded: boolean;
  onToggleChemblExpanded: () => void;
  allChemblSelected: boolean;
  someChemblSelected: boolean;
  onToggleAllChembl: (enabled: boolean) => void;
  onScreenAlerts: () => void;
  /** !molecule.trim() || isAnyLoading || selectedCatalogs.length === 0 */
  screenDisabled: boolean;
  alertsLoading: boolean;
  /** canonicalSmiles || '' */
  safetySmiles: string;
  safetyAlertResult: ComponentProps<typeof SafetyAccordion>['alertResult'];
  safetyResult: ComponentProps<typeof SafetyAccordion>['safetyResult'];
  safetyLoading: boolean;
  safetyError: ComponentProps<typeof SafetyAccordion>['error'];
}

export function AlertsTabPanel({
  chemblCatalogs,
  selectedCatalogs,
  onToggleCatalog,
  chemblExpanded,
  onToggleChemblExpanded,
  allChemblSelected,
  someChemblSelected,
  onToggleAllChembl,
  onScreenAlerts,
  screenDisabled,
  alertsLoading,
  safetySmiles,
  safetyAlertResult,
  safetyResult,
  safetyLoading,
  safetyError,
}: AlertsTabPanelProps) {
  return (
    <div className="space-y-4">
      <div className="flex items-start gap-3">
        <div className="w-8 h-8 rounded-lg bg-amber-500/10 flex items-center justify-center text-amber-500 flex-shrink-0">
          <Info className="w-4 h-4" />
        </div>
        <div className="flex-1 min-w-0">
          <div className="flex items-start gap-1.5">
            <p className="text-[var(--color-text-secondary)] text-sm flex-1">
              PAINS and BRENK preselected. Toggle additional catalogs (NIH, ZINC, ChEMBL filters) below.
            </p>
            <InfoTooltip
              title="Safety Assessment Metrics"
              position="bottom"
              content={
                <div className="text-xs space-y-2">
                  <div>
                    <p className="font-semibold text-white">CYP Soft-Spots</p>
                    <p className="text-white/70">SMARTS-based cytochrome P450 metabolism site prediction with atom highlighting.</p>
                    <p className="text-white/50 italic">Rydberg et al. ACS Med Chem Lett (2010)</p>
                    <DoiLink doi="10.1021/ml100016x" />
                  </div>
                  <div>
                    <p className="font-semibold text-white">hERG Risk</p>
                    <p className="text-white/70">Rule-based hERG channel liability assessment using amphiphilic properties.</p>
                    <p className="text-white/50 italic">Aronov. Drug Discov Today (2005)</p>
                    <DoiLink doi="10.1016/S1359-6446(04)03278-7" />
                  </div>
                  <div>
                    <p className="font-semibold text-white">bRo5 (Beyond Rule of 5)</p>
                    <p className="text-white/70">Relaxed thresholds for macrocycles, PROTACs, and natural products (MW &gt; 500).</p>
                    <p className="text-white/50 italic">Doak et al. Chem Biol (2014)</p>
                    <DoiLink doi="10.1016/j.chembiol.2014.08.013" />
                  </div>
                  <div>
                    <p className="font-semibold text-white">REOS Filter</p>
                    <p className="text-white/70">Rapid Elimination of Swill — 7 physicochemical property range filters.</p>
                    <p className="text-white/50 italic">Walters &amp; Murcko. Curr Opin Chem Biol (1999)</p>
                    <DoiLink doi="10.1016/S1367-5931(99)80058-1" />
                  </div>
                  <div>
                    <p className="font-semibold text-white">Complexity Analysis</p>
                    <p className="text-white/70">Bertz complexity index percentile vs commercial compound distributions.</p>
                    <p className="text-white/50 italic">Bertz. J Am Chem Soc (1981)</p>
                    <DoiLink doi="10.1021/ja00402a071" />
                  </div>
                </div>
              }
            />
          </div>
        </div>
      </div>

      {/* Catalog selector */}
      <div>
        <p className="text-xs text-[var(--color-text-muted)] mb-2">Select catalogs to screen:</p>
        {/* Core catalogs */}
        <div className="flex flex-wrap gap-2">
          {['PAINS', 'BRENK', 'NIH', 'ZINC'].map((catalog) => (
            <button
              key={catalog}
              onClick={() => onToggleCatalog(catalog)}
              className={cn(
                'px-3 py-1.5 text-sm rounded-lg transition-all',
                selectedCatalogs.includes(catalog)
                  ? 'bg-[var(--color-primary)]/15 text-[var(--color-primary)] border border-[var(--color-primary)]/30'
                  : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] border border-transparent hover:border-[var(--color-border)]'
              )}
            >
              {catalog}
            </button>
          ))}
        </div>

        {/* ChEMBL catalogs — collapsible group */}
        <div className="mt-3 border border-[var(--color-border)] rounded-lg overflow-hidden">
          <button
            onClick={onToggleChemblExpanded}
            className="w-full flex items-center justify-between px-3 py-2 text-sm bg-[var(--color-surface-sunken)] hover:bg-[var(--color-surface-sunken)]/80 transition-colors"
          >
            <div className="flex items-center gap-2">
              <input
                type="checkbox"
                checked={allChemblSelected}
                ref={(el) => { if (el) el.indeterminate = someChemblSelected && !allChemblSelected; }}
                onChange={(e) => { e.stopPropagation(); onToggleAllChembl(e.target.checked); }}
                onClick={(e) => e.stopPropagation()}
                className="rounded border-gray-300 text-[var(--color-primary)] focus:ring-[var(--color-primary)]"
              />
              <span className="font-medium text-[var(--color-text-primary)]">ChEMBL Filters</span>
              {someChemblSelected && (
                <span className="text-xs text-[var(--color-text-muted)]">
                  ({chemblCatalogs.filter((c) => selectedCatalogs.includes(c.id)).length}/{chemblCatalogs.length})
                </span>
              )}
            </div>
            <ChevronDown className={cn('w-4 h-4 text-[var(--color-text-muted)] transition-transform', chemblExpanded && 'rotate-180')} />
          </button>
          {chemblExpanded && (
            <div className="px-3 py-2 flex flex-wrap gap-2 border-t border-[var(--color-border)]">
              {chemblCatalogs.map((catalog) => (
                <button
                  key={catalog.id}
                  onClick={() => onToggleCatalog(catalog.id)}
                  className={cn(
                    'px-3 py-1.5 text-sm rounded-lg transition-all',
                    selectedCatalogs.includes(catalog.id)
                      ? 'bg-[var(--color-primary)]/15 text-[var(--color-primary)] border border-[var(--color-primary)]/30'
                      : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] border border-transparent hover:border-[var(--color-border)]'
                  )}
                >
                  {catalog.label}
                </button>
              ))}
            </div>
          )}
        </div>
      </div>

      <ClayButton
        variant="primary"
        onClick={onScreenAlerts}
        disabled={screenDisabled}
        loading={alertsLoading}
        leftIcon={<AlertTriangle className="w-4 h-4" />}
      >
        Screen Alerts
      </ClayButton>

      {/* Safety Assessment (merged into this tab) */}
      <SafetyAccordion
        smiles={safetySmiles}
        alertResult={safetyAlertResult}
        safetyResult={safetyResult}
        isLoading={safetyLoading}
        error={safetyError}
      />
    </div>
  );
}

interface StandardizeTabPanelProps {
  includeTautomer: boolean;
  onToggleTautomer: (checked: boolean) => void;
  onStandardize: () => void;
  /** !molecule.trim() || isAnyLoading */
  standardizeDisabled: boolean;
  standardizeLoading: boolean;
}

export function StandardizeTabPanel({
  includeTautomer,
  onToggleTautomer,
  onStandardize,
  standardizeDisabled,
  standardizeLoading,
}: StandardizeTabPanelProps) {
  return (
    <div className="space-y-4">
      <div className="flex items-start gap-3">
        <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
          <Info className="w-4 h-4" />
        </div>
        <div>
          <p className="text-[var(--color-text-secondary)] text-sm mb-3">
            Click Standardize to apply the ChEMBL structure pipeline:
          </p>
          <ul className="list-none space-y-1 text-sm text-[var(--color-text-secondary)]">
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)]"></span>
              Salt and solvent removal
            </li>
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)]"></span>
              Charge neutralization
            </li>
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)]"></span>
              Stereochemistry normalization
            </li>
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)]"></span>
              Tautomer canonicalization (optional)
            </li>
          </ul>
        </div>
      </div>

      <label className="flex items-center gap-2 text-sm text-[var(--color-text-secondary)] cursor-pointer">
        <input
          type="checkbox"
          checked={includeTautomer}
          onChange={(e) => onToggleTautomer(e.target.checked)}
          className="rounded border-[var(--color-border-strong)] text-[var(--color-primary)] focus:ring-[var(--color-primary)]/30"
        />
        Enable tautomer canonicalization
      </label>

      <ClayButton
        variant="primary"
        onClick={onStandardize}
        disabled={standardizeDisabled}
        loading={standardizeLoading}
        leftIcon={<Layers className="w-4 h-4" />}
      >
        Standardize
      </ClayButton>
    </div>
  );
}
