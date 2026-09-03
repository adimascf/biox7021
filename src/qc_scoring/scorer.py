from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
from typing import Dict, List, Optional
import pandas as pd
from qc_scoring.criteria import calculate_cohort_criteria, CohortCriteriaResult
from qc_scoring.models import Scenario
from qc_scoring.preferences import GateConfig, WeightsConfig, validate_weights
from qc_scoring.validation import filter_scenario_data

SCORING_VERSION = "1.0"
SOURCE_DATA_BASELINE_COMMIT = "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562"


@dataclass
class Provenance:
    scoring_version: str
    source_data_hash: str
    source_data_commit: Optional[str]
    scenario_model: str
    scenario_depth: str
    weights: Dict[str, float]
    gates: Dict[str, bool]
    generated_at: str


@dataclass
class CombinationRecommendation:
    combo: str
    rank: Optional[int]
    is_eligible: bool
    ineligible_reason: Optional[str]
    is_near_tie: bool
    overall_score: float
    display_score: Optional[float]

    # 4 Criterion scores
    score_contiguity: float
    score_accuracy: float
    score_residual: float
    score_replicon: float

    # Raw evidence summaries
    mean_auNGA_ratio: float
    mean_error_rate: float
    mean_mismatches: float
    mean_indels: float

    residual_total_hits: int
    residual_clean_isolates: int
    residual_affected_isolates: int

    replicon_total_missed: int
    replicon_full_missed: int
    replicon_partial_missed: int
    replicon_affected_isolates: int

    warnings: List[str] = field(default_factory=list)
    isolate_contiguity_scores: Dict[str, float] = field(default_factory=dict)
    isolate_accuracy_scores: Dict[str, float] = field(default_factory=dict)


@dataclass
class ScoringResult:
    scoring_version: str
    scenario: Scenario
    weights: WeightsConfig
    gates: GateConfig
    recommendations: List[CombinationRecommendation]
    provenance: Provenance

    @property
    def eligible_count(self) -> int:
        return sum(1 for r in self.recommendations if r.is_eligible)

    @property
    def ineligible_count(self) -> int:
        return sum(1 for r in self.recommendations if not r.is_eligible)

    @property
    def leading_recommendation(self) -> Optional[CombinationRecommendation]:
        eligible = [r for r in self.recommendations if r.is_eligible]
        return eligible[0] if eligible else None


def _calculate_content_hash(df: pd.DataFrame) -> str:
    # Deterministic content hash of the dataframe
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    return hashlib.sha256(csv_bytes).hexdigest()


def score_benchmark(
    df: pd.DataFrame,
    scenario: Scenario,
    weights: WeightsConfig,
    gates: Optional[GateConfig] = None,
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
) -> ScoringResult:
    if gates is None:
        gates = GateConfig()

    validate_weights(weights)
    filtered_df, incomplete_reasons = filter_scenario_data(df, scenario)

    eligible_recs: List[CombinationRecommendation] = []
    ineligible_recs: List[CombinationRecommendation] = []

    # Process complete combinations
    grouped = filtered_df.groupby("combo", observed=True)
    for combo, combo_df in grouped:
        combo_name = str(combo)
        cohort: CohortCriteriaResult = calculate_cohort_criteria(combo_df)

        # Preference alignment score (unrounded full precision)
        overall = (
            (weights.accuracy * cohort.score_accuracy)
            + (weights.contiguity * cohort.score_contiguity)
            + (weights.residual * cohort.score_residual)
            + (weights.replicon * cohort.score_replicon)
        ) / 100.0

        # Evaluate eligibility gates
        is_eligible = True
        fail_reasons: List[str] = []

        if gates.complete_recovery and not cohort.all_replicons_complete:
            is_eligible = False
            fail_reasons.append("fails complete-recovery gate (<95% coverage on at least one replicon)")

        if gates.zero_residual_hits and not cohort.zero_residual_hits:
            is_eligible = False
            fail_reasons.append(
                f"fails zero-residual-hits gate ({cohort.residual_total_hits} residual hit(s) detected)"
            )

        # Generate concrete warnings
        warnings: List[str] = []
        if cohort.residual_total_hits > 0:
            warnings.append(
                f"{cohort.residual_total_hits} residual adapter/barcode hit(s) across "
                f"{cohort.residual_affected_isolates} isolate(s)"
            )

        if cohort.replicon_total_missed > 0:
            warnings.append(
                f"{cohort.replicon_total_missed} missed or severely incomplete replicon(s) across "
                f"{cohort.replicon_affected_isolates} isolate(s) "
                f"({cohort.replicon_full_missed} full, {cohort.replicon_partial_missed} partial)"
            )

        # Isolate variability warnings (>= 25 points between mean and min isolate)
        min_acc = min(cohort.isolate_accuracy_scores.values()) if cohort.isolate_accuracy_scores else cohort.score_accuracy
        if (cohort.score_accuracy - min_acc) >= 25.0:
            warnings.append("Variable sequence accuracy across isolates (lowest isolate score is >=25 points below cohort mean)")

        min_cont = min(cohort.isolate_contiguity_scores.values()) if cohort.isolate_contiguity_scores else cohort.score_contiguity
        if (cohort.score_contiguity - min_cont) >= 25.0:
            warnings.append("Variable contiguity across isolates (lowest isolate score is >=25 points below cohort mean)")

        rec = CombinationRecommendation(
            combo=combo_name,
            rank=None,
            is_eligible=is_eligible,
            ineligible_reason="; ".join(fail_reasons) if fail_reasons else None,
            is_near_tie=False,
            overall_score=overall,
            display_score=round(overall, 1) if is_eligible else round(overall, 1),
            score_contiguity=cohort.score_contiguity,
            score_accuracy=cohort.score_accuracy,
            score_residual=cohort.score_residual,
            score_replicon=cohort.score_replicon,
            mean_auNGA_ratio=cohort.mean_auNGA_ratio,
            mean_error_rate=cohort.mean_error_rate,
            mean_mismatches=cohort.mean_mismatches,
            mean_indels=cohort.mean_indels,
            residual_total_hits=cohort.residual_total_hits,
            residual_clean_isolates=cohort.residual_clean_isolates,
            residual_affected_isolates=cohort.residual_affected_isolates,
            replicon_total_missed=cohort.replicon_total_missed,
            replicon_full_missed=cohort.replicon_full_missed,
            replicon_partial_missed=cohort.replicon_partial_missed,
            replicon_affected_isolates=cohort.replicon_affected_isolates,
            warnings=warnings,
            isolate_contiguity_scores=cohort.isolate_contiguity_scores,
            isolate_accuracy_scores=cohort.isolate_accuracy_scores,
        )

        if is_eligible:
            eligible_recs.append(rec)
        else:
            ineligible_recs.append(rec)

    # Process combinations with insufficient data
    for combo_name, reason in incomplete_reasons.items():
        rec = CombinationRecommendation(
            combo=combo_name,
            rank=None,
            is_eligible=False,
            ineligible_reason=reason,
            is_near_tie=False,
            overall_score=0.0,
            display_score=None,
            score_contiguity=0.0,
            score_accuracy=0.0,
            score_residual=0.0,
            score_replicon=0.0,
            mean_auNGA_ratio=0.0,
            mean_error_rate=0.0,
            mean_mismatches=0.0,
            mean_indels=0.0,
            residual_total_hits=0,
            residual_clean_isolates=0,
            residual_affected_isolates=0,
            replicon_total_missed=0,
            replicon_full_missed=0,
            replicon_partial_missed=0,
            replicon_affected_isolates=0,
            warnings=["Insufficient benchmark data"],
        )
        ineligible_recs.append(rec)

    # Rank eligible combinations: higher overall_score first, stable alphabetical combo sort
    eligible_recs.sort(key=lambda r: (-r.overall_score, r.combo))

    # Competition ranking for exact ties (1, 1, 3, etc.)
    for idx, rec in enumerate(eligible_recs):
        if idx == 0:
            rec.rank = 1
        else:
            prev = eligible_recs[idx - 1]
            if abs(rec.overall_score - prev.overall_score) < 1e-9:
                rec.rank = prev.rank
            else:
                rec.rank = idx + 1

    # Near-tie calculation (< 1.0 point difference between adjacent eligible combinations)
    for idx in range(len(eligible_recs) - 1):
        diff = abs(eligible_recs[idx].overall_score - eligible_recs[idx + 1].overall_score)
        if diff < 1.0:
            eligible_recs[idx].is_near_tie = True
            eligible_recs[idx + 1].is_near_tie = True

    # Order ineligible combinations by overall score descending (or combo name)
    ineligible_recs.sort(key=lambda r: (-r.overall_score, r.combo))

    all_recommendations = eligible_recs + ineligible_recs

    now_iso = datetime.now(timezone.utc).isoformat()
    content_hash = _calculate_content_hash(df)

    provenance = Provenance(
        scoring_version=SCORING_VERSION,
        source_data_hash=content_hash,
        source_data_commit=source_data_commit,
        scenario_model=scenario.model,
        scenario_depth=scenario.depth,
        weights=weights.as_dict(),
        gates=gates.as_dict(),
        generated_at=now_iso,
    )

    return ScoringResult(
        scoring_version=SCORING_VERSION,
        scenario=scenario,
        weights=weights,
        gates=gates,
        recommendations=all_recommendations,
        provenance=provenance,
    )
