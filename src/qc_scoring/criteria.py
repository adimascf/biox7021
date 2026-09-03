from dataclasses import dataclass, field
from typing import Dict, List
import numpy as np
import pandas as pd
from qc_scoring.parser import validate_and_parse_replicon_row

def score_isolate_contiguity(auNGA_ratio: float) -> float:
    return float(max(0.0, 1.0 - abs(auNGA_ratio - 1.0)) * 100.0)

def score_isolate_accuracy(mismatches_per_100kbp: float, indels_per_100kbp: float) -> float:
    event_rate = mismatches_per_100kbp + indels_per_100kbp
    val = 1.0 - (event_rate / 10.0)
    clipped = min(1.0, max(0.0, val))
    return float(clipped * 100.0)

def score_isolate_residual_clean(contamination_count: int | float) -> float:
    return 100.0 if contamination_count == 0 else 0.0

def score_replicon_recovery(total_missed_cohort: int | float) -> float:
    val = 1.0 - (total_missed_cohort / 3.0)
    return float(max(0.0, val) * 100.0)

@dataclass
class CohortCriteriaResult:
    # 4 fixed 0-100 criterion scores
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

    # Gate evaluations
    all_replicons_complete: bool
    zero_residual_hits: bool

    # Per-isolate details (for warnings and visualisations)
    isolate_contiguity_scores: Dict[str, float] = field(default_factory=dict)
    isolate_accuracy_scores: Dict[str, float] = field(default_factory=dict)

def calculate_cohort_criteria(combo_df: pd.DataFrame) -> CohortCriteriaResult:
    n_isolates = len(combo_df)
    if n_isolates == 0:
        raise ValueError("Cannot calculate cohort criteria on empty DataFrame")

    contiguity_scores: Dict[str, float] = {}
    accuracy_scores: Dict[str, float] = {}
    clean_count = 0
    total_hits = 0
    affected_residual_count = 0

    full_missed_sum = 0
    partial_missed_sum = 0
    total_missed_sum = 0
    affected_replicon_count = 0

    all_complete_recovery = True

    mismatches_sum = 0.0
    indels_sum = 0.0
    aunga_sum = 0.0

    for _, row in combo_df.iterrows():
        sample_name = str(row["sample"])

        # Contiguity
        aunga = float(row["auNGA_ratio"])
        aunga_sum += aunga
        cont_score = score_isolate_contiguity(aunga)
        contiguity_scores[sample_name] = cont_score

        # Accuracy
        mism = float(row["Mismatches per 100kbp"])
        indel = float(row["Indels per 100kbp"])
        mismatches_sum += mism
        indels_sum += indel
        acc_score = score_isolate_accuracy(mism, indel)
        accuracy_scores[sample_name] = acc_score

        # Residual hits
        hits = int(row["contamination_count"])
        total_hits += hits
        if hits == 0:
            clean_count += 1
        else:
            affected_residual_count += 1

        # Replicons
        f_mis = int(row["full_missed"])
        p_mis = int(row["partial_missed"])
        t_mis = int(row["total_missed"])
        full_missed_sum += f_mis
        partial_missed_sum += p_mis
        total_missed_sum += t_mis
        if t_mis > 0:
            affected_replicon_count += 1

        # Parse replicons and verify completeness gate (>= 95% for all replicons)
        _, is_complete = validate_and_parse_replicon_row(dict(row))
        if not is_complete:
            all_complete_recovery = False

    # Cohort scores
    score_contiguity = float(np.mean(list(contiguity_scores.values())))
    score_accuracy = float(np.mean(list(accuracy_scores.values())))
    score_residual = float((clean_count / n_isolates) * 100.0)
    score_replicon = score_replicon_recovery(total_missed_sum)

    return CohortCriteriaResult(
        score_contiguity=score_contiguity,
        score_accuracy=score_accuracy,
        score_residual=score_residual,
        score_replicon=score_replicon,
        mean_auNGA_ratio=aunga_sum / n_isolates,
        mean_error_rate=(mismatches_sum + indels_sum) / n_isolates,
        mean_mismatches=mismatches_sum / n_isolates,
        mean_indels=indels_sum / n_isolates,
        residual_total_hits=total_hits,
        residual_clean_isolates=clean_count,
        residual_affected_isolates=affected_residual_count,
        replicon_total_missed=total_missed_sum,
        replicon_full_missed=full_missed_sum,
        replicon_partial_missed=partial_missed_sum,
        replicon_affected_isolates=affected_replicon_count,
        all_replicons_complete=all_complete_recovery,
        zero_residual_hits=(total_hits == 0),
        isolate_contiguity_scores=contiguity_scores,
        isolate_accuracy_scores=accuracy_scores,
    )
