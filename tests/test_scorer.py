import pytest
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.preferences import (
    WeightsConfig,
    GateConfig,
    PresetName,
    get_preset,
    validate_weights,
    InvalidWeightsError,
)
from qc_scoring.scorer import (
    score_benchmark,
    ScoringResult,
    CombinationRecommendation,
)
from qc_scoring.validation import EXPECTED_SAMPLES

def test_presets_weights_and_gates():
    # Community-balanced
    cb = get_preset(PresetName.COMMUNITY_BALANCED)
    assert cb.weights == WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0)
    assert cb.gates == GateConfig(complete_recovery=False, zero_residual_hits=False)

    # Complete-replicon-recovery
    cr = get_preset(PresetName.COMPLETE_REPLICON_RECOVERY)
    assert cr.weights == WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0)
    assert cr.gates == GateConfig(complete_recovery=True, zero_residual_hits=False)

    # Sequence-accurate-assembly
    sa = get_preset(PresetName.SEQUENCE_ACCURATE_ASSEMBLY)
    assert sa.weights == WeightsConfig(accuracy=50.0, contiguity=14.0, residual=12.0, replicon=24.0)
    assert sa.gates == GateConfig(complete_recovery=False, zero_residual_hits=False)

def test_weight_validation_rules():
    # Valid
    validate_weights(WeightsConfig(accuracy=25.0, contiguity=25.0, residual=25.0, replicon=25.0))
    validate_weights(WeightsConfig(accuracy=100.0, contiguity=0.0, residual=0.0, replicon=0.0))
    
    # Negative weight
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=-5.0, contiguity=35.0, residual=35.0, replicon=35.0))
        
    # All zero
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=0.0, contiguity=0.0, residual=0.0, replicon=0.0))
        
    # Total below 100
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=20.0, contiguity=20.0, residual=20.0, replicon=20.0))
        
    # Total above 100
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=30.0, contiguity=30.0, residual=30.0, replicon=30.0))

def test_all_ineligible_scenario():
    # Load benchmark data and force strict gates that cause every combination to fail
    df = pd.read_csv("logbook/assembly_metrics.csv")
    scenario = Scenario(model="hac", depth="20x")
    # Complete recovery gate + zero residual hits gate
    result = score_benchmark(
        df,
        scenario,
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(complete_recovery=True, zero_residual_hits=True)
    )
    # Check that even if combinations fail, gates were not relaxed
    assert result.scenario == scenario
    for rec in result.recommendations:
        if not rec.is_eligible:
            assert rec.rank is None
            assert rec.ineligible_reason is not None

def test_near_ties_and_exact_ties():
    # Construct scenario with exact tie and near tie using expected isolate names
    sample_list = sorted(EXPECTED_SAMPLES)
    records = []
    for c in ["chopper-barbell", "chopper-dorado", "chopper-untrimmed"]:
        for sample in sample_list:
            if c in ("chopper-barbell", "chopper-dorado"):
                aunga = 1.0
                err = 0.0
            else:  # chopper-untrimmed slightly worse
                aunga = 1.0
                err = 0.05  # rate 0.05 -> score 99.5 instead of 100.0 (diff 0.5 < 1.0)
            records.append({
                "combo": c,
                "depth": "100x",
                "sample": sample,
                "model": "hac",
                "Mismatches per 100kbp": err,
                "Indels per 100kbp": 0.0,
                "auNGA_ratio": aunga,
                "contamination_count": 0,
                "full_missed": 0,
                "partial_missed": 0,
                "total_missed": 0,
                "all_contigs_coverage": "chr (1000000bp, 100.0% cov)",
            })
    df_tie = pd.DataFrame(records)
    res = score_benchmark(
        df_tie,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=50.0, contiguity=50.0, residual=0.0, replicon=0.0),
        gates=GateConfig()
    )
    # Filter to the 3 evaluated combos in eligible list
    eval_recs = [r for r in res.recommendations if r.combo in ("chopper-barbell", "chopper-dorado", "chopper-untrimmed")]
    # chopper-barbell and chopper-dorado should have exact same rank 1
    assert eval_recs[0].combo == "chopper-barbell"
    assert eval_recs[0].rank == 1
    assert eval_recs[1].combo == "chopper-dorado"
    assert eval_recs[1].rank == 1
    
    # chopper-untrimmed should have rank 3
    assert eval_recs[2].combo == "chopper-untrimmed"
    assert eval_recs[2].rank == 3
    
    # Near tie flags: adjacent difference is 0.25 (< 1.0)
    assert eval_recs[1].is_near_tie is True
    assert eval_recs[2].is_near_tie is True

def test_warning_25_point_gap():
    # Test accuracy variability warning: mean - min >= 25.0
    sample_list = sorted(EXPECTED_SAMPLES)
    records = []
    for idx, sample in enumerate(sample_list):
        # 12 isolates score 100 (err=0), 1 isolate scores 70 (err=3.0) -> mean = 97.69
        # min = 70.0. Gap = 97.69 - 70.0 = 27.69 >= 25.0 -> should warn!
        err = 3.0 if idx == 0 else 0.0
        records.append({
            "combo": "chopper-barbell",
            "depth": "100x",
            "sample": sample,
            "model": "hac",
            "Mismatches per 100kbp": err,
            "Indels per 100kbp": 0.0,
            "auNGA_ratio": 1.0,
            "contamination_count": 0,
            "full_missed": 0,
            "partial_missed": 0,
            "total_missed": 0,
            "all_contigs_coverage": "chr (1000000bp, 100.0% cov)",
        })
    df_warn = pd.DataFrame(records)
    res = score_benchmark(
        df_warn,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=100.0, contiguity=0.0, residual=0.0, replicon=0.0),
        gates=GateConfig()
    )
    rec = [r for r in res.recommendations if r.combo == "chopper-barbell"][0]
    has_var_warn = any("Variable sequence accuracy" in w for w in rec.warnings)
    assert has_var_warn is True

def test_insufficient_data_missing_isolate():
    # Omit 1 isolate for chopper-barbell
    df = pd.read_csv("logbook/assembly_metrics.csv")
    mask = ~((df["combo"] == "chopper-barbell") & (df["sample"] == "AJ292__202310"))
    df_missing = df[mask].copy()

    res = score_benchmark(
        df_missing,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(),
    )
    rec = [r for r in res.recommendations if r.combo == "chopper-barbell"][0]
    assert rec.is_eligible is False
    assert rec.rank is None
    assert "insufficient benchmark data" in (rec.ineligible_reason or "")
    assert "AJ292__202310" in (rec.ineligible_reason or "")

def test_pinned_real_data_acceptance_fixture():
    df = pd.read_csv("logbook/assembly_metrics.csv")
    cb_preset = get_preset(PresetName.COMMUNITY_BALANCED)

    for model in ["hac", "sup"]:
        for depth in ["20x", "100x"]:
            sc = Scenario(model=model, depth=depth)
            result = score_benchmark(
                df,
                sc,
                weights=cb_preset.weights,
                gates=cb_preset.gates,
            )
            assert len(result.recommendations) == 17
            assert result.scoring_version == "1.0"
            assert result.scenario == sc
            assert result.provenance is not None
            # Every recommendation must have 4 criterion scores in [0, 100]
            for rec in result.recommendations:
                assert 0.0 <= rec.score_contiguity <= 100.0
                assert 0.0 <= rec.score_accuracy <= 100.0
                assert 0.0 <= rec.score_residual <= 100.0
                assert 0.0 <= rec.score_replicon <= 100.0
                if rec.is_eligible:
                    assert rec.rank is not None
                    assert 0.0 <= rec.overall_score <= 100.0
                    assert rec.display_score == round(rec.overall_score, 1)

def test_export_recommendations_csv():
    from qc_scoring.scorer import export_recommendations_csv, recommendations_to_dataframe
    df = pd.read_csv("logbook/assembly_metrics.csv")
    res = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(),
    )
    csv_text = export_recommendations_csv(res)
    assert len(csv_text) > 0
    lines = csv_text.strip().split("\n")
    assert len(lines) == 18  # 1 header + 17 rows
    export_df = recommendations_to_dataframe(res)
    assert len(export_df) == 17
    assert "rank" in export_df.columns
    assert "combo" in export_df.columns
    assert "overall_score" in export_df.columns
    assert "scoring_version" in export_df.columns
    assert "source_data_commit" in export_df.columns


def test_complete_replicon_recovery_preset_and_gate():
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    preset = get_preset(PresetName.COMPLETE_REPLICON_RECOVERY)
    assert preset.weights == WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0)
    assert preset.gates == GateConfig(complete_recovery=True, zero_residual_hits=False)

    res = score_benchmark(df, Scenario(model="hac", depth="100x"), weights=preset.weights, gates=preset.gates)
    assert res.eligible_count == 7
    assert res.ineligible_count == 10
    leader = res.leading_recommendation
    assert leader is not None
    assert leader.combo == "seqkit-barbell"
    assert leader.rank == 1
    assert leader.display_score == 98.0

    # Ineligible combinations must receive rank=None and exact failure reasons
    for rec in res.recommendations:
        if not rec.is_eligible:
            assert rec.rank is None
            assert "fails complete-recovery gate (<95% coverage on at least one replicon)" in (rec.ineligible_reason or "")


def test_sequence_accurate_assembly_preset_and_ranking():
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    preset = get_preset(PresetName.SEQUENCE_ACCURATE_ASSEMBLY)
    assert preset.weights == WeightsConfig(accuracy=50.0, contiguity=14.0, residual=12.0, replicon=24.0)
    assert preset.gates == GateConfig(complete_recovery=False, zero_residual_hits=False)

    res = score_benchmark(df, Scenario(model="hac", depth="100x"), weights=preset.weights, gates=preset.gates)
    assert res.eligible_count == 17
    assert res.ineligible_count == 0
    leader = res.leading_recommendation
    assert leader is not None
    assert leader.combo == "seqkit-barbell"
    assert leader.rank == 1
    assert leader.display_score == 97.6


def test_custom_decimal_percentages_and_edge_cases():
    # Valid decimal percentages totaling exactly 100.0
    w_dec = WeightsConfig(accuracy=33.3, contiguity=33.3, residual=16.7, replicon=16.7)
    validate_weights(w_dec)

    # Valid with single non-zero weight
    w_single = WeightsConfig(accuracy=100.0, contiguity=0.0, residual=0.0, replicon=0.0)
    validate_weights(w_single)

    # Valid with zeros in other fields
    w_zero = WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0)
    validate_weights(w_zero)

    # Invalid: total is 99.9 (fails exact 100)
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=33.3, contiguity=33.3, residual=16.6, replicon=16.7))

    # Invalid: total is 100.1
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=33.3, contiguity=33.4, residual=16.7, replicon=16.7))

    # Invalid: negative decimal
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=-0.5, contiguity=50.5, residual=25.0, replicon=25.0))


def test_multiple_exclusion_reasons_and_reasons_format():
    # Load benchmark data and apply both gates
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    res = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(complete_recovery=True, zero_residual_hits=True),
    )
    # Find unprocessed-dorado which has both residual hits and sub-95% replicons
    rec = [r for r in res.recommendations if r.combo == "unprocessed-dorado"][0]
    assert rec.is_eligible is False
    assert rec.rank is None
    reason = rec.ineligible_reason or ""
    assert "fails complete-recovery gate (<95% coverage on at least one replicon)" in reason
    assert "fails zero-residual-hits gate" in reason


def test_complete_recovery_sub_95_gate_edge_case():
    # Construct a case where total_missed == 0 (no replicon < 50%), but one replicon has 94.0% coverage (< 95%)
    sample_list = sorted(EXPECTED_SAMPLES)
    records = []
    for idx, sample in enumerate(sample_list):
        # 12 isolates have 100% cov, 1 isolate has 94.0% cov on plasmid
        cov_str = "chr (1000000bp, 100.0% cov); plas (50000bp, 94.0% cov)" if idx == 0 else "chr (1000000bp, 100.0% cov)"
        records.append({
            "combo": "seqkit-dorado",
            "depth": "100x",
            "sample": sample,
            "model": "hac",
            "Mismatches per 100kbp": 1.0,
            "Indels per 100kbp": 0.0,
            "auNGA_ratio": 1.0,
            "contamination_count": 0,
            "full_missed": 0,
            "partial_missed": 0,
            "total_missed": 0,
            "all_contigs_coverage": cov_str,
        })
    df_edge = pd.DataFrame(records)
    res = score_benchmark(
        df_edge,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0),
        gates=GateConfig(complete_recovery=True, zero_residual_hits=False),
    )
    rec = [r for r in res.recommendations if r.combo == "seqkit-dorado"][0]
    # Even though total_missed == 0, 94.0% is < 95.0%, so it must FAIL the complete recovery gate!
    assert rec.is_eligible is False
    assert rec.rank is None
    assert "fails complete-recovery gate (<95% coverage on at least one replicon)" in (rec.ineligible_reason or "")

