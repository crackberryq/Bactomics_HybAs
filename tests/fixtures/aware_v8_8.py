# Frozen controller from b14ff7688b72f1ae4a39bd99a4305d16a0e53ec1.
import os

def normalize_taxid(x):
    s = str(x or "").strip()
    if not s:
        return ""
    if s.endswith(".0"):
        s = s[:-2]
    return s


def safe_float(x, default=0.0):
    try:
        return float(str(x).strip())
    except Exception:
        return default


def combine_platform_values(val_illumina, val_ont, use_both=True):
    vals = []
    try:
        if val_illumina is not None:
            vals.append(("illumina", float(val_illumina)))
    except Exception:
        pass
    try:
        if val_ont is not None:
            vals.append(("ont", float(val_ont)))
    except Exception:
        pass

    # keep only meaningful non-negative values from present platforms
    vals = [(name, v) for name, v in vals if v >= 0]

    if not vals:
        return 0.0, "no_platform_values"

    if len(vals) == 1:
        return vals[0][1], f"single_platform_{vals[0][0]}"

    if use_both:
        return min(v for _, v in vals), "both_platforms_min"
    return max(v for _, v in vals), "both_platforms_max"


def choose_platform_cleaning_decision(mode, selected_taxid, platform_pct, platform_name):
    selected_taxid = normalize_taxid(selected_taxid)

    if mode == "off":
        return ("passthrough", f"{platform_name}_mode_off")

    if mode == "manual":
        if not selected_taxid:
            return ("passthrough", f"{platform_name}_manual_mode_no_target_taxid")
        return ("whitelist_filter", f"{platform_name}_manual_mode_target_taxid_supplied")

    # auto / aware: decide per platform
    if not selected_taxid:
        return ("passthrough", f"{platform_name}_{mode}_mode_no_selected_target_taxid")

    try:
        pct = float(platform_pct)
    except Exception:
        pct = 0.0

    if pct >= KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:
        return (
            "passthrough",
            f"{platform_name}_{mode}_mode_selected_target_pct={pct:.2f}; "
            f"threshold={KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:.2f}"
        )

    return (
        "whitelist_filter",
        f"{platform_name}_{mode}_mode_selected_target_pct={pct:.2f}; "
        f"threshold={KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:.2f}"
    )



def parse_kraken_report_rows(report_path):
    rows = []
    if not os.path.exists(report_path):
        return rows

    with open(report_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rows.append({
                "pct": safe_float(parts[0], 0.0),
                "clade_reads": parts[1].strip(),
                "direct_reads": parts[2].strip(),
                "rank": parts[3].strip(),
                "taxid": normalize_taxid(parts[4]),
                "name": parts[5].strip(),
            })
    return rows


def parse_kraken_report_target_fraction(report_path, target_taxid):
    target_taxid = normalize_taxid(target_taxid)
    if not target_taxid:
        return 0.0

    for row in parse_kraken_report_rows(report_path):
        if row["taxid"] == target_taxid:
            return row["pct"]
    return 0.0


def is_meaningful_kraken_candidate(row):
    rank = row["rank"]
    name = row["name"].strip().lower()

    if not row["taxid"] or row["taxid"] == "0":
        return False
    if name.startswith("unclassified"):
        return False
    if rank in {"U", "R", "R1", "R2", "K", "P", "C", "O", "F"}:
        return False

    # allow only genus/species level primary nodes for auto target discovery
    return rank in {"G", "S"}


def infer_target_from_kraken_report(report_path):
    rows = [r for r in parse_kraken_report_rows(report_path) if is_meaningful_kraken_candidate(r)]

    if not rows:
        return {
            "taxid": "",
            "name": "NA",
            "rank": "NA",
            "pct": 0.0,
            "source_reason": "no_meaningful_candidate"
        }

    best_species = None
    best_genus = None

    for row in rows:
        if row["rank"] == "S":
            if best_species is None or row["pct"] > best_species["pct"]:
                best_species = row
        elif row["rank"] == "G":
            if best_genus is None or row["pct"] > best_genus["pct"]:
                best_genus = row

    if best_species is not None and best_species["pct"] >= KRAKEN_AUTO_PREFER_SPECIES_MIN_PCT:
        chosen = best_species
        why = f"species_preferred_at_{best_species['pct']:.2f}_percent"
    elif best_genus is not None:
        chosen = best_genus
        why = f"genus_preferred_at_{best_genus['pct']:.2f}_percent"
    elif best_species is not None:
        chosen = best_species
        why = f"species_fallback_at_{best_species['pct']:.2f}_percent"
    else:
        return {
            "taxid": "",
            "name": "NA",
            "rank": "NA",
            "pct": 0.0,
            "source_reason": "no_candidate_after_selection"
        }

    return {
        "taxid": chosen["taxid"],
        "name": chosen["name"],
        "rank": chosen["rank"],
        "pct": chosen["pct"],
        "source_reason": why
    }


def choose_selected_kraken_target(mode, manual_target_taxid, illumina_report, ont_report):
    manual_target_taxid = normalize_taxid(manual_target_taxid)

    if mode == "off":
        return {
            "selected_taxid": "",
            "selected_name": "NA",
            "selected_rank": "NA",
            "selected_source": "off_mode",
            "illumina_selected_pct": 0.0,
            "ont_selected_pct": 0.0,
            "selected_pct_for_decision": 0.0,
            "inferred_taxid": "",
            "inferred_name": "NA",
            "inferred_rank": "NA",
            "inferred_pct": 0.0,
            "selection_reason": "kraken_mode_off"
        }

    if mode == "manual":
        if not manual_target_taxid:
            return {
                "selected_taxid": "",
                "selected_name": "NA",
                "selected_rank": "NA",
                "selected_source": "manual_empty",
                "illumina_selected_pct": 0.0,
                "ont_selected_pct": 0.0,
                "selected_pct_for_decision": 0.0,
                "inferred_taxid": "",
                "inferred_name": "NA",
                "inferred_rank": "NA",
                "inferred_pct": 0.0,
                "selection_reason": "manual_mode_no_target_taxid"
            }

        illumina_pct = parse_kraken_report_target_fraction(illumina_report, manual_target_taxid)
        ont_pct = parse_kraken_report_target_fraction(ont_report, manual_target_taxid)

        metric, metric_mode = combine_platform_values(
            illumina_pct if os.path.getsize(illumina_report) > 0 else None,
            ont_pct if os.path.getsize(ont_report) > 0 else None,
            use_both=KRAKEN_AUTO_USE_BOTH_PLATFORMS,
        )

        return {
            "selected_taxid": manual_target_taxid,
            "selected_name": "MANUAL_TARGET",
            "selected_rank": "MANUAL",
            "selected_source": "manual_target_taxid",
            "illumina_selected_pct": illumina_pct,
            "ont_selected_pct": ont_pct,
            "selected_pct_for_decision": metric,
            "inferred_taxid": "",
            "inferred_name": "NA",
            "inferred_rank": "NA",
            "inferred_pct": 0.0,
            "selection_reason": "manual_mode_target_taxid_supplied"
        }

    # auto and aware share the same target-selection logic:
    # 1) manual override if TSV taxid exists
    # 2) otherwise infer from Kraken
    if mode in {"auto", "aware"}:
        if manual_target_taxid:
            illumina_pct = parse_kraken_report_target_fraction(illumina_report, manual_target_taxid)
            ont_pct = parse_kraken_report_target_fraction(ont_report, manual_target_taxid)
            metric, metric_mode = combine_platform_values(
                illumina_pct if os.path.getsize(illumina_report) > 0 else None,
                ont_pct if os.path.getsize(ont_report) > 0 else None,
                use_both=KRAKEN_AUTO_USE_BOTH_PLATFORMS,
            )

            return {
                "selected_taxid": manual_target_taxid,
                "selected_name": "MANUAL_TARGET_IN_AUTO",
                "selected_rank": "MANUAL",
                "selected_source": "manual_override_in_auto" if mode == "auto" else "manual_override_in_aware",
                "illumina_selected_pct": illumina_pct,
                "ont_selected_pct": ont_pct,
                "selected_pct_for_decision": metric,
                "inferred_taxid": "",
                "inferred_name": "NA",
                "inferred_rank": "NA",
                "inferred_pct": 0.0,
                "selection_reason": "auto_mode_manual_override_used" if mode == "auto" else "aware_mode_manual_override_used"
            }

        illumina_inferred = infer_target_from_kraken_report(illumina_report)
        ont_inferred = infer_target_from_kraken_report(ont_report)

        chosen = illumina_inferred if illumina_inferred["pct"] >= ont_inferred["pct"] else ont_inferred
        selected_taxid = normalize_taxid(chosen["taxid"])

        if not selected_taxid:
            return {
                "selected_taxid": "",
                "selected_name": "NA",
                "selected_rank": "NA",
                "selected_source": "auto_inference_failed" if mode == "auto" else "aware_inference_failed",
                "illumina_selected_pct": 0.0,
                "ont_selected_pct": 0.0,
                "selected_pct_for_decision": 0.0,
                "inferred_taxid": "",
                "inferred_name": "NA",
                "inferred_rank": "NA",
                "inferred_pct": 0.0,
                "selection_reason": "auto_mode_no_inferred_target" if mode == "auto" else "aware_mode_no_inferred_target"
            }

        illumina_pct = parse_kraken_report_target_fraction(illumina_report, selected_taxid)
        ont_pct = parse_kraken_report_target_fraction(ont_report, selected_taxid)
        metric = min(illumina_pct, ont_pct) if KRAKEN_AUTO_USE_BOTH_PLATFORMS else max(illumina_pct, ont_pct)

        return {
            "selected_taxid": selected_taxid,
            "selected_name": chosen["name"],
            "selected_rank": chosen["rank"],
            "selected_source": "inferred_from_kraken",
            "illumina_selected_pct": illumina_pct,
            "ont_selected_pct": ont_pct,
            "selected_pct_for_decision": metric,
            "inferred_taxid": chosen["taxid"],
            "inferred_name": chosen["name"],
            "inferred_rank": chosen["rank"],
            "inferred_pct": chosen["pct"],
            "selection_reason": chosen["source_reason"]
        }

    return {
        "selected_taxid": "",
        "selected_name": "NA",
        "selected_rank": "NA",
        "selected_source": "unhandled_mode",
        "illumina_selected_pct": 0.0,
        "ont_selected_pct": 0.0,
        "selected_pct_for_decision": 0.0,
        "inferred_taxid": "",
        "inferred_name": "NA",
        "inferred_rank": "NA",
        "inferred_pct": 0.0,
        "selection_reason": f"unhandled_mode_{mode}"
    }


def choose_kraken_decision(mode, selected_taxid, selected_pct_for_decision):
    selected_taxid = normalize_taxid(selected_taxid)

    if mode == "off":
        return ("passthrough", "kraken_mode_off")

    if mode == "manual":
        if not selected_taxid:
            return ("passthrough", "manual_mode_no_target_taxid")
        return ("whitelist_filter", "manual_mode_target_taxid_supplied")

    # auto / aware base decision
    if not selected_taxid:
        return ("passthrough", f"{mode}_mode_no_selected_target_taxid")

    if selected_pct_for_decision >= KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:
        return (
            "passthrough",
            f"{mode}_mode_selected_target_pct={selected_pct_for_decision:.2f}; "
            f"threshold={KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:.2f}"
        )

    return (
        "whitelist_filter",
        f"{mode}_mode_selected_target_pct={selected_pct_for_decision:.2f}; "
        f"threshold={KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:.2f}"
    )


def top_two_meaningful_candidates(report_path):
    rows = [r for r in parse_kraken_report_rows(report_path) if is_meaningful_kraken_candidate(r)]
    rows = sorted(rows, key=lambda x: x["pct"], reverse=True)
    top1 = rows[0] if len(rows) >= 1 else None
    top2 = rows[1] if len(rows) >= 2 else None
    return top1, top2


def compute_platform_agreement(illumina_report, ont_report):
    illumina_present = os.path.exists(illumina_report) and os.path.getsize(illumina_report) > 0
    ont_present = os.path.exists(ont_report) and os.path.getsize(ont_report) > 0

    if illumina_present and not ont_present:
        return 1.0, "single_platform_illumina_only"
    if ont_present and not illumina_present:
        return 1.0, "single_platform_ont_only"
    if not illumina_present and not ont_present:
        return 0.0, "no_inference"

    i = infer_target_from_kraken_report(illumina_report)
    o = infer_target_from_kraken_report(ont_report)

    if not i["taxid"] or not o["taxid"]:
        return 0.0, "no_inference"

    if i["taxid"] == o["taxid"]:
        return 1.0, "exact_taxid_match"

    if i["rank"] == "S" and o["rank"] == "S":
        i_name = i["name"].split()[0] if i["name"] else ""
        o_name = o["name"].split()[0] if o["name"] else ""
        if i_name and i_name == o_name:
            return 0.5, "same_genus_different_species"

    return 0.0, "discordant_targets"


def compute_taxonomic_confidence(illumina_report, ont_report, selected_taxid):
    selected_taxid = normalize_taxid(selected_taxid)

    if not selected_taxid:
        return {
            "runner_up_taxid": "",
            "runner_up_name": "NA",
            "runner_up_pct": 0.0,
            "delta_pct": 0.0,
            "platform_agreement_score": 0.0,
            "platform_agreement_reason": "no_selected_target",
            "decision_confidence": "low",
            "awareness_state": "fail_no_target"
        }

    top1_i, top2_i = top_two_meaningful_candidates(illumina_report)
    top1_o, top2_o = top_two_meaningful_candidates(ont_report)

    illumina_present = os.path.exists(illumina_report) and os.path.getsize(illumina_report) > 0
    ont_present = os.path.exists(ont_report) and os.path.getsize(ont_report) > 0

    selected_pct_i = parse_kraken_report_target_fraction(illumina_report, selected_taxid) if illumina_present else None
    selected_pct_o = parse_kraken_report_target_fraction(ont_report, selected_taxid) if ont_present else None

    selected_pct, selected_pct_mode = combine_platform_values(
        selected_pct_i,
        selected_pct_o,
        use_both=KRAKEN_AUTO_USE_BOTH_PLATFORMS,
    )

    runner_candidates = []
    for x in [top1_i, top2_i, top1_o, top2_o]:
        if x and normalize_taxid(x["taxid"]) != selected_taxid:
            runner_candidates.append(x)

    runner_candidates = sorted(runner_candidates, key=lambda x: x["pct"], reverse=True)
    runner = runner_candidates[0] if runner_candidates else None

    runner_pct = runner["pct"] if runner else 0.0
    delta_pct = max(selected_pct - runner_pct, 0.0)

    agreement_score, agreement_reason = compute_platform_agreement(illumina_report, ont_report)

    if (
        selected_pct >= KRAKEN_AWARE_HIGH_TARGET_MIN_PCT
        and delta_pct >= KRAKEN_AWARE_DELTA_MIN_PCT
        and agreement_score >= 0.5
    ):
        confidence = "high"
        state = "clean_stable"
    elif (
        selected_pct >= KRAKEN_AWARE_MODERATE_TARGET_MIN_PCT
        and delta_pct >= KRAKEN_AWARE_MODERATE_DELTA_MIN_PCT
    ):
        confidence = "moderate"
        state = "clean_filterable"
    elif selected_pct > 0:
        confidence = "low"
        state = "taxonomically_ambiguous"
    else:
        confidence = "low"
        state = "fail_no_target"

    return {
        "runner_up_taxid": normalize_taxid(runner["taxid"]) if runner else "",
        "runner_up_name": runner["name"] if runner else "NA",
        "runner_up_pct": runner_pct,
        "delta_pct": delta_pct,
        "platform_agreement_score": agreement_score,
        "platform_agreement_reason": agreement_reason,
        "decision_confidence": confidence,
        "awareness_state": state
    }


def top_two_meaningful_candidates(report_path):
    rows = [r for r in parse_kraken_report_rows(report_path) if is_meaningful_kraken_candidate(r)]
    rows = sorted(rows, key=lambda x: x["pct"], reverse=True)
    top1 = rows[0] if len(rows) >= 1 else None
    top2 = rows[1] if len(rows) >= 2 else None
    return top1, top2


def compute_platform_agreement(illumina_report, ont_report):
    i = infer_target_from_kraken_report(illumina_report)
    o = infer_target_from_kraken_report(ont_report)

    if not i["taxid"] or not o["taxid"]:
        return 0.0, "no_inference"

    if i["taxid"] == o["taxid"]:
        return 1.0, "exact_taxid_match"

    if i["rank"] == "S" and o["rank"] == "S":
        i_name = i["name"].split()[0] if i["name"] else ""
        o_name = o["name"].split()[0] if o["name"] else ""
        if i_name and i_name == o_name:
            return 0.5, "same_genus_different_species"

    return 0.0, "discordant_targets"


def compute_taxonomic_confidence(illumina_report, ont_report, selected_taxid):
    selected_taxid = normalize_taxid(selected_taxid)
    if not selected_taxid:
        return {
            "runner_up_taxid": "",
            "runner_up_name": "NA",
            "runner_up_pct": 0.0,
            "delta_pct": 0.0,
            "platform_agreement_score": 0.0,
            "platform_agreement_reason": "no_selected_target",
            "decision_confidence": "low",
            "awareness_state": "fail_no_target"
        }

    top1_i, top2_i = top_two_meaningful_candidates(illumina_report)
    top1_o, top2_o = top_two_meaningful_candidates(ont_report)

    selected_pct_i = parse_kraken_report_target_fraction(illumina_report, selected_taxid)
    selected_pct_o = parse_kraken_report_target_fraction(ont_report, selected_taxid)
    selected_pct = min(selected_pct_i, selected_pct_o) if KRAKEN_AUTO_USE_BOTH_PLATFORMS else max(selected_pct_i, selected_pct_o)

    runner_candidates = []
    for x in [top1_i, top2_i, top1_o, top2_o]:
        if x and normalize_taxid(x["taxid"]) != selected_taxid:
            runner_candidates.append(x)

    runner_candidates = sorted(runner_candidates, key=lambda x: x["pct"], reverse=True)
    runner = runner_candidates[0] if runner_candidates else None
    runner_pct = runner["pct"] if runner else 0.0
    delta_pct = max(selected_pct - runner_pct, 0.0)

    agreement_score, agreement_reason = compute_platform_agreement(illumina_report, ont_report)

    if selected_pct >= 95.0 and delta_pct >= 20.0 and agreement_score >= 0.5:
        confidence = "high"
        state = "clean_stable"
    elif selected_pct >= 80.0 and delta_pct >= 10.0:
        confidence = "moderate"
        state = "clean_filterable"
    elif selected_pct > 0:
        confidence = "low"
        state = "taxonomically_ambiguous"
    else:
        confidence = "low"
        state = "fail_no_target"

    return {
        "runner_up_taxid": normalize_taxid(runner["taxid"]) if runner else "",
        "runner_up_name": runner["name"] if runner else "NA",
        "runner_up_pct": runner_pct,
        "delta_pct": delta_pct,
        "platform_agreement_score": agreement_score,
        "platform_agreement_reason": agreement_reason,
        "decision_confidence": confidence,
        "awareness_state": state
    }

