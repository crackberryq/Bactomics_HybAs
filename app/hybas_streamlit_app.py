import os
import shlex
import subprocess
import signal
from pathlib import Path

import pandas as pd
import psutil
import streamlit as st
import streamlit.components.v1 as components
import yaml
import time

# -----------------------------------------------------------------------------
# 0. SHUTDOWN LOGIC
# -----------------------------------------------------------------------------
if st.query_params.get("shutdown") == "true":
    os._exit(0)

# -----------------------------------------------------------------------------
# 1. PAGE CONFIG & STYLING
# -----------------------------------------------------------------------------
st.set_page_config(page_title="HybAs Control Tower", page_icon="🧬", layout="wide")

if "logs" not in st.session_state:
    st.session_state.logs = ""
if "val_logs" not in st.session_state:
    st.session_state.val_logs = ""
if "proc_pid" not in st.session_state:
    st.session_state.proc_pid = None
if "val_proc_pid" not in st.session_state:
    st.session_state.val_proc_pid = None
if "prepared_cmd" not in st.session_state:
    st.session_state.prepared_cmd = ""
if "confirm_step" not in st.session_state:
    st.session_state.confirm_step = 0

CURRENT_DIR = Path(__file__).resolve().parent
if (CURRENT_DIR / "HybAs.smk").exists() or (CURRENT_DIR / "Hybas.smk").exists():
    ROOT = CURRENT_DIR
else:
    ROOT = CURRENT_DIR.parent

VAL_DIR = ROOT / "validation"
LOG_FILE = ROOT / ".hybas_runtime.log"
VAL_LOG_FILE = ROOT / ".validation_runtime.log"
VAL_SNAKEFILE_NAME = "HybAs_validation.smk"
VAL_CONFIG_NAME = "validation_config.yaml"

st.markdown(
    """
    <style>
    .hybas-header {
        background-color: #2e76b7;
        padding: 0.8rem;
        border-radius: 10px;
        color: white;
        text-align: center;
        margin-bottom: 0.8rem;
        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
    }
    .provenance-badge {
        background-color: #e0f0ff;
        color: #004a99;
        padding: 2px 10px;
        border-radius: 15px;
        font-size: 0.8rem;
        font-weight: bold;
        margin-left: 10px;
        border: 1px solid #b3d7ff;
    }
    .soft-card {
        background: #0f172a;
        border: 1px solid #334155;
        border-radius: 8px;
        padding: 1rem;
        margin-bottom: 0.8rem;
        color: #e2e8f0;
    }
    .command-box {
        background: #1e1e1e;
        color: #76e08d;
        border-radius: 8px;
        padding: 8px;
        font-family: 'Courier New', Courier, monospace;
        font-size: 0.85rem;
        border: 1px solid #333;
        margin-bottom: 10px;
        overflow-x: auto;
    }
    .stButton > button { width: 100%; height: 2.2rem; font-weight: bold; }
    </style>
    <div class="hybas-header">
        <h1 style='margin:0; font-size: 1.5rem;'>🧬 HybAs: Hybrid Assembly & Analysis</h1>
        <p style='margin:5px 0 0 0; font-size: 0.9rem;'>Version: 8.7 | <span class="provenance-badge">Provenance: Version: 1.0 (Active)</span></p>
    </div>
    """,
    unsafe_allow_html=True,
)

# -----------------------------------------------------------------------------
# 2. HELPERS
# -----------------------------------------------------------------------------
def load_yaml(path: Path):
    if not path.exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}
    except Exception:
        return {}


def kill_process_tree(pid):
    try:
        parent = psutil.Process(pid)
        for child in parent.children(recursive=True):
            child.send_signal(signal.SIGKILL)
        parent.send_signal(signal.SIGKILL)
        return True
    except Exception:
        return False


def resolve_threads(profile: str, cpu_count: int) -> int:
    profile_map = {"conservative": 0.25, "balanced": 0.50, "aggressive": 0.75, "max": 0.90}
    frac = profile_map.get(profile, 0.50)
    return max(1, min(cpu_count, int(round(cpu_count * frac))))


def quoted(pathlike) -> str:
    return shlex.quote(str(pathlike))


def cmd_path(pathlike, cwd: Path) -> str:
    return shlex.quote(os.path.relpath(str(Path(pathlike).resolve()), start=str(Path(cwd).resolve())))


def is_pipeline_running():
    for pid_key in ["proc_pid", "val_proc_pid"]:
        pid = st.session_state.get(pid_key)
        if pid:
            try:
                p = psutil.Process(pid)
                if p.is_running() and p.status() != psutil.STATUS_ZOMBIE:
                    return True
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
    return False


# -----------------------------------------------------------------------------
# 3. SETTINGS & PATHS
# -----------------------------------------------------------------------------
DEFAULTS = {
    "project_root": ROOT,
    "main_snakefile": ROOT / "HybAs.smk" if (ROOT / "HybAs.smk").exists() else ROOT / "Hybas.smk",
    "main_config": ROOT / "config.yaml",
    "samples_tsv": ROOT / "samples.tsv",
    "ci_config": ROOT / "ci" / "config.ci.yaml",
    "ci_samples": ROOT / "ci" / "samples.ci.tsv",
}

# -----------------------------------------------------------------------------
# 4. SIDEBAR
# -----------------------------------------------------------------------------
currently_running = is_pipeline_running()
sys_cpu = os.cpu_count() or 1


def sync_threads():
    st.session_state.threads = resolve_threads(st.session_state.res_profile, sys_cpu)


with st.sidebar:
    st.header("1. Control Hub")
    if currently_running:
        st.warning("🔒 System Locked: Workflow Active")

    op_mode = st.selectbox("Operation Mode", ["standard", "admin", "ci"], key="op_mode", disabled=currently_running)

    with st.expander("⚙️ System Scaling", expanded=True):
        st.selectbox(
            "Resource Profile",
            ["conservative", "balanced", "aggressive", "max"],
            index=1,
            key="res_profile",
            on_change=sync_threads,
            disabled=currently_running,
        )
        if "threads" not in st.session_state:
            st.session_state.threads = resolve_threads(st.session_state.get("res_profile", "balanced"), sys_cpu)

        threads_val = st.slider("Max Cores", 1, sys_cpu, key="threads", disabled=(op_mode == "ci" or currently_running))
        use_conda = st.toggle("Use Conda", value=True, disabled=currently_running)
        dry_run = st.toggle("Dry Run (-n)", value=False, disabled=currently_running)
        rerun_incomplete = st.toggle("Rerun Incomplete", value=True, disabled=currently_running)

# -----------------------------------------------------------------------------
# 5. MAIN INTERFACE
# -----------------------------------------------------------------------------
tab_setup, tab_config, tab_run, tab_comp = st.tabs(["📊 Setup", "📝 Config", "🚀 Run & Monitor", "🧮 Computational"])

with tab_run:
    st.markdown('<div class="soft-card">', unsafe_allow_html=True)

    if op_mode == "admin":
        st.subheader("🛠️ Admin Maintenance")
        c1, c2 = st.columns(2)
        with c1:
            if st.button("🔍 Check DB Integrity", disabled=currently_running):
                with st.status("Verifying Databases...", expanded=True) as status:
                    cfg = load_yaml(DEFAULTS["main_config"])
                    db_p = Path(cfg.get("base_dir", ROOT)) / cfg.get("db_root", "db")
                    k_dir = db_p / cfg.get("kraken_db_name", "kraken2_std_db")
                    k_ok = all((k_dir / f).exists() for f in ["hash.k2d", "opts.k2d", "taxo.k2d"])
                    if k_ok:
                        st.success("✅ Kraken2 DB FOUND")
                    else:
                        st.error("❌ Kraken2 DB Files Missing")

                    b_lineage = cfg.get("busco_lineage", "bacteria_odb10")
                    b_dir = db_p / "busco" / b_lineage
                    if b_dir.exists():
                        st.success(f"✅ BUSCO {b_lineage} FOUND")
                    else:
                        st.error(f"❌ BUSCO {b_lineage} MISSING")
                    status.update(label="Integrity Check Complete", state="complete")

            if st.button("📦 Install Conda Envs", disabled=currently_running):
                st.session_state.prepared_cmd = f"snakemake -s {cmd_path(DEFAULTS['main_snakefile'], ROOT)} setup_envs --use-conda -c 4"
                st.session_state.confirm_step = 0

        with c2:
            if st.button("🚀 Update Databases", disabled=currently_running):
                st.session_state.confirm_step = 1

        if st.session_state.confirm_step == 1:
            st.warning("⚠️ WARNING: You are about to download more than 125 GB of data. Proceed?")
            ac1, ac2 = st.columns(2)
            if ac1.button("✅ OK (Step 1/2)"):
                st.session_state.confirm_step = 2
                st.rerun()
            if ac2.button("❌ Cancel"):
                st.session_state.confirm_step = 0
                st.rerun()
        elif st.session_state.confirm_step == 2:
            st.info("ℹ️ Final Confirmation: This will synchronize all IEEE-grade databases.")
            ac3, ac4 = st.columns(2)
            if ac3.button("✅ Final OK (Confirm Download)"):
                st.session_state.prepared_cmd = (
                    f"snakemake -s {cmd_path(DEFAULTS['main_snakefile'], ROOT)} "
                    f"update_databases --config admin_mode=true setup_only=true --use-conda --cores 4"
                )
                st.session_state.confirm_step = 0
                st.rerun()
            if ac4.button("❌ Cancel"):
                st.session_state.confirm_step = 0
                st.rerun()

        cmd = st.session_state.prepared_cmd if st.session_state.prepared_cmd else "Select admin task..."

    else:
        samples_path = DEFAULTS["ci_samples"] if op_mode == "ci" else DEFAULTS["samples_tsv"]
        iso_df = pd.read_csv(samples_path, sep="\t") if samples_path.exists() else pd.DataFrame()
        if not iso_df.empty:
            st.markdown(f"##### Sample Scope: {len(iso_df)} isolate(s) defined in {samples_path.name}")
            st.caption("Run & Monitor will execute using the isolates currently defined in the active sample sheet.")
        else:
            st.warning(f"No sample sheet rows found in {samples_path.name}. Add isolates in Setup before launching runs.")

        c_bits = []
        config_str = f"--config {' '.join(c_bits)}" if c_bits else ""

        if op_mode == "ci":
            cmd = (
                f"snakemake -s {cmd_path(DEFAULTS['main_snakefile'], ROOT)} "
                f"--configfile {quoted(DEFAULTS['ci_config'])} --use-conda --cores 4 "
                f"{'-n' if dry_run else ''} {'--rerun-incomplete' if rerun_incomplete else ''} {config_str}"
            ).strip()
        else:
            cmd = (
                f"snakemake -s {quoted(DEFAULTS['main_snakefile'])} -c {threads_val} "
                f"{'--use-conda' if use_conda else ''} {'-n' if dry_run else ''} "
                f"{'--rerun-incomplete' if rerun_incomplete else ''} {config_str}"
            ).strip()

    st.markdown(f'<div class="command-box">CONDARC=$(pwd)/.condarc {cmd}</div>', unsafe_allow_html=True)

    ctrl_c1, ctrl_c2, ctrl_c3 = st.columns(3)
    b_start = ctrl_c1.button("🚀 START WORKFLOW", type="primary", disabled=currently_running or not cmd or "Select" in cmd)
    b_stop = ctrl_c2.button("🛑 STOP WORKFLOW", disabled=not currently_running)
    if ctrl_c3.button("🔓 Unlock Directory", disabled=currently_running):
        subprocess.run(
            f"CONDARC={quoted(ROOT / '.condarc')} snakemake -s {cmd_path(DEFAULTS['main_snakefile'], ROOT)} --unlock",
            shell=True,
            cwd=str(ROOT),
            executable="/bin/bash",
        )
        st.success("Unlocked.")

    st.markdown("##### Runtime Progress & Output")
    log_view = st.empty()

    if b_start:
        st.session_state.logs = ""
        proc_env = os.environ.copy()
        proc_env["CONDARC"] = str(ROOT / ".condarc")
        with open(LOG_FILE, "w") as f:
            p = subprocess.Popen(
                cmd,
                cwd=str(ROOT),
                shell=True,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True,
                executable="/bin/bash",
                env=proc_env,
                preexec_fn=os.setsid,
            )
            st.session_state.proc_pid = p.pid
        st.rerun()

    if b_stop and st.session_state.proc_pid:
        kill_process_tree(st.session_state.proc_pid)
        st.session_state.proc_pid = None
        st.warning("Workflow Terminated by User.")
        st.rerun()

    if LOG_FILE.exists():
        with open(LOG_FILE, "r") as f:
            lines = f.readlines()
            st.session_state.logs = "".join(lines[-100:])

    log_view.code(st.session_state.logs if st.session_state.logs else "Dashboard Standby.", language="bash")
    st.markdown('</div>', unsafe_allow_html=True)

    if currently_running:
        time.sleep(1.5)
        st.rerun()

with tab_setup:
    path = DEFAULTS["ci_samples"] if op_mode == "ci" else DEFAULTS["samples_tsv"]
    st.subheader("📊 Active Sample Sheet")

    if path.exists():
        setup_df = pd.read_csv(path, sep="\t").fillna("")
    else:
        setup_df = pd.DataFrame(columns=["isolate", "target_taxid", "busco_lineage", "medaka_model"])

    edited_df = st.data_editor(
        setup_df,
        width="stretch",
        disabled=currently_running,
        num_rows="dynamic",
        key="setup_samples_editor",
    )

    c_setup1, c_setup2 = st.columns(2)
    if c_setup1.button("💾 Save Sample Sheet", disabled=currently_running):
        edited_df = edited_df.fillna("")
        if "isolate" not in edited_df.columns:
            st.error("The sample sheet must include an 'isolate' column.")
        else:
            edited_df["isolate"] = edited_df["isolate"].astype(str).str.strip()
            edited_df = edited_df[edited_df["isolate"] != ""].copy()
            for col in ["target_taxid", "busco_lineage", "medaka_model"]:
                if col not in edited_df.columns:
                    edited_df[col] = ""
                edited_df[col] = edited_df[col].astype(str).replace("nan", "").str.strip()
            path.parent.mkdir(parents=True, exist_ok=True)
            edited_df.to_csv(path, sep="\t", index=False)
            st.success(f"Saved {len(edited_df)} row(s) to {path.name}")

    if c_setup2.button("🔄 Reload Sample Sheet", disabled=currently_running):
        st.rerun()

    st.caption("You can edit cells directly, add rows, and save the active sample sheet used by the selected operation mode.")

with tab_config:
    cfg_run_path = DEFAULTS["main_config"]
    cfg_val_path = VAL_DIR / VAL_CONFIG_NAME
    cfg_run = load_yaml(cfg_run_path)
    cfg_val = load_yaml(cfg_val_path)

    cfg_tab_run, cfg_tab_val = st.tabs(["⚙️ Run Configuration", "🧮 Validation Configuration"])

    with cfg_tab_run:
        st.subheader("⚙️ Main HybAs Run Configuration")
        st.caption("Edit config.yaml through structured controls, including Kraken-aware settings.")

        medaka_models = [
            "",
            "r1041_e82_400bps_sup_g615",
            "r1041_e82_400bps_hac_g615",
            "r941_min_sup_g507",
            "r941_min_hac_g507",
        ]

        existing_busco_downloads = cfg_run.get(
            "busco_download_lineages",
            [cfg_run.get("busco_lineage", "bacteria_odb10")],
        )
        if not isinstance(existing_busco_downloads, list):
            existing_busco_downloads = [str(existing_busco_downloads)]

        rc1, rc2 = st.columns(2)

        with rc1:
            run_base_dir = st.text_input(
                "Base directory",
                value=str(cfg_run.get("base_dir", ".")),
                disabled=currently_running,
                key="cfg_run_base_dir",
            )
            run_sample_sheet = st.text_input(
                "Sample sheet",
                value=str(cfg_run.get("sample_sheet", "samples.tsv")),
                disabled=currently_running,
                key="cfg_run_sample_sheet",
            )
            run_threads = st.number_input(
                "Default threads",
                min_value=1,
                value=int(cfg_run.get("threads", 8)),
                step=1,
                disabled=currently_running,
                key="cfg_run_threads",
            )
            run_keep_percent = st.number_input(
                "Filtlong keep_percent",
                min_value=1,
                max_value=100,
                value=int(cfg_run.get("keep_percent", 95)),
                step=1,
                disabled=currently_running,
                key="cfg_run_keep_percent",
            )
            run_racon_rounds = st.number_input(
                "Racon rounds",
                min_value=0,
                value=int(cfg_run.get("racon_rounds", 2)),
                step=1,
                disabled=currently_running,
                key="cfg_run_racon_rounds",
            )

            current_polish_mode = str(cfg_run.get("polish_mode", "tripolish")).strip().lower()
            if current_polish_mode not in ["tripolish", "medaka", "racon", "none"]:
                current_polish_mode = "tripolish"

            run_polish_mode = st.selectbox(
                "Default polish mode",
                ["tripolish", "medaka", "racon", "none"],
                index=["tripolish", "medaka", "racon", "none"].index(current_polish_mode),
                disabled=currently_running,
                key="cfg_run_polish_mode",
            )

            run_workflow_version = st.text_input(
                "Workflow version",
                value=str(cfg_run.get("workflow_version", "HybAs-8.8-ieee-batch")),
                disabled=currently_running,
                key="cfg_run_workflow_version",
            )

        with rc2:
            run_busco_lineage = st.text_input(
                "Default BUSCO lineage",
                value=str(cfg_run.get("busco_lineage", "bacteria_odb10")),
                disabled=currently_running,
                key="cfg_run_busco_lineage",
            )

            current_medaka_model = str(cfg_run.get("medaka_model", ""))
            medaka_options = medaka_models[:] if current_medaka_model in medaka_models else medaka_models + [current_medaka_model]

            run_medaka_model = st.selectbox(
                "Default Medaka model",
                medaka_options,
                index=medaka_options.index(current_medaka_model),
                disabled=currently_running,
                key="cfg_run_medaka_model",
            )

            run_target_taxid = st.text_input(
                "Default target_taxid",
                value=str(cfg_run.get("target_taxid", "")),
                disabled=currently_running,
                key="cfg_run_target_taxid",
            )
            run_db_root = st.text_input(
                "DB root",
                value=str(cfg_run.get("db_root", "db")),
                disabled=currently_running,
                key="cfg_run_db_root",
            )
            run_kraken_db_name = st.text_input(
                "Kraken DB name",
                value=str(cfg_run.get("kraken_db_name", "kraken2_std_db")),
                disabled=currently_running,
                key="cfg_run_kraken_db_name",
            )
            run_busco_download_lineages = st.text_input(
                "BUSCO download lineages (comma-separated)",
                value=", ".join(str(x) for x in existing_busco_downloads),
                disabled=currently_running,
                key="cfg_run_busco_download_lineages",
            )

        st.markdown("##### Module Toggles")
        rm1, rm2, rm3, rm4, rm5 = st.columns(5)
        run_run_kraken = rm1.toggle("run_kraken", value=bool(cfg_run.get("run_kraken", True)), disabled=currently_running, key="cfg_run_run_kraken")
        run_run_busco = rm2.toggle("run_busco", value=bool(cfg_run.get("run_busco", True)), disabled=currently_running, key="cfg_run_run_busco")
        run_run_prokka = rm3.toggle("run_prokka", value=bool(cfg_run.get("run_prokka", True)), disabled=currently_running, key="cfg_run_run_prokka")
        run_run_multiqc = rm4.toggle("run_multiqc", value=bool(cfg_run.get("run_multiqc", True)), disabled=currently_running, key="cfg_run_run_multiqc")
        run_run_polypolish = rm5.toggle("run_polypolish", value=bool(cfg_run.get("run_polypolish", True)), disabled=currently_running, key="cfg_run_run_polypolish")

        st.markdown("##### Workflow Gates")
        rg1, rg2, rg3, rg4, rg5 = st.columns(5)
        run_setup_only = rg1.toggle("setup_only", value=bool(cfg_run.get("setup_only", False)), disabled=currently_running, key="cfg_run_setup_only")
        run_admin_mode = rg2.toggle("admin_mode", value=bool(cfg_run.get("admin_mode", False)), disabled=currently_running, key="cfg_run_admin_mode")
        run_ci_mode = rg3.toggle("ci_mode", value=bool(cfg_run.get("ci_mode", False)), disabled=currently_running, key="cfg_run_ci_mode")
        run_allow_ont_only = rg4.toggle("allow_ont_only", value=bool(cfg_run.get("allow_ont_only", True)), disabled=currently_running, key="cfg_run_allow_ont_only")
        run_allow_illumina_only = rg5.toggle("allow_illumina_only", value=bool(cfg_run.get("allow_illumina_only", True)), disabled=currently_running, key="cfg_run_allow_illumina_only")
        run_require_taxid = st.toggle(
            "require_target_taxid_for_kraken",
            value=bool(cfg_run.get("require_target_taxid_for_kraken", False)),
            disabled=currently_running,
            key="cfg_run_require_taxid",
        )

        st.markdown("##### Kraken-aware Decision Settings")
        rk1, rk2 = st.columns(2)

        with rk1:
            current_kmode = str(cfg_run.get("kraken_decision_mode", "aware")).strip().lower()
            if current_kmode not in ["off", "manual", "auto", "aware"]:
                current_kmode = "aware"

            run_kraken_decision_mode = st.selectbox(
                "kraken_decision_mode",
                ["off", "manual", "auto", "aware"],
                index=["off", "manual", "auto", "aware"].index(current_kmode),
                disabled=currently_running,
                key="cfg_run_kraken_decision_mode",
            )
            run_kraken_auto_passthrough_min_percent = st.number_input(
                "kraken_auto_passthrough_min_percent",
                min_value=0.0,
                max_value=100.0,
                value=float(cfg_run.get("kraken_auto_passthrough_min_percent", 95.0)),
                step=0.5,
                disabled=currently_running,
                key="cfg_run_kraken_auto_passthrough_min_percent",
            )
            run_kraken_auto_prefer_species_min_percent = st.number_input(
                "kraken_auto_prefer_species_min_percent",
                min_value=0.0,
                max_value=100.0,
                value=float(cfg_run.get("kraken_auto_prefer_species_min_percent", 90.0)),
                step=0.5,
                disabled=currently_running,
                key="cfg_run_kraken_auto_prefer_species_min_percent",
            )
            run_kraken_aware_delta_min_percent = st.number_input(
                "kraken_aware_delta_min_percent",
                min_value=0.0,
                max_value=100.0,
                value=float(cfg_run.get("kraken_aware_delta_min_percent", 20.0)),
                step=0.5,
                disabled=currently_running,
                key="cfg_run_kraken_aware_delta_min_percent",
            )
            run_kraken_aware_moderate_delta_min_percent = st.number_input(
                "kraken_aware_moderate_delta_min_percent",
                min_value=0.0,
                max_value=100.0,
                value=float(cfg_run.get("kraken_aware_moderate_delta_min_percent", 10.0)),
                step=0.5,
                disabled=currently_running,
                key="cfg_run_kraken_aware_moderate_delta_min_percent",
            )

        with rk2:
            run_kraken_auto_use_both_platforms = st.toggle(
                "kraken_auto_use_both_platforms",
                value=bool(cfg_run.get("kraken_auto_use_both_platforms", True)),
                disabled=currently_running,
                key="cfg_run_kraken_auto_use_both_platforms",
            )
            run_kraken_include_children = st.toggle(
                "kraken_include_children",
                value=bool(cfg_run.get("kraken_include_children", True)),
                disabled=currently_running,
                key="cfg_run_kraken_include_children",
            )
            run_kraken_aware_high_target_min_percent = st.number_input(
                "kraken_aware_high_target_min_percent",
                min_value=0.0,
                max_value=100.0,
                value=float(cfg_run.get("kraken_aware_high_target_min_percent", 95.0)),
                step=0.5,
                disabled=currently_running,
                key="cfg_run_kraken_aware_high_target_min_percent",
            )
            run_kraken_aware_moderate_target_min_percent = st.number_input(
                "kraken_aware_moderate_target_min_percent",
                min_value=0.0,
                max_value=100.0,
                value=float(cfg_run.get("kraken_aware_moderate_target_min_percent", 80.0)),
                step=0.5,
                disabled=currently_running,
                key="cfg_run_kraken_aware_moderate_target_min_percent",
            )

        if st.button("💾 Save run config.yaml", disabled=currently_running, key="save_run_config"):
            new_run_cfg = dict(cfg_run)
            new_run_cfg.update({
                "base_dir": run_base_dir,
                "sample_sheet": run_sample_sheet,
                "threads": int(run_threads),
                "keep_percent": int(run_keep_percent),
                "racon_rounds": int(run_racon_rounds),
                "busco_lineage": run_busco_lineage,
                "medaka_model": run_medaka_model,
                "target_taxid": run_target_taxid,
                "run_kraken": bool(run_run_kraken),
                "run_busco": bool(run_run_busco),
                "run_prokka": bool(run_run_prokka),
                "run_multiqc": bool(run_run_multiqc),
                "run_polypolish": bool(run_run_polypolish),
                "setup_only": bool(run_setup_only),
                "admin_mode": bool(run_admin_mode),
                "ci_mode": bool(run_ci_mode),
                "allow_ont_only": bool(run_allow_ont_only),
                "allow_illumina_only": bool(run_allow_illumina_only),
                "require_target_taxid_for_kraken": bool(run_require_taxid),
                "polish_mode": run_polish_mode,
                "workflow_version": run_workflow_version,
                "db_root": run_db_root,
                "kraken_db_name": run_kraken_db_name,
                "busco_download_lineages": [x.strip() for x in run_busco_download_lineages.split(",") if x.strip()],
                "kraken_decision_mode": run_kraken_decision_mode,
                "kraken_auto_passthrough_min_percent": float(run_kraken_auto_passthrough_min_percent),
                "kraken_auto_use_both_platforms": bool(run_kraken_auto_use_both_platforms),
                "kraken_include_children": bool(run_kraken_include_children),
                "kraken_auto_prefer_species_min_percent": float(run_kraken_auto_prefer_species_min_percent),
                "kraken_aware_delta_min_percent": float(run_kraken_aware_delta_min_percent),
                "kraken_aware_moderate_delta_min_percent": float(run_kraken_aware_moderate_delta_min_percent),
                "kraken_aware_high_target_min_percent": float(run_kraken_aware_high_target_min_percent),
                "kraken_aware_moderate_target_min_percent": float(run_kraken_aware_moderate_target_min_percent),
            })
            cfg_run_path.write_text(yaml.safe_dump(new_run_cfg, sort_keys=False), encoding="utf-8")
            st.success(f"Saved {cfg_run_path.name}")

    with cfg_tab_val:
        st.subheader("🧮 Validation Workflow Configuration")
        st.caption("Edit validation_config.yaml through the same kind of structured controls used in the Computational tab.")

        vv1, vv2 = st.columns(2)

        with vv1:
            val_isolate = st.text_input(
                "Default isolate",
                value=str(cfg_val.get("isolate", "")),
                disabled=currently_running,
                key="cfg_val_isolate",
            )
            val_phase1_lineage = st.text_input(
                "Phase 1 lineage",
                value=str(cfg_val.get("phase1_lineage", "bacteria_odb10")),
                disabled=currently_running,
                key="cfg_val_phase1_lineage",
            )
            val_phase2_window = st.number_input(
                "Phase 2 window (bp)",
                min_value=1,
                value=int(cfg_val.get("phase2_window_bp", 10000)),
                step=1000,
                disabled=currently_running,
                key="cfg_val_phase2_window",
            )
            val_phase3_window = st.number_input(
                "Phase 3 window (bp)",
                min_value=1,
                value=int(cfg_val.get("phase3_window_bp", 10000)),
                step=1000,
                disabled=currently_running,
                key="cfg_val_phase3_window",
            )
            val_phase3_top_n = st.number_input(
                "Phase 3 top N",
                min_value=1,
                value=int(cfg_val.get("phase3_top_n", 20)),
                step=1,
                disabled=currently_running,
                key="cfg_val_phase3_top_n",
            )

        with vv2:
            val_run_phase1 = st.toggle("run_phase1", value=bool(cfg_val.get("run_phase1", True)), disabled=currently_running, key="cfg_val_run_phase1")
            val_run_phase2 = st.toggle("run_phase2", value=bool(cfg_val.get("run_phase2", True)), disabled=currently_running, key="cfg_val_run_phase2")
            val_run_phase3 = st.toggle("run_phase3", value=bool(cfg_val.get("run_phase3", True)), disabled=currently_running, key="cfg_val_run_phase3")
            val_run_plots = st.toggle("run_plots", value=bool(cfg_val.get("run_plots", True)), disabled=currently_running, key="cfg_val_run_plots")
            val_enable_multiqc = st.toggle("enable_multiqc", value=bool(cfg_val.get("enable_multiqc", True)), disabled=currently_running, key="cfg_val_enable_multiqc")

        if st.button("💾 Save validation_config.yaml", disabled=currently_running, key="save_val_config"):
            new_val_cfg = {
                "isolate": val_isolate,
                "run_phase1": bool(val_run_phase1),
                "run_phase2": bool(val_run_phase2),
                "run_phase3": bool(val_run_phase3),
                "run_plots": bool(val_run_plots),
                "enable_multiqc": bool(val_enable_multiqc),
                "phase1_lineage": val_phase1_lineage,
                "phase2_window_bp": int(val_phase2_window),
                "phase3_window_bp": int(val_phase3_window),
                "phase3_top_n": int(val_phase3_top_n),
            }
            cfg_val_path.write_text(yaml.safe_dump(new_val_cfg, sort_keys=False), encoding="utf-8")
            st.success(f"Saved {cfg_val_path.name}")

with tab_comp:
    st.markdown('<div class="soft-card">', unsafe_allow_html=True)
    st.subheader("🧮 Validation Suite Hub")

    val_samples_path = DEFAULTS["ci_samples"] if op_mode == "ci" else DEFAULTS["samples_tsv"]
    val_iso_df = pd.read_csv(val_samples_path, sep="\t") if val_samples_path.exists() else pd.DataFrame()
    if not val_iso_df.empty:
        st.caption(f"Computational validation will execute using the isolates currently listed in {val_samples_path.name} and the parameters saved in {VAL_CONFIG_NAME}.")
    else:
        st.warning(f"No sample sheet rows found in {val_samples_path.name}. Add isolates in Setup before launching validation.")

    st.markdown("##### Execution Controls")
    vc1, vc2, vc3 = st.columns(3)
    val_dry = vc1.toggle("Dry Run (-n)", value=False, key="val_dry", disabled=currently_running)
    val_conda = vc2.toggle("Use Conda", value=True, key="val_conda", disabled=currently_running)
    val_rerun = vc3.toggle("Rerun Incomplete", value=True, key="val_rerun", disabled=currently_running)

    v_c_bits = []

    val_flags = f"{'-n' if val_dry else ''} {'--use-conda' if val_conda else ''} {'--rerun-incomplete' if val_rerun else ''}"
    val_cmd = (
        f"snakemake -s {cmd_path(VAL_DIR / VAL_SNAKEFILE_NAME, VAL_DIR)} -c {threads_val} {val_flags} "
        f"--configfile {cmd_path(VAL_DIR / VAL_CONFIG_NAME, VAL_DIR)} --config {' '.join(v_c_bits)}"
    ).strip()

    st.markdown(f'<div class="command-box">CONDARC=$(pwd)/.condarc {val_cmd}</div>', unsafe_allow_html=True)

    v_ctrl1, v_ctrl2, v_ctrl3 = st.columns(3)
    val_start = v_ctrl1.button("▶ START VALIDATION", type="primary", disabled=currently_running, key="b_val_start")
    val_stop = v_ctrl2.button("⏹ STOP VALIDATION", disabled=not st.session_state.val_proc_pid, key="b_val_stop")
    if v_ctrl3.button("🔓 Unlock Validation", disabled=currently_running, key="b_val_unlock"):
        subprocess.run(
            f"CONDARC={quoted(ROOT / '.condarc')} snakemake -s {cmd_path(VAL_DIR / VAL_SNAKEFILE_NAME, VAL_DIR)} --unlock",
            shell=True,
            cwd=str(VAL_DIR),
            executable="/bin/bash",
        )
        st.success("Validation Lock Released.")

    st.markdown("##### Validation Progress & Telemetry")
    v_log_view = st.empty()

    if val_start:
        st.session_state.val_logs = ""
        v_proc_env = os.environ.copy()
        v_proc_env["CONDARC"] = str(ROOT / ".condarc")
        with open(VAL_LOG_FILE, "w") as f:
            p = subprocess.Popen(
                val_cmd,
                cwd=str(VAL_DIR),
                shell=True,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True,
                executable="/bin/bash",
                env=v_proc_env,
                preexec_fn=os.setsid,
            )
            st.session_state.val_proc_pid = p.pid
        st.rerun()

    if val_stop and st.session_state.val_proc_pid:
        kill_process_tree(st.session_state.val_proc_pid)
        st.session_state.val_proc_pid = None
        st.warning("Validation Stopped.")
        st.rerun()

    if VAL_LOG_FILE.exists():
        with open(VAL_LOG_FILE, "r") as f:
            v_lines = f.readlines()
            st.session_state.val_logs = "".join(v_lines[-100:])

    v_log_view.code(st.session_state.val_logs if st.session_state.val_logs else "Validation Dashboard Standby.", language="bash")
    st.markdown('</div>', unsafe_allow_html=True)

components.html(
    """
    <script>
    window.addEventListener('unload', function() {
        navigator.sendBeacon('/?shutdown=true');
    });
    </script>
    """,
    height=0,
)