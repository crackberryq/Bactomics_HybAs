"""Explicit AWARE+ modes; absent/legacy mode preserves older configuration."""
MODES = ('legacy', 'off', 'taxonomy', 'auto', 'split')

def resolve(config):
    cfg = dict(config)
    mode = cfg.get('aware_mode', 'legacy')
    if mode == 'reconstruct':
        raise ValueError('AWARE+ Reconstruct is not implemented.')
    if mode not in MODES:
        raise ValueError('aware_mode must be ' + ', '.join(MODES))
    if mode == 'legacy':
        return cfg
    cfg.update(run_kraken=mode != 'off', kraken_decision_mode='off' if mode == 'off' else 'aware',
               aware_preserve_groups=mode == 'split')
    if mode == 'off':
        cfg.update(kraken_review_before_filtering=False, aware_run_until='assembly', aware_statistics=False)
    if mode == 'split':
        # Researcher selects eligible derived groups before separate assembly drafts.
        cfg['aware_run_until'] = 'aware'
    return cfg


def kraken_skipped(decision, summary):
    """Interpret an explicit off mode in current and historical run records."""
    mode = decision.get('kraken_decision_mode') or summary.get('kraken_decision_mode')
    if mode:
        return mode == 'off'
    return (decision.get('kraken_decision_reason') == 'kraken_mode_off'
            or summary.get('kraken_decision_reason') == 'kraken_mode_off')
