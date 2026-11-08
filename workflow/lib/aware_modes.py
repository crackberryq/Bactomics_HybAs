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
