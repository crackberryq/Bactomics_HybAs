# Medaka model selection

A new HybAs configuration template sets
`medaka_model_selection_policy: documented_default`. Existing configurations
that omit the setting retain `strict`; saved source copies, settings, receipts
and outputs are not migrated in place. Bactologix must approve the exact source
with model-review capability version 2 before offering this policy for new drafts.

Selection happens before costly processing, using every raw ONT read, and is
checked again against every filtered read before Medaka execution:

1. A researcher-supplied `medaka_model` has priority. An unmatched choice still
   requires the recorded acknowledgement, reason and exact model binding.
2. Otherwise, one basecaller identity recorded on every read is matched exactly
   against the installed Medaka mapping.
3. With `documented_default`, wholly absent metadata or one unsupported
   basecaller identity uses the pinned Medaka 2.2.0 consensus default
   `r1041_e82_400bps_sup_v5.2.0`. The installed version and default identity must
   agree with the reviewed source; its `medaka/options.py` SHA256 must also match. This is labelled **fallback**, not a match.

`strict` disables step 3. Mixed identities, partial annotations, an unavailable
installed mapping, and malformed/empty FASTQ inputs require attention. No common
chemistry is inferred from mixed model names. Researchers may provide an explicit
reviewed model; the workflow does not make that research decision for them.

A failure resolving or downloading a chosen model remains a failure. A damaged
model, unavailable device, failed inference or other runtime error never causes
another model to be selected or Medaka to be omitted. Native compatibility checks
remain active. Earlier attempts and generated intermediate outputs are retained.

## Upstream basis and limits

The reviewed upstream source is Medaka 2.2.0, commit
`019a2ab597218a1c239fb7a5315f271a4577ed6b`:

- [options.py, default consensus model](https://github.com/nanoporetech/medaka/blob/019a2ab597218a1c239fb7a5315f271a4577ed6b/medaka/options.py#L9-L22)
  defines the default model above.
- [medaka_consensus, unsuccessful automatic selection](https://github.com/nanoporetech/medaka/blob/019a2ab597218a1c239fb7a5315f271a4577ed6b/scripts/medaka_consensus#L115-L129)
  warns and retains that default when automatic detection returns no model in
  its ordinary consensus mode. Its bacterial mode has a separate error branch;
  this HybAs policy does not enable or emulate `--bacteria`.
- [README, model selection](https://github.com/nanoporetech/medaka/blob/019a2ab597218a1c239fb7a5315f271a4577ed6b/README.md#models)
  explains matching to basecaller metadata and supplying an explicit model when
  the basecaller identity is known but missing from the input.

HybAs implements a narrower policy than catching arbitrary upstream resolver
failure: it distinguishes eligible metadata conditions from installation,
download, input and runtime errors before resolving the selected model. The
source-defined default is not evidence of suitability for an unknown basecaller,
mixed chemistry, this assembly route or a particular research endpoint. No model
accuracy or historical benchmark equivalence is established by this change.

## Recorded evidence

Receipts preserve the requested model (empty when no researcher chose one),
selection policy, original metadata assessment, selection reason, selected model
ID, resolved file, model SHA256, installed Medaka version, resolution command and
upstream default source commit. A fallback records `selection: fallback`,
`metadata_match: false` and `scientific_suitability: not_assessed`. The execution
receipt binds the decision, policy, software and model hash to preflight.

Native summaries, reports and Methods drafts identify fallback use and its
reason. Missing historical evidence is not reconstructed from current defaults.
Software fixtures check precedence, failure handling and receipt binding; they
do not establish model performance or identical results to earlier runs.
