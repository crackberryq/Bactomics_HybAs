# Security policy

## Supported versions

Security fixes are prioritized for the current supported HybAs release branch. Historical branches are retained for reproducibility but are not automatically backported.

## Reporting a vulnerability

Please **do not open a public issue** for a vulnerability involving arbitrary command execution, path traversal, unsafe input handling, credential disclosure, or other security-sensitive behavior. Use the repository's GitHub **Security → Report a vulnerability** private-reporting feature if it is enabled. If private reporting is unavailable, contact the repository maintainer through an existing private channel to request a secure disclosure route; do not post exploit details publicly.

Include the affected commit/version, operating system, a minimal non-sensitive reproduction, expected versus observed behavior, and potential impact. Do not include sequencing data, sample identifiers, private paths, passwords, or API tokens.

The project processes locally supplied sequencing files and invokes external bioinformatics tools. Only run trusted configurations, software environments, and databases. Security reporting does not constitute clinical validation or certification.
