
# Pre-publication checklist

- [ ] Replace placeholder author and repository fields in `CITATION.cff`.
- [ ] Add data archive accessions and DOI to `README.md`.
- [ ] Select and add approved code and metadata licenses.
- [ ] Remove or parameterize cluster-specific absolute paths and account names.
- [ ] Confirm no credentials, tokens, or protected data are present.
- [ ] Run `python -m py_compile scripts/*.py`.
- [ ] Run a shell linter on `scripts/*.sh`.
- [ ] Verify `sha256sum -c checksums/SHA256SUMS`.
- [ ] Test representative workflows against archived inputs.
- [ ] Create annotated tag `v1.0.0` after all blockers are closed.
