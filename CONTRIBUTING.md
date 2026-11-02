# Contributing

## Layout

```
main.nf                       entry point (help, validation, input discovery)
workflows/annovar_nf.nf       the workflow: channels and process wiring
modules/local/<process>/      one process per module, each with a stub block
subworkflows/local/utils_*    helper functions (params, inputs, protocols)
bin/                          Python helpers (stdlib only), on PATH inside tasks
conf/                         resources, publishing, protocols, profiles
tests/                        mock ANNOVAR, test VCFs, Python unit tests
```

## Before opening a pull request

```bash
ruff check && pytest                         # Python helpers
nextflow lint .                              # Nextflow strict syntax
nextflow run main.nf -profile test,conda     # end-to-end with the mock ANNOVAR
nextflow run main.nf -profile test -stub     # wiring only, seconds
```

Guidelines:

- Comments and documentation in English.
- A new process gets its own module directory, a `stub:` block and a `versions.yml` output.
- Python helpers stay dependency-free, stream their input and get unit tests in `tests/python/`.
- A new input edge case gets a small file in `tests/data/` and a row in `tests/data/samplesheet.csv`.
- Update `CHANGELOG.md` and, for user-visible changes, `README.md` / `docs/`.
- Never commit ANNOVAR code or databases (licensed).

On Windows, keep LF line endings (enforced by `.gitattributes`) and mark scripts executable in git:

```bash
git update-index --chmod=+x bin/*.py tests/mock_annovar/table_annovar.pl benchmarks/benchmark.sh examples/*.sh
```

## Releasing

Versions follow [Semantic Versioning](https://semver.org/): new features raise the minor
version, fixes the patch version; from 1.0.0 on, a change to parameters, output columns or the
JSON schema that breaks existing use raises the major version.

1. Set the version in `nextflow.config` (`manifest.version` and the `container` tag),
   `nextflow_schema.json` (`container` default), `CITATION.cff` (`version`, `date-released`) and
   a new `## [x.y.z] - date` section on top of `CHANGELOG.md`; `pytest` checks they agree.
2. Merge to `main` with CI green.
3. Create the GitHub release with tag `vx.y.z` and the CHANGELOG section as notes. The tag
   builds and publishes `ghcr.io/assetdaniyarov/annovar-nf:x.y.z` and `:latest`
   (`.github/workflows/docker.yml`); with the Zenodo integration enabled, the release also gets a
   DOI.
