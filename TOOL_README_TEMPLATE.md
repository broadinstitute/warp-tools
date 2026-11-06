<!--
WARP-Tools image README template.

Copy this file to <image-dir>/README.md and fill in every section. Delete the
HTML comments. Keep it accurate and runnable: every command should work if
copy-pasted. Required sections: Synopsis, Inputs, Outputs, Usage, Task
descriptions. The others are recommended; drop a section only if it genuinely
does not apply (and say why). See 3rd-party-tools/scvi-scanvi/README.md for a
worked example, and BUILDING.md for build/versioning conventions.
-->

# <Image Name>

## Quick reference

<!-- One copy-pasteable pull command (current tag from docker_versions.tsv) and a
one-line description of what this image is. -->

```bash
docker pull us.gcr.io/broad-gotc-prod/<image>:<image-version>-<tool-version>-<unix-timestamp>
```

- **What is this image:** <one-sentence description of the image and its base OS>
- **What is `<tool>`:** <one-sentence description of the wrapped tool, with a link to its upstream docs>

## Synopsis

<!-- 2-4 sentences: what the tool/scripts in this image do, and what role the
image plays in WARP (which pipeline/step consumes it, and why it exists). -->

## Image contents

<!-- The base image and the key tools/libraries with their pinned versions.
A short bullet list is fine. -->

- Base image: `<base>`
- `<tool>` `<version>`
- `<library>` `<version>`
- Scripts: `<script1>`, `<script2>` (see [Scripts](#scripts))

## Versioning

This image follows the WARP-Tools tag convention
`us.gcr.io/broad-gotc-prod/<image>:<image-version>-<tool-version>-<unix-timestamp>`
(see [BUILDING.md → Versioning strategy](../../BUILDING.md#versioning-strategy)).
Past versions are tracked in [docker_versions.tsv](docker_versions.tsv); the last
line is the version currently used by WARP. Inspect an image with:

```bash
docker inspect us.gcr.io/broad-gotc-prod/<image>:<tag>
```

## Scripts

<!-- OPTIONAL: list each independent script (Python/C++/etc.) baked into the image
and what it does. Omit if the image wraps a single binary with no custom scripts. -->

| Script | Language | Purpose |
| --- | --- | --- |
| `<script.py>` | Python | <what it does> |

## Inputs

<!-- The inputs the tool/scripts consume: command-line arguments, input files and
their formats. If there are multiple scripts/entrypoints, group by script. Be
explicit about file formats (e.g. AnnData .h5ad, BAM, GTF) and required columns/fields. -->

| Input | Type / format | Required | Description |
| --- | --- | --- | --- |
| `<--arg>` | `<type>` | yes/no | <description> |
| `<input file>` | `<format>` | yes/no | <description> |

## Outputs

<!-- The files/artifacts produced, with formats and (if useful) where they are written. -->

| Output | Format | Description |
| --- | --- | --- |
| `<output file>` | `<format>` | <description> |

## Usage

<!-- Runnable examples. Show the docker run invocation and the script/tool command(s).
Include a minimal end-to-end example with representative arguments. -->

Run an interactive shell in the image:

```bash
docker run --rm -it us.gcr.io/broad-gotc-prod/<image>:<tag> bash
```

Run the tool/script:

```bash
<command with representative arguments>
```

## Task descriptions

<!-- How WARP consumes this image: which pipeline(s)/WDL task(s) use it, and how
(e.g. "imports run_multi_model from /usr/local", or "invoked as `<binary> ...`").
Link to the consuming pipeline. This is the bridge back to WARP. -->

- **Pipeline:** [`<Pipeline>`](https://github.com/broadinstitute/warp/tree/develop/pipelines/wdl/<path>)
- **Task(s):** <task name(s)> — <how the image is invoked / what functions are imported>

## Requirements / dependencies

<!-- OPTIONAL: notable runtime requirements (e.g. GPU, network access at runtime,
large memory/disk), and how to install deps if running scripts standalone. -->

## Troubleshooting

To run the image independently of a WDL for testing, start a shell explicitly
(`docker run -it --rm <image-url> bash`). See
[BUILDING.md → Troubleshooting and running standalone](../../BUILDING.md#troubleshooting-and-running-standalone).
