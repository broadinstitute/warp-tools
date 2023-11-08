# DOCKER PATCH BY MSFT

## Quick reference

Copy and paste either commands to pull this image (from GCR or ACR)

#### `docker pull us.gcr.io/broad-gotc-prod/docker-path-patch:1.0`
#### `docker pull dsppipelinedev.azurecr.io/docker-path-patch:1.0`
`

- __What is this image:__ This image is a python image for running docker path patch script created by Microsoft to help us transition pipelines to Azure
- __What is Docker Path Patch:__ This program will parse a WDL (using miniwdl) to find instances where the docker runtime parameter is hard-coded as a string URL or specified as a variable with a default in the caller, and replace those instances with variables that are passed up to the root workflow. This allows the user to easily change docker locations without having to modify the WDL.


## Usage

### Docker Path Patch (replace placeholder variables accordingly)

```bash
$ docker run --rm -it \
    -v /replace-with-path-to-WDL-dir/:/data \
    us.gcr.io/broad-gotc-prod/docker-path-patch:1.0 bash;

$ python3 patch_docker_runtime.py -i /data/replace-with-wdl-filename -d /data/optional-hints.json -o /data
```

### Example Run
```bash
$ docker run --rm -it -v /Users/palis/Work/Repositories/smartseq2-single-nucleus-multisample:/data \
    us.gcr.io/broad-gotc-prod/docker-path-patch:kp_containerize_patching_script_pd-2400 bash

$ python3 patch_docker_runtime.py -i /data/SmartSeq2SingleSample.wdl -d /data/docker_hints.json -o /data
``````