```
function test() {
  console.log("This code will have a copy button to the right of it");
}
```

# CliMush Bioinformatics Pipeline

## Getting Started

### Docker + Singularity

Everything needed to run the bioinformatics pipeline is contained within a Docker image. This includes the CliMush Python package, custom pipeline scripts, and version-controlled open-source software like VSEARCH. The image also sets up the file structure required to running the bioinformatics pipeline.

If you are planning to run the pipeline primarily on a high-performance computing cluster, then consult the Singularity subsection. If you plan to run the pipeline on your local machine, then see the Docker subsection.

#### Singularity

#### Docker

In order to access and run the pipeline in version-controlled, platform-independent manner, you first need to pull the Docker image from Docker Hub and run a container from this image.

The Docker image is located on [Docker Hub](https://hub.docker.com/r/cdelevich/climush-bioinfo). In order to run  a Docker container from this image, copy the code below and run it in Terminal.
```
docker image pull cdelevic/climush-bioinfo:latest
```
## Description

A collection of scripts for processing raw reads from Sanger, Illumina, and PacBio sequencing technologies.

The steps of the pipeline are ordered, but can be adjusted based on user needs.
1. sort files
2. demultiplex
3. pre-filter
3. trim primers
3. dereplicate full-length reads
4. quality filter
5. merge paired-end reads
6. etc.

This README is incomplete:
-[ ] finish writing the README file

> [!CAUTION]
> Not ready for use.

```
git status
```