# Docker files

**_TMvis_** currently supports membrane embeddings of the methods PPM3 and ANVIL. These can be accessed via the two separate docker files in this folder.

## First run
### PPM3 
The docker image `Dockerfile_ppm3` can be build by running ```sudo docker image build . --rm -t ppm3 --file Dockerfile_ppm3```
### ANVIL
To run ANVIL you need the software [naccess](http://www.bioinf.manchester.ac.uk/naccess/). While the Dockerfile will install the software itself, you need a key that must be copied added to `Dockefile_anvil`. In brief, you should send an email to simon.hubbard@manchester.ac.uk with the signed CONFIDENTIALITY AGREEMENT (can be found in [here](http://www.bioinf.manchester.ac.uk/naccess/)) and enter the key into the placeholder in `Dockerfile_ANVIL`. Then you can build the Docker image by running ```sudo docker image build . --rm -t anvil --file Dockerfile_anvil```


## Charliecloud

We will demonstrate the usage of Dockerfile with Charliecloud setup (https://hpc.github.io/charliecloud/index.html). However, from this step the native docker can be used, too. For this, you need to adjust nextflow config file `custom.config` by enabling Docker like `docker.enabled = true` and changing pathes to container to docker register. See Nextflow documentation pages for more https://www.nextflow.io/docs/latest/config.html

```
sudo ch-convert -i docker -o tar anvil:latest anvil.tar.gz
sudo ch-convert -i docker -o tar ppm3:latest ppm3.tar.gz 

ch-tar2dir anvil.tar.gz charliecloud/
ch-tar2dir ppm3.tar.gz charliecloud/
```
