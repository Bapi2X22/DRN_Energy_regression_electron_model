singularity exec --nv /cvmfs/unpacked.cern.ch/registry.hub.docker.com/fastml/triton-torchgeo:22.07-py3-geometric \
tritonserver --model-repository /tmp/models --http-port 9000 --grpc-port 9001 --metrics-port 9002 --allow-http=1 > triton.log 2>&1 & 
