FROM python:3.7-slim-buster

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update --fix-missing && \
    apt-get install --no-install-recommends -y \      
	procps \
	git \
        make \
        gcc \
        libhts-dev \
        libbz2-dev \
        python3-dev \
        python3-pip && \
    apt-get clean && \
    rm -rf /tmp/downloaded_packages/* && \
    rm -rf /var/lib/apt/lists/*

# RUN pip install git+https://gitlab.com/milothepsychic/filter_sam
ADD ./* /opt/filter_sam/
RUN pip install -e /opt/filter_sam

CMD [ "filter_sam" ]