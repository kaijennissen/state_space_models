FROM julia:1.3.1-buster

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    apt-utils \
    gcc \
    g++ \
    openssh-server \
    cmake \
    build-essential \
    gdb \
    gdbserver \
    rsync \
    vim \
    locales

RUN apt-get install -y \
    bzip2 \
    wget \
    gnupg \
    dirmngr \
    apt-transport-https \
    ca-certificates \
    openssh-server \
    tmux \
    fortune \
    cowsay \
    figlet && \
        apt-get clean

#setup ssh
RUN mkdir /var/run/sshd && \
    echo 'root:root_pwd' |chpasswd && \
        sed -ri 's/^#?PermitRootLogin\s+.*/PermitRootLogin yes/' \
            /etc/ssh/sshd_config && \
                sed -ri 's/UsePAM yes/#UsePAM yes/g' /etc/ssh/sshd_config && \
                    mkdir /root/.ssh

#remove leftovers
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Expose 22 for ssh server. 7777 for gdb server.
EXPOSE 22 7777


########################################################
# Add custom packages and development environment here
########################################################
#RUN julia -O3 -e 'using Pkg;Pkg.REPLMode.pkgstr("add CSV   ;precompile");using CSV'
#RUN julia -O3 -e 'using Pkg;Pkg.REPLMode.pkgstr("add https://github.com/JuliaComputing/MKL.jl ;precompile");using CSV'
#RUN julia -O3 -e 'using Pkg;Pkg.REPLMode.pkgstr("add Distributions   ;precompile");using Distributions'
#RUN julia -O3 -e 'using Pkg;Pkg.REPLMode.pkgstr("add Plots   ;precompile");using Plots'

########################################################
COPY startup.sh / 

ENV PATH $PATH:/usr/games

CMD ["/bin/sh", "startup.sh"]

