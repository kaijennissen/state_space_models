FROM julia:1.4.0-buster

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
    locales \
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
RUN mkdir /var/run/sshd
RUN echo 'root:root_pwd' |chpasswd
RUN sed -ri 's/^#?PermitRootLogin\s+.*/PermitRootLogin yes/' \
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
RUN julia -O3 -e 'using Pkg; Pkg.add("CSV"); using CSV'
RUN julia -O3 -e 'using Pkg; Pkg.add("Distributions"); using Distributions'
RUN julia -O3 -e 'using Pkg; Pkg.add("Plots"); using Plots'
RUN julia -O3 -e 'using Pkg; Pkg.add("JuliaInterpreter"); using JuliaInterpreter'
RUN julia -O3 -e 'using Pkg; Pkg.add("Infiltrator"); using Infiltrator'
RUN julia -O3 -e 'using Pkg; Pkg.add("BenchmarkTools"); using BenchmarkTools'
RUN julia -O3 -e 'using Pkg; Pkg.add("StaticArrays"); using StaticArrays'
#RUN julia -O3 -e 'using Pkg; Pkg.add("Infiltrator"); using Infiltrator'
#RUN julia -O3 -e 'using Pkg; Pkg.REPLMode.pkgstr("add https://github.com/JuliaComputing/MKL.jl")'

########################################################
COPY startup.sh / 

ENV PATH $PATH:/usr/games

CMD ["/bin/sh", "startup.sh"]

