FROM mambaorg/micromamba:1.5.0

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:/usr/bin:$PATH \
    MAMBA_EXTRACT_THREADS=4

SHELL ["bash", "-lc"]
USER root

# Update micromamba and all packages to latest versions
RUN micromamba update --all -y

# Install system dependencies in a single layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-pip \
    python3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    ca-certificates \
    curl \
    git \
    unzip \
    zsh \
    bash \
    libcairo2-dev \
    libbz2-dev \
    liblzma-dev \
    wget \
    software-properties-common \
    dirmngr \
    lsb-release \
    gnupg2 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install core system packages first (conda-forge only)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba install --channel-priority strict -c conda-forge \
    libstdcxx-ng \
    python=3.10 \
    starship \
    cmake \
    make \
    gcc \
    gxx \
    datamash \
    openjdk=17 \
    pip \
    -y && micromamba clean --all --yes

# Install lsd manually (not available in conda-forge)
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
    https://github.com/lsd-rs/lsd/releases/download/v1.1.5/lsd-v1.1.5-x86_64-unknown-linux-gnu.tar.gz && \
    tar -xzf lsd-v1.1.5-x86_64-unknown-linux-gnu.tar.gz && \
    mv lsd-v1.1.5-x86_64-unknown-linux-gnu/lsd /usr/local/bin/ && \
    rm -rf lsd-v1.1.5-*

# Install Ruby and colorls
RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba install --channel-priority strict -c conda-forge \
    ruby \
    -y && micromamba clean --all --yes

# Set up Ruby gem environment and install colorls
RUN export GEM_HOME="/opt/conda/share/rubygems" && \
    export GEM_PATH="/opt/conda/share/rubygems" && \
    export PATH="/opt/conda/share/rubygems/bin:$PATH" && \
    gem install colorls

# Create colorls configuration directory and file
RUN mkdir -p ~/.config/colorls && \
    echo "unrecognized_file: white" > ~/.config/colorls/dark_colors.yaml && \
    echo "recognized_file: white" >> ~/.config/colorls/dark_colors.yaml && \
    echo "executable_file: red" >> ~/.config/colorls/dark_colors.yaml && \
    echo "dir: blue" >> ~/.config/colorls/dark_colors.yaml && \
    echo "user: magenta" >> ~/.config/colorls/dark_colors.yaml && \
    echo "group: cyan" >> ~/.config/colorls/dark_colors.yaml && \
    echo "date: yellow" >> ~/.config/colorls/dark_colors.yaml && \
    echo "time: darkgreen" >> ~/.config/colorls/dark_colors.yaml && \
    echo "file_size: palegreen" >> ~/.config/colorls/dark_colors.yaml && \
    echo "read: darkgreen" >> ~/.config/colorls/dark_colors.yaml && \
    echo "write: yellow" >> ~/.config/colorls/dark_colors.yaml && \
    echo "exec: red" >> ~/.config/colorls/dark_colors.yaml && \
    echo "no_access: gray" >> ~/.config/colorls/dark_colors.yaml && \
    echo "image: magenta" >> ~/.config/colorls/dark_colors.yaml && \
    echo "video: blue" >> ~/.config/colorls/dark_colors.yaml && \
    echo "music: cyan" >> ~/.config/colorls/dark_colors.yaml && \
    echo "log: yellow" >> ~/.config/colorls/dark_colors.yaml

# Install Jupyter ecosystem (conda-forge only)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba install --channel-priority strict -c conda-forge \
    jupyter \
    jupyterlab \
    notebook \
    ipykernel \
    -y && micromamba clean --all --yes

# Install bioinformatics tools (conda-forge + bioconda)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba install --channel-priority strict -c conda-forge -c bioconda \
    bcftools \
    samtools \
    bedtools \
    tabix \
    snakemake \
    mamba \
    -y && micromamba clean --all --yes

# Install R base and Bioconductor packages (conda-forge + bioconda)
RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba install --channel-priority strict -c conda-forge -c bioconda \
    r-base \
    bioconductor-deseq2 \
    bioconductor-tximport \
    bioconductor-annotationdbi \
    bioconductor-biomart \
    bioconductor-sva \
    sra-tools \
    entrez-direct \
    -y && micromamba clean --all --yes

# Set LD_LIBRARY_PATH to use conda's libstdc++
ENV LD_LIBRARY_PATH=/opt/conda/lib

# Install essential Python packages only
RUN pip3 install --no-cache-dir \
    pandas \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    requests \
    beautifulsoup4 \
    xmltodict \
    lxml \
    biopython \
    scikit-learn

# Install ALL R packages in a single layer
RUN R -e "install.packages(c('here', 'data.table', 'metafor', 'tidyverse', 'ggplot2', 'qqman', 'qqplotr', 'reticulate', 'broom', 'readxl', 'writexl', 'knitr', 'rmarkdown'), repos='https://cloud.r-project.org/', Ncpus=4)"

# Install and configure shell environment in a single layer
RUN git clone --depth=1 https://github.com/ohmyzsh/ohmyzsh.git ~/.oh-my-zsh && \
    git clone --depth=1 https://github.com/romkatv/powerlevel10k.git ~/.oh-my-zsh/custom/themes/powerlevel10k && \
    git clone --depth=1 https://github.com/dracula/zsh.git ~/.oh-my-zsh/themes/dracula && \
    git clone --depth=1 https://github.com/zsh-users/zsh-completions.git ~/.oh-my-zsh/custom/plugins/zsh-completions && \
    git clone --depth=1 https://github.com/zsh-users/zsh-autosuggestions.git ~/.oh-my-zsh/custom/plugins/zsh-autosuggestions && \
    git clone --depth=1 https://github.com/zsh-users/zsh-syntax-highlighting.git ~/.oh-my-zsh/custom/plugins/zsh-syntax-highlighting && \
    cp ~/.oh-my-zsh/templates/zshrc.zsh-template ~/.zshrc && \
    sed -i 's/ZSH_THEME=".*"/ZSH_THEME="dracula\/dracula"/' ~/.zshrc && \
    sed -i 's/plugins=(git)/plugins=(git zsh-completions zsh-autosuggestions zsh-syntax-highlighting)/' ~/.zshrc && \
    echo 'export DISABLE_AUTO_UPDATE="true"' >> ~/.zshrc && \
    echo 'export DISABLE_UPDATE_PROMPT="true"' >> ~/.zshrc && \
    echo 'POWERLEVEL9K_DISABLE_CONFIGURATION_WIZARD=true' >> ~/.zshrc && \
    sed -i 's/typeset -g POWERLEVEL9K_TIME_BACKGROUND=.*/typeset -g POWERLEVEL9K_TIME_BACKGROUND=magenta/' ~/.p10k.zsh || echo 'typeset -g POWERLEVEL9K_TIME_BACKGROUND=magenta' >> ~/.p10k.zsh && \
    echo 'if command -v lsd > /dev/null; then' >> ~/.zshrc && \
    echo '  sed -i "/alias ls=/d" ~/.zshrc' >> ~/.zshrc && \
    echo '  sed -i "/LS_COLORS=/d" ~/.zshrc' >> ~/.zshrc && \
    echo '  export LS_COLORS="di=1;34:ln=1;36:so=1;35:pi=1;94:ex=1;31:bd=1;95:cd=1;96:ur=0;32:uw=0;33:ux=0;31:ue=0;32:gr=0;32:gw=0;33:gx=0;31:tr=0;90:tw=0;93:tx=0;92"' >> ~/.zshrc && \
    echo '  alias ls="lsd --color=always --header"' >> ~/.zshrc && \
    echo '  alias ll="colorls --long --almost-all --sort-dirs --git-status"' >> ~/.zshrc && \
    echo '  alias la="lsd -la --color=always --header"' >> ~/.zshrc && \
    echo '  alias lt="lsd --tree --color=always --header"' >> ~/.zshrc && \
    echo 'fi' >> ~/.zshrc && \
    echo '# Match local colorls environment' >> ~/.zshrc && \
    echo 'export LS_COLORS="di=1;34:ln=1;36:so=1;35:pi=1;94:ex=1;31:bd=1;95:cd=1;96:ur=0;32:uw=0;33:ux=0;31:ue=0;32:gr=0;32:gw=0;33:gx=0;31:tr=0;90:tw=0;93:tx=0;92"' >> ~/.zshrc && \
    echo 'export LSCOLORS="Gxfxcxdxbxegedabagacad"' >> ~/.zshrc && \
    echo 'export TERM="xterm-256color"' >> ~/.zshrc && \
    echo 'export COLORTERM="truecolor"' >> ~/.zshrc && \
    echo 'export EZA_COLORS="ur=0:uw=0:ux=0:ue=0:gr=0:gw=0:gx=0:tr=0:tw=0:tx=0:su=0:sf=0:xa=0"' >> ~/.zshrc && \
    echo 'export COLORFGBG="15;0"' >> ~/.zshrc && \
    echo '# Ruby gem environment' >> ~/.zshrc && \
    echo 'export GEM_HOME="/opt/conda/share/rubygems"' >> ~/.zshrc && \
    echo 'export GEM_PATH="/opt/conda/share/rubygems"' >> ~/.zshrc && \
    echo 'export PATH="/opt/conda/share/rubygems/bin:$PATH"' >> ~/.zshrc && \
    echo '# Initialize conda' >> ~/.zshrc && \
    echo 'export PATH="/opt/conda/bin:${PATH}"' >> ~/.zshrc && \
    echo 'eval "$(starship init zsh)"' >> ~/.zshrc && \
    echo 'export PATH="/opt/conda/bin:${PATH}"' >> ~/.bashrc

# Set working directory and create directory structure
WORKDIR /proj
RUN mkdir -p /proj/data/{raw,metadata,references} \
    /proj/scripts/{01_download_data,02_run_rnaseq,03_qc_analysis,04_differential_expression,05_visualization,utils} \
    /proj/logs

# Set environment variables
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# Copy project files
COPY . /proj/

# Final environment setup
ENTRYPOINT ["zsh", "-l", "-c", "source ~/.zshrc && exec zsh"]
CMD []
