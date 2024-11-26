# Используем базовый образ Miniconda
FROM continuumio/miniconda3

# Устанавливаем рабочую директорию
WORKDIR /app

# Устанавливаем Mamba
RUN conda install -n base -c conda-forge mamba

# Копируем только файл environment.yml
COPY environment.yml .

# Создаем окружение с помощью Mamba
RUN mamba env create -f environment.yml

# Активируем окружение
ENV PATH /opt/conda/envs/lpce/bin:$PATH
ENV CONDA_DEFAULT_ENV=lpce

# Устанавливаем дополнительные системные зависимости
RUN apt-get update && \
    apt-get install -y --no-install-recommends make tmux && \
    rm -rf /var/lib/apt/lists/*

# Копируем .dockerignore для исключения ненужных файлов
COPY .dockerignore .

# Копируем оставшиеся файлы проекта
COPY . /app

# Устанавливаем SHELL для использования окружения
SHELL ["conda", "run", "-n", "lpce", "/bin/bash", "-c"]

# Устанавливаем команду по умолчанию
CMD ["make", "all", "CONFIG_NAME=config"]
